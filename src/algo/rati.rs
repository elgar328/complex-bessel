//! I Bessel function ratios by backward recurrence.
//!
//! Translation of Fortran ZRATI from TOMS 644 (zbsubs.f lines 3104-3236).
//! Uses the Sookne algorithm (J. Res. NBS-B, Vol 77B, 1973) to determine
//! the starting index for backward recurrence.

// Exact Fortran constants — preserve verbatim.
#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]

use num_complex::Complex;

use crate::machine::BesselFloat;
use crate::utils::{mul_add, reciprocal_z, zabs, zdiv};

/// Compute I-function ratios by backward recurrence.
///
/// Writes `cy[k]` = I(fnu+k+1, z) / I(fnu+k, z) for k = 0, 1, ..., n-1.
///
/// Equivalent to Fortran ZRATI in TOMS 644 (zbsubs.f lines 3104-3236).
///
/// # Parameters
/// - `z`: complex argument (must be nonzero)
/// - `fnu`: starting order (>= 0)
/// - `cy`: output slice of length n (>= 1)
/// - `tol`: convergence tolerance
pub(crate) fn zrati<T: BesselFloat>(z: Complex<T>, fnu: T, cy: &mut [Complex<T>], tol: T) {
    let zero = T::zero();
    let one = T::one();
    // sqrt(2), Fortran DATA constant (zbsubs.f line 3126)
    let rt2 = T::from_f64(1.41421356237309505);

    let n = cy.len();

    // ── Initialization (Fortran lines 3127-3148) ──
    let az = zabs(z);
    // Safety: fnu is finite and < ~1e15 per upper-interface checks
    let inu = fnu.to_i32().unwrap();
    let idnu = inu + n as i32 - 1;
    // Safety: az (modulus) is finite and bounded
    let magz = az.to_i32().unwrap();
    let amagz = T::from_f64((magz + 1) as f64);
    let fdnu = T::from_f64(idnu as f64);
    let fnup = amagz.max(fdnu);
    let mut id = idnu - magz - 1;
    if id > 0 {
        id = 0;
    }
    let mut refined = false;
    let mut k: i32 = 1;

    // RZ = 2/z (overflow-safe, Fortran lines 3137-3139)
    let rz = reciprocal_z(z);

    // T1 = RZ * FNUP (Fortran lines 3140-3141)
    let mut t1 = rz * fnup;
    // P2 = -T1 (Fortran lines 3142-3143)
    let mut p2 = -t1;
    // P1 = CONE (Fortran lines 3144-3145)
    let mut p1 = Complex::from(one);
    // T1 += RZ (Fortran lines 3146-3147)
    t1 = t1 + rz;

    let mut ap2 = zabs(p2);
    let mut ap1 = zabs(p1);

    // Scale test to prevent premature overflow (Fortran lines 3157-3165)
    let arg = (ap2 + ap2) / (ap1 * tol);
    let test1 = arg.sqrt();
    let mut test = test1;
    let rap1 = one / ap1;
    p1 = p1 * rap1;
    p2 = p2 * rap1;
    ap2 = ap2 * rap1;

    // ── Stage 1: Forward recurrence for starting index (Fortran label 10) ──
    loop {
        k += 1;
        ap1 = ap2;
        let pt = p2;
        p2 = mul_add(-t1, pt, p1);
        p1 = pt;
        t1 = t1 + rz;
        ap2 = zabs(p2);

        if ap1 <= test {
            continue;
        }
        if refined {
            break;
        }
        // First pass: refine convergence test (Fortran lines 3180-3184)
        let ak = zabs(t1) * T::from_f64(0.5);
        let flam = ak + (ak * ak - one).sqrt();
        let rho = (ap2 / ap1).min(flam);
        test = test1 * (rho / (rho * rho - one)).sqrt();
        refined = true;
    }

    // ── Stage 2: Backward recurrence (Fortran label 20, lines 3186-3212) ──
    let kk = (k + 1 - id) as usize;
    let dfnu = fnu + T::from_f64((n - 1) as f64);
    let mut p1 = Complex::from(one / ap2);
    let mut p2 = Complex::new(zero, zero);
    let mut t1r_bk = T::from_f64(kk as f64);

    for _ in 0..kk {
        let pt = p1;
        let rap1 = dfnu + t1r_bk;
        let tt = rz * rap1;
        p1 = mul_add(pt, tt, p2);
        p2 = pt;
        t1r_bk = t1r_bk - one;
    }

    // Protect against zero denominator (Fortran lines 3208-3210)
    if p1.re == zero && p1.im == zero {
        p1 = Complex::new(tol, tol);
    }

    // CY(N) = P2 / P1 (Fortran line 3212)
    cy[n - 1] = zdiv(p2, p1);

    if n == 1 {
        return;
    }

    // ── Stage 3: Forward ratio computation (Fortran lines 3214-3234) ──
    let cdfnu = rz * fnu;

    for k_idx in (0..=(n - 2)).rev() {
        // T1R = k_idx + 1 corresponds to Fortran K (1-based)
        let t1r_val = T::from_f64((k_idx + 1) as f64);
        let mut pt = cdfnu + rz * t1r_val + cy[k_idx + 1];
        let mut ak = zabs(pt);

        if ak == zero {
            pt = Complex::new(tol, tol);
            ak = tol * rt2;
        }
        // CY(K) = conj(PT)/|PT|^2 = 1/PT (Fortran lines 3229-3231)
        cy[k_idx] = zdiv(Complex::from(one), pt);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn rati_self_consistency() {
        // Verify the recurrence relation: 1/cy[k] = (fnu+k+1)*RZ + cy[k+1]
        // where RZ = 2/z
        let z = Complex64::new(2.0, 1.0);
        let fnu = 0.5;
        let n = 5;
        let tol = f64::tol();
        let mut cy = [Complex64::new(0.0, 0.0); 5];
        zrati(z, fnu, &mut cy, tol);

        let rz = Complex64::new(2.0, 0.0) / z;
        for k in 0..(n - 1) {
            let v = fnu + (k + 1) as f64; // order index for the recurrence
            let lhs = Complex64::new(1.0, 0.0) / cy[k];
            let rhs = Complex64::from(v) * rz + cy[k + 1];
            let err = (lhs - rhs).norm() / lhs.norm();
            assert!(err < 1e-12, "Recurrence failed at k={k}: err={err:.2e}");
        }
    }

    #[test]
    fn rati_n1_returns_single_ratio() {
        let z = Complex64::new(3.0, 0.0);
        let fnu = 1.0;
        let tol = f64::tol();
        let mut cy = [Complex64::new(0.0, 0.0); 1];
        zrati(z, fnu, &mut cy, tol);
        // cy[0] = I(2,3)/I(1,3)
        // I(1,3) ≈ 3.95337, I(2,3) ≈ 2.24521 → ratio ≈ 0.5680
        let ratio = cy[0].re;
        assert!(
            (ratio - 0.568).abs() < 0.01,
            "I(2,3)/I(1,3) ≈ {ratio}, expected ~0.568"
        );
        assert!(cy[0].im.abs() < 1e-14);
    }

    #[test]
    fn rati_real_positive_z() {
        // For real positive z, ratios should be real and positive
        let z = Complex64::new(5.0, 0.0);
        let fnu = 0.0;
        let n = 4;
        let tol = f64::tol();
        let mut cy = [Complex64::new(0.0, 0.0); 4];
        zrati(z, fnu, &mut cy, tol);

        for (k, c) in cy.iter().enumerate().take(n) {
            assert!(c.re > 0.0, "cy[{k}].re = {} should be positive", c.re);
            assert!(c.im.abs() < 1e-14, "cy[{k}].im = {} should be ~0", c.im);
        }
    }
}
