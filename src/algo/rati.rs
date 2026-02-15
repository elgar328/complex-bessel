//! I Bessel function ratios by backward recurrence.
//!
//! Translation of Fortran ZRATI from TOMS 644.
//! Uses the Sookne algorithm (J. Res. NBS-B, Vol 77B, 1973) to determine
//! the starting index for backward recurrence.

// Exact Fortran constants — preserve verbatim.
#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]

#[cfg(not(feature = "std"))]
use alloc::vec;
#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

use num_complex::Complex;

use crate::machine::BesselFloat;
use crate::utils::{zabs, zdiv};

/// Compute I-function ratios by backward recurrence.
///
/// Returns `cy[k]` = I(fnu+k+1, z) / I(fnu+k, z) for k = 0, 1, ..., n-1.
///
/// Equivalent to Fortran ZRATI in TOMS 644 (zbsubs.f lines 3104-3236).
///
/// # Parameters
/// - `z`: complex argument (must be nonzero)
/// - `fnu`: starting order (>= 0)
/// - `n`: number of ratios to compute (>= 1)
/// - `tol`: convergence tolerance
pub(crate) fn zrati<T: BesselFloat>(z: Complex<T>, fnu: T, n: usize, tol: T) -> Vec<Complex<T>> {
    let zero = T::zero();
    let one = T::one();
    // sqrt(2), Fortran DATA constant (zbsubs.f line 3126)
    let rt2 = T::from(1.41421356237309505).unwrap();

    // ── Initialization (Fortran lines 3127-3148) ──
    let az = zabs(z);
    let inu = fnu.to_i32().unwrap();
    let idnu = inu + n as i32 - 1;
    let magz = az.to_i32().unwrap();
    let amagz = T::from((magz + 1) as f64).unwrap();
    let fdnu = T::from(idnu as f64).unwrap();
    let fnup = amagz.max(fdnu);
    let mut id = idnu - magz - 1;
    if id > 0 {
        id = 0;
    }
    let mut itime: i32 = 1;
    let mut k: i32 = 1;

    // RZ = 2/z (overflow-safe, Fortran lines 3137-3139)
    let ptr = one / az;
    let rzr = ptr * (z.re + z.re) * ptr;
    let rzi = -ptr * (z.im + z.im) * ptr;

    // T1 = RZ * FNUP (Fortran lines 3140-3141)
    let mut t1r = rzr * fnup;
    let mut t1i = rzi * fnup;
    // P2 = -T1 (Fortran lines 3142-3143)
    let mut p2r = -t1r;
    let mut p2i = -t1i;
    // P1 = CONE (Fortran lines 3144-3145)
    let mut p1r = one;
    let mut p1i = zero;
    // T1 += RZ (Fortran lines 3146-3147)
    t1r = t1r + rzr;
    t1i = t1i + rzi;

    let mut ap2 = zabs(Complex::new(p2r, p2i));
    let mut ap1 = zabs(Complex::new(p1r, p1i));

    // Scale test to prevent premature overflow (Fortran lines 3157-3165)
    let arg = (ap2 + ap2) / (ap1 * tol);
    let test1 = arg.sqrt();
    let mut test = test1;
    let rap1 = one / ap1;
    p1r = p1r * rap1;
    p1i = p1i * rap1;
    p2r = p2r * rap1;
    p2i = p2i * rap1;
    ap2 = ap2 * rap1;

    // ── Stage 1: Forward recurrence for starting index (Fortran label 10) ──
    loop {
        k += 1;
        ap1 = ap2;
        let pt_r = p2r;
        let pt_i = p2i;
        p2r = p1r - (t1r * pt_r - t1i * pt_i);
        p2i = p1i - (t1r * pt_i + t1i * pt_r);
        p1r = pt_r;
        p1i = pt_i;
        t1r = t1r + rzr;
        t1i = t1i + rzi;
        ap2 = zabs(Complex::new(p2r, p2i));

        if ap1 <= test {
            continue;
        }
        if itime == 2 {
            break;
        }
        // First pass: refine convergence test (Fortran lines 3180-3184)
        let ak = zabs(Complex::new(t1r, t1i)) * T::from(0.5).unwrap();
        let flam = ak + (ak * ak - one).sqrt();
        let rho = (ap2 / ap1).min(flam);
        test = test1 * (rho / (rho * rho - one)).sqrt();
        itime = 2;
    }

    // ── Stage 2: Backward recurrence (Fortran label 20, lines 3186-3212) ──
    let kk = (k + 1 - id) as usize;
    let dfnu = fnu + T::from((n - 1) as f64).unwrap();
    let mut p1r = one / ap2;
    let mut p1i = zero;
    let mut p2r = zero;
    let mut p2i = zero;
    let mut t1r_bk = T::from(kk as f64).unwrap();

    for _ in 0..kk {
        let pt_r = p1r;
        let pt_i = p1i;
        let rap1 = dfnu + t1r_bk;
        let ttr = rzr * rap1;
        let tti = rzi * rap1;
        p1r = (pt_r * ttr - pt_i * tti) + p2r;
        p1i = (pt_r * tti + pt_i * ttr) + p2i;
        p2r = pt_r;
        p2i = pt_i;
        t1r_bk = t1r_bk - one;
    }

    // Protect against zero denominator (Fortran lines 3208-3210)
    if p1r == zero && p1i == zero {
        p1r = tol;
        p1i = tol;
    }

    // CY(N) = P2 / P1 (Fortran line 3212)
    let mut cy = vec![Complex::new(zero, zero); n];
    cy[n - 1] = zdiv(Complex::new(p2r, p2i), Complex::new(p1r, p1i));

    if n == 1 {
        return cy;
    }

    // ── Stage 3: Forward ratio computation (Fortran lines 3214-3234) ──
    let cdfnur = fnu * rzr;
    let cdfnui = fnu * rzi;

    for k_idx in (0..=(n - 2)).rev() {
        // T1R = k_idx + 1 corresponds to Fortran K (1-based)
        let t1r_val = T::from((k_idx + 1) as f64).unwrap();
        let mut pt_r = cdfnur + (t1r_val * rzr) + cy[k_idx + 1].re;
        let mut pt_i = cdfnui + (t1r_val * rzi) + cy[k_idx + 1].im;
        let mut ak = zabs(Complex::new(pt_r, pt_i));

        if ak == zero {
            pt_r = tol;
            pt_i = tol;
            ak = tol * rt2;
        }
        // CY(K) = conj(PT)/|PT|^2 = 1/PT (Fortran lines 3229-3231)
        let rak = one / ak;
        cy[k_idx] = Complex::new(rak * pt_r * rak, -rak * pt_i * rak);
    }

    cy
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
        let cy = zrati(z, fnu, n, tol);

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
        let cy = zrati(z, fnu, 1, tol);
        assert_eq!(cy.len(), 1);
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
        let cy = zrati(z, fnu, n, tol);

        for (k, c) in cy.iter().enumerate() {
            assert!(c.re > 0.0, "cy[{k}].re = {} should be positive", c.re);
            assert!(c.im.abs() < 1e-14, "cy[{k}].im = {} should be ~0", c.im);
        }
    }
}
