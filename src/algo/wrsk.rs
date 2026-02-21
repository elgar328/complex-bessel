//! I Bessel function via Wronskian normalization.
//!
//! Translation of Fortran ZWRSK from TOMS 644 (zbsubs.f lines 3527-3621).
//! Computes I(fnu+k, z) for Re(z) >= 0 by normalizing the I-function
//! ratios from ZRATI using the Wronskian with K(fnu, z) and K(fnu+1, z).

use num_complex::Complex;

use crate::algo::bknu::zbknu;
use crate::algo::rati::zrati;
use crate::machine::BesselFloat;
use crate::types::{Error, Scaling};
use crate::utils::{zabs, zdiv};

/// Compute I Bessel function for Re(z) >= 0 via Wronskian normalization.
///
/// Writes results into `y` and returns `nz` (underflow count, always 0 on success).
///
/// Equivalent to Fortran ZWRSK in TOMS 644 (zbsubs.f lines 3527-3621).
///
/// # Algorithm
/// 1. Get K(fnu, z), K(fnu+1, z) from zbknu
/// 2. Get I-function ratios from zrati (written directly into y)
/// 3. Normalize via Wronskian: I(v)·K(v+1) + I(v+1)·K(v) = 1/z
/// 4. Forward recurrence for remaining values
pub(crate) fn zwrsk<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    y: &mut [Complex<T>],
    tol: T,
    elim: T,
    alim: T,
) -> Result<usize, Error> {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);
    let nz: usize = 0;
    let n = y.len();

    // ── Step 1: K(fnu, z) and K(fnu+1, z) via zbknu (Fortran line 3550) ──
    let mut cw = [czero; 2];
    let nw = zbknu(z, fnu, kode, &mut cw, tol, elim, alim)?;
    if nw != 0 {
        // Any nonzero NW (including underflows) is fatal for normalization
        return Err(Error::Overflow);
    }

    // ── Step 2: I-function ratios via zrati (written directly into y) ──
    zrati(z, fnu, y, tol);

    // ── Step 3: Normalization constant (Fortran lines 3557-3605) ──
    // KODE=1: CINU = 1; KODE=2: CINU = exp(i·Im(z))
    let mut cinu = match kode {
        Scaling::Unscaled => Complex::from(one),
        Scaling::Exponential => Complex::new(z.im.cos(), z.im.sin()),
    };

    // 3-level scaling to prevent under/overflow (Fortran lines 3569-3579)
    let acw = zabs(cw[1]);
    let ascle = T::from_f64(1.0e3) * T::MACH_TINY / tol;
    let csclr;
    if acw <= ascle {
        csclr = one / tol; // scale up: K is very small
    } else {
        let ascle_inv = one / ascle;
        if acw >= ascle_inv {
            csclr = tol; // scale down: K is very large
        } else {
            csclr = one; // no scaling needed
        }
    }

    // Scale K values (Fortran lines 3580-3583)
    let c1 = cw[0] * csclr;
    let c2 = cw[1] * csclr;

    // Save first ratio before overwriting (Fortran lines 3584-3585)
    let mut ratio_saved = y[0];

    // PT = ratio * C1 + C2 = (ratio * K(fnu) + K(fnu+1)) * csclr
    // (Fortran lines 3590-3593)
    let pt = ratio_saved * c1 + c2;

    // CT = z * PT (Fortran lines 3594-3595)
    let ct = z * pt;

    // CINU = CINU / CT (overflow-safe) (Fortran lines 3596-3603)
    cinu = zdiv(cinu, ct);

    // Y(1) = CINU * csclr = I(fnu, z) (Fortran lines 3604-3605)
    y[0] = cinu * csclr;

    if n == 1 {
        return Ok(nz);
    }

    // ── Step 4: Forward recurrence (Fortran lines 3607-3615) ──
    for item in y.iter_mut().skip(1) {
        cinu = ratio_saved * cinu;
        ratio_saved = *item; // save next ratio before overwriting
        *item = cinu * csclr;
    }

    Ok(nz)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn wrsk_wronskian_identity() {
        // Verify: I(v,z)*K(v+1,z) + I(v+1,z)*K(v,z) = 1/z
        let z = Complex64::new(2.0, 1.0);
        let fnu = 0.5;
        let tol = f64::tol();
        let elim = f64::elim();
        let alim = f64::alim();
        let czero = Complex64::new(0.0, 0.0);

        let mut i_vals = [czero; 2];
        let mut k_vals = [czero; 2];
        zwrsk(z, fnu, Scaling::Unscaled, &mut i_vals, tol, elim, alim).unwrap();
        zbknu(z, fnu, Scaling::Unscaled, &mut k_vals, tol, elim, alim).unwrap();

        let wronskian = i_vals[0] * k_vals[1] + i_vals[1] * k_vals[0];
        let expected = Complex64::new(1.0, 0.0) / z;
        let err = (wronskian - expected).norm() / expected.norm();
        assert!(err < 1e-13, "Wronskian error = {err:.2e}");
    }

    #[test]
    fn wrsk_i0_real_positive() {
        // I(0, 1.0) ≈ 1.2660658777520084
        let z = Complex64::new(1.0, 0.0);
        let tol = f64::tol();
        let elim = f64::elim();
        let alim = f64::alim();

        let mut y = [Complex64::new(0.0, 0.0)];
        let nz = zwrsk(z, 0.0, Scaling::Unscaled, &mut y, tol, elim, alim).unwrap();
        assert_eq!(nz, 0);
        let err = (y[0].re - 1.2660658777520084).abs();
        assert!(err < 1e-13, "I(0,1) error = {err:.2e}");
        assert!(y[0].im.abs() < 1e-14);
    }

    #[test]
    fn wrsk_sequence() {
        // I(0,2), I(1,2), I(2,2)
        let z = Complex64::new(2.0, 0.0);
        let tol = f64::tol();
        let elim = f64::elim();
        let alim = f64::alim();

        let mut y = [Complex64::new(0.0, 0.0); 3];
        zwrsk(z, 0.0, Scaling::Unscaled, &mut y, tol, elim, alim).unwrap();
        // I(0,2) ≈ 2.2795853023360673
        let err0 = (y[0].re - 2.2795853023360673).abs();
        assert!(err0 < 1e-13, "I(0,2) error = {err0:.2e}");
        // I(1,2) ≈ 1.5906368546373291
        let err1 = (y[1].re - 1.590_636_854_637_329).abs();
        assert!(err1 < 1e-13, "I(1,2) error = {err1:.2e}");
        // I(2,2) ≈ 0.688_948_447_698_738_1
        let err2 = (y[2].re - 0.688_948_447_698_738_1).abs();
        assert!(err2 < 1e-13, "I(2,2) error = {err2:.2e}");
    }
}
