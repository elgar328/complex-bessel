//! K function dispatcher for large orders via uniform asymptotics.
//!
//! Translation of Fortran ZBUNK from TOMS 644 / SLATEC (zbsubs.f lines 3287-3321).
//! Dispatches to ZUNK1 (region 1: |arg(z)| <= pi/3) or ZUNK2 (region 2).

#![allow(clippy::too_many_arguments)]

use num_complex::Complex;

use crate::algo::unk1::zunk1;
use crate::algo::unk2::zunk2;
use crate::machine::BesselFloat;
use crate::types::Scaling;

/// Dispatch K function computation to ZUNK1 or ZUNK2 based on argument region.
///
/// Equivalent to Fortran ZBUNK in TOMS 644 (zbsubs.f lines 3287-3321).
///
/// # Returns
/// `(y, nz)` where nz = -1 indicates overflow.
pub(crate) fn zbunk<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    mr: i32,
    n: usize,
    tol: T,
    elim: T,
    alim: T,
) -> (Vec<Complex<T>>, i32) {
    let ax = z.re.abs() * T::from(1.7321).unwrap();
    let ay = z.im.abs();

    if ay > ax {
        // Region 2: |Im(z)| > |Re(z)|*sqrt(3)
        zunk2(z, fnu, kode, mr, n, tol, elim, alim)
    } else {
        // Region 1: |arg(z)| <= pi/3
        zunk1(z, fnu, kode, mr, n, tol, elim, alim)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    const TOL: f64 = 2.220446049250313e-16;
    const ELIM: f64 = 700.9217936944459;
    const ALIM: f64 = 664.8716455337102;

    #[test]
    fn zbunk_region1_dispatch() {
        // |Im(z)| < |Re(z)|*sqrt(3) → ZUNK1
        let z = Complex64::new(5.0, 1.0);
        let (y, nz) = zbunk(z, 90.0, Scaling::Unscaled, 0, 1, TOL, ELIM, ALIM);
        assert!(nz >= 0);
        assert!(y[0].re.is_finite());
    }

    #[test]
    fn zbunk_region2_dispatch() {
        // |Im(z)| > |Re(z)|*sqrt(3) → ZUNK2
        let z = Complex64::new(1.0, 10.0);
        let (y, nz) = zbunk(z, 90.0, Scaling::Unscaled, 0, 1, TOL, ELIM, ALIM);
        assert!(nz >= 0);
        assert!(y[0].re.is_finite());
    }

    #[test]
    fn zbunk_analytic_continuation() {
        let z = Complex64::new(-3.0, 5.0);
        let (y, nz) = zbunk(z, 90.0, Scaling::Unscaled, 1, 1, TOL, ELIM, ALIM);
        assert!(nz >= 0, "nz = {}", nz);
        assert!(y[0].re.is_finite());
    }
}
