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
/// `nz` where nz = -1 indicates overflow.
#[inline]
pub(crate) fn zbunk<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    mr: i32,
    y: &mut [Complex<T>],
    tol: T,
    elim: T,
    alim: T,
) -> i32 {
    let ax = z.re.abs() * T::from_f64(1.7321);
    let ay = z.im.abs();

    if ay > ax {
        // Region 2: |Im(z)| > |Re(z)|*sqrt(3)
        zunk2(z, fnu, kode, mr, y, tol, elim, alim)
    } else {
        // Region 1: |arg(z)| <= pi/3
        zunk1(z, fnu, kode, mr, y, tol, elim, alim)
    }
}
