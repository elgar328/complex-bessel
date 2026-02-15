//! J Bessel function upper interface.
//!
//! Translation of Fortran ZBESJ from TOMS 644 (zbsubs.f lines 626-893).
//! J(fnu, z) = exp(fnu*pi*i/2) * I(fnu, -i*z) for Im(z) >= 0
//! J(fnu, z) = exp(-fnu*pi*i/2) * I(fnu, i*z) for Im(z) < 0

#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]

use num_complex::Complex;

use crate::algo::binu::zbinu;
use crate::algo::constants::HPI;
use crate::machine::BesselFloat;
use crate::types::{BesselError, BesselResult, Scaling};
use crate::utils::zabs;

/// Compute J_{fnu+j}(z) for j = 0, 1, ..., n-1.
///
/// Equivalent to Fortran ZBESJ in TOMS 644.
pub(crate) fn zbesj<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    scaling: Scaling,
    n: usize,
) -> Result<BesselResult<T>, BesselError> {
    let zero = T::zero();
    let one = T::one();
    let hpi_t = T::from(HPI).unwrap();

    // Input validation (Fortran lines 783-788)
    if n < 1 {
        return Err(BesselError::InvalidInput);
    }
    if fnu < zero {
        return Err(BesselError::InvalidInput);
    }

    // Machine constants (Fortran lines 800-812)
    let tol = T::tol();
    let elim = T::elim();
    let alim = T::alim();
    let rl = T::rl();
    let fnul = T::fnul();

    // Range check (Fortran lines 816-825)
    let az = zabs(z);
    let fn_val = fnu + T::from((n - 1) as f64).unwrap();
    let aa_tol = T::from(0.5).unwrap() / tol;
    let bb = T::from(2147483647.0 * 0.5).unwrap();
    let aa = aa_tol.min(bb);

    if az > aa || fn_val > aa {
        return Err(BesselError::TotalPrecisionLoss);
    }

    let aa_sqrt = aa.sqrt();
    let mut _precision_warning = false;
    if az > aa_sqrt || fn_val > aa_sqrt {
        _precision_warning = true;
    }

    // CSGN = exp(fnu*HPI*i) with precision preservation (Fortran lines 830-839)
    let mut cii = one;
    let inu = fnu.to_i32().unwrap();
    let inuh = inu / 2;
    let ir = inu - 2 * inuh;
    let arg = (fnu - T::from((inu - ir) as f64).unwrap()) * hpi_t;
    let mut csgnr = arg.cos();
    let mut csgni = arg.sin();
    if inuh % 2 != 0 {
        csgnr = -csgnr;
        csgni = -csgni;
    }

    // ZN is in the right half-plane (Fortran lines 844-850)
    // J(fnu, z) = exp(fnu*pi*i/2) * I(fnu, -i*z) for Im(z) >= 0
    let mut znr = z.im;
    let mut zni = -z.re;
    if z.im < zero {
        znr = -znr;
        zni = -zni;
        csgni = -csgni;
        cii = -cii;
    }
    let zn = Complex::new(znr, zni);

    // Call ZBINU (Fortran lines 852-853)
    let (mut cy, nz_raw) = zbinu(zn, fnu, scaling, n, rl, fnul, tol, elim, alim)?;
    let nz = nz_raw;

    // Apply phase factor (Fortran lines 855-878)
    let nl = n - nz;
    if nl == 0 {
        return Ok(BesselResult {
            values: cy,
            underflow_count: nz,
        });
    }

    let rtol = one / tol;
    let ascle = T::MACH_TINY * rtol * T::from(1.0e3).unwrap();

    for cy_item in cy.iter_mut().take(nl) {
        let mut aa_val = cy_item.re;
        let mut bb_val = cy_item.im;
        let mut atol = one;
        if aa_val.abs().max(bb_val.abs()) <= ascle {
            aa_val = aa_val * rtol;
            bb_val = bb_val * rtol;
            atol = tol;
        }
        let str = aa_val * csgnr - bb_val * csgni;
        let sti = aa_val * csgni + bb_val * csgnr;
        *cy_item = Complex::new(str * atol, sti * atol);

        // Advance CSGN: multiply by (0, cii) → CSGN *= i*CII
        // (Fortran lines 875-877)
        let str_new = -csgni * cii;
        csgni = csgnr * cii;
        csgnr = str_new;
    }

    Ok(BesselResult {
        values: cy,
        underflow_count: nz,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn besj_input_validation() {
        let z = Complex64::new(1.0, 0.0);
        assert!(zbesj(z, -1.0, Scaling::Unscaled, 1).is_err());
        assert!(zbesj(z, 0.0, Scaling::Unscaled, 0).is_err());
    }
}
