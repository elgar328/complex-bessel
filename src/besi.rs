//! I Bessel function upper interface.
//!
//! Translation of Fortran ZBESI from TOMS 644 (zbsubs.f lines 355-624).

#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]

use num_complex::Complex;

use crate::algo::binu::zbinu;
use crate::algo::constants::PI;
use crate::machine::BesselFloat;
use crate::types::{BesselError, BesselStatus, Scaling};
use crate::utils::zabs;

/// Compute I_{fnu+j}(z) for j = 0, 1, ..., n-1.
///
/// Equivalent to Fortran ZBESI in TOMS 644.
pub(crate) fn zbesi<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    scaling: Scaling,
    y: &mut [Complex<T>],
) -> Result<(usize, BesselStatus), BesselError> {
    let zero = T::zero();
    let one = T::one();
    let pi_t = T::from(PI).unwrap();
    let czero = Complex::new(zero, zero);

    let n = y.len();

    // Zero the output buffer
    for v in y.iter_mut() {
        *v = czero;
    }

    // Input validation (Fortran IERR=1, lines 518-523)
    if n < 1 {
        return Err(BesselError::InvalidInput);
    }
    if fnu < zero {
        return Err(BesselError::InvalidInput);
    }

    // Machine constants (Fortran lines 535-547)
    let tol = T::tol();
    let elim = T::elim();
    let alim = T::alim();
    let rl = T::rl();
    let fnul = T::fnul();

    // Range check (Fortran lines 551-560)
    let az = zabs(z);
    let fn_val = fnu + T::from((n - 1) as f64).unwrap();
    let aa_tol = T::from(0.5).unwrap() / tol;
    let bb = T::from(2147483647.0 * 0.5).unwrap();
    let aa = aa_tol.min(bb);

    if az > aa || fn_val > aa {
        return Err(BesselError::TotalPrecisionLoss);
    }

    let aa_sqrt = aa.sqrt();
    let mut precision_warning = false;
    if az > aa_sqrt || fn_val > aa_sqrt {
        precision_warning = true;
    }

    // Compute in right half-plane, continue to left if needed
    // (Fortran lines 561-610)
    let mut znr = z.re;
    let mut zni = z.im;
    let mut csgnr = one;
    let mut csgni = zero;

    if z.re < zero {
        // Left half-plane: use I(fnu, -z) * exp(fnu*pi*i)
        znr = -z.re;
        zni = -z.im;

        // CSGN = exp(fnu*pi*i) with precision preservation (Fortran lines 572-579)
        let inu = fnu.to_i32().unwrap();
        let mut arg = (fnu - T::from(inu as f64).unwrap()) * pi_t;
        if z.im < zero {
            arg = -arg;
        }
        csgnr = arg.cos();
        csgni = arg.sin();
        if inu % 2 != 0 {
            csgnr = -csgnr;
            csgni = -csgni;
        }
    }

    let zn = Complex::new(znr, zni);

    // Call ZBINU (Fortran lines 581-583)
    let nz = zbinu(zn, fnu, scaling, y, rl, fnul, tol, elim, alim)?;

    let status = if precision_warning {
        BesselStatus::ReducedPrecision
    } else {
        BesselStatus::Normal
    };

    if z.re >= zero {
        // Right half-plane: done
        return Ok((nz, status));
    }

    // Analytic continuation to left half-plane (Fortran lines 586-610)
    let nn = n - nz;
    if nn == 0 {
        return Ok((nz, status));
    }

    let rtol = one / tol;
    let ascle = T::MACH_TINY * rtol * T::from(1.0e3).unwrap();

    for cy_item in y.iter_mut().take(nn) {
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
        // CSGN alternates sign each order (Fortran lines 608-609)
        csgnr = -csgnr;
        csgni = -csgni;
    }

    Ok((nz, status))
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn besi_input_validation() {
        let z = Complex64::new(1.0, 0.0);
        assert!(zbesi(z, -1.0, Scaling::Unscaled, &mut [Complex64::new(0.0, 0.0)]).is_err());
        assert!(zbesi(z, 0.0, Scaling::Unscaled, &mut []).is_err());
    }
}
