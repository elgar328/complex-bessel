//! I Bessel function upper interface.
//!
//! Translation of Fortran ZBESI from TOMS 644 (zbsubs.f lines 355-625).

use num_complex::Complex;

use crate::algo::binu::zbinu;
use crate::algo::constants::PI;
use crate::machine::BesselFloat;
use crate::types::{Accuracy, Error, Scaling};
use crate::utils::zabs;

/// Compute I_{fnu+j}(z) for j = 0, 1, ..., n-1.
///
/// Equivalent to Fortran ZBESI in TOMS 644.
#[inline]
pub(crate) fn zbesi<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    scaling: Scaling,
    y: &mut [Complex<T>],
) -> Result<(usize, Accuracy), Error> {
    let zero = T::zero();
    let one = T::one();
    let pi_t = T::from_f64(PI);

    let n = y.len();

    // Input validation (Fortran IERR=1, lines 518-523)
    if n < 1 {
        return Err(Error::InvalidInput);
    }
    if fnu < zero {
        return Err(Error::InvalidInput);
    }

    // Machine constants (Fortran lines 535-547)
    let tol = T::tol();
    let elim = T::elim();
    let alim = T::alim();
    let rl = T::rl();
    let fnul = T::fnul();

    // Range check (Fortran lines 551-560)
    let az = zabs(z);
    let fn_val = fnu + T::from_f64((n - 1) as f64);
    let aa_tol = T::from_f64(0.5) / tol;
    let bb = T::from_f64(2147483647.0 * 0.5);
    let aa = aa_tol.min(bb);

    if az > aa || fn_val > aa {
        return Err(Error::TotalPrecisionLoss);
    }

    let aa_sqrt = aa.sqrt();
    let precision_warning = az > aa_sqrt || fn_val > aa_sqrt;

    // Compute in right half-plane, continue to left if needed
    // (Fortran lines 561-610)
    let (zn, mut csgn) = if z.re < zero {
        // Left half-plane: use I(fnu, -z) * exp(fnu*pi*i)
        // CSGN = exp(fnu*pi*i) with precision preservation (Fortran lines 572-579)
        // Safety: fnu is finite and < ~1e15 per upper-interface checks
        let inu = fnu.to_i32().unwrap();
        let mut arg = (fnu - T::from_f64(inu as f64)) * pi_t;
        if z.im < zero {
            arg = -arg;
        }
        let mut c = Complex::new(arg.cos(), arg.sin());
        if inu % 2 != 0 {
            c = -c;
        }
        (-z, c)
    } else {
        (z, Complex::from(one))
    };

    // Call ZBINU (Fortran lines 581-583)
    let nz = zbinu(zn, fnu, scaling, y, rl, fnul, tol, elim, alim)?;

    let status = if precision_warning {
        Accuracy::Reduced
    } else {
        Accuracy::Normal
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
    let ascle = T::MACH_TINY * rtol * T::from_f64(1.0e3);

    for cy_item in y.iter_mut().take(nn) {
        let (scaled, atol) = if cy_item.re.abs().max(cy_item.im.abs()) <= ascle {
            (*cy_item * rtol, tol)
        } else {
            (*cy_item, one)
        };
        *cy_item = scaled * csgn * atol;
        // CSGN alternates sign each order (Fortran lines 608-609)
        csgn = -csgn;
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
