//! J Bessel function upper interface.
//!
//! Translation of Fortran ZBESJ from TOMS 644 (zbsubs.f lines 626-893).
//! J(fnu, z) = exp(fnu*pi*i/2) * I(fnu, -i*z) for Im(z) >= 0
//! J(fnu, z) = exp(-fnu*pi*i/2) * I(fnu, i*z) for Im(z) < 0

use num_complex::Complex;

use crate::algo::binu::zbinu;
use crate::algo::constants::HPI;
use crate::machine::BesselFloat;
use crate::types::{Accuracy, Error, Scaling};
use crate::utils::{mul_i, mul_neg_i, zabs};

/// Compute J_{fnu+j}(z) for j = 0, 1, ..., n-1.
///
/// Equivalent to Fortran ZBESJ in TOMS 644.
#[inline]
pub(crate) fn zbesj<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    scaling: Scaling,
    y: &mut [Complex<T>],
) -> Result<(usize, Accuracy), Error> {
    let zero = T::zero();
    let one = T::one();
    let hpi_t = T::from_f64(HPI);

    let n = y.len();

    // Input validation (Fortran lines 783-788)
    if n < 1 {
        return Err(Error::InvalidInput);
    }
    if fnu < zero {
        return Err(Error::InvalidInput);
    }

    // Machine constants (Fortran lines 800-812)
    let tol = T::tol();
    let elim = T::elim();
    let alim = T::alim();
    let rl = T::rl();
    let fnul = T::fnul();

    // Range check (Fortran lines 816-825)
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

    // CSGN = exp(fnu*HPI*i) with precision preservation (Fortran lines 830-839)
    let mut cii = one;
    // Safety: fnu is finite and < ~1e15 per upper-interface checks
    let inu = fnu.to_i32().unwrap();
    let inuh = inu / 2;
    let ir = inu - 2 * inuh;
    let arg = (fnu - T::from_f64((inu - ir) as f64)) * hpi_t;
    let mut csgn = Complex::new(arg.cos(), arg.sin());
    if inuh % 2 != 0 {
        csgn = -csgn;
    }

    // ZN is in the right half-plane (Fortran lines 844-850)
    // J(fnu, z) = exp(fnu*pi*i/2) * I(fnu, -i*z) for Im(z) >= 0
    let mut zn = Complex::new(z.im, -z.re);
    if z.im < zero {
        zn = -zn;
        csgn = csgn.conj();
        cii = -cii;
    }

    // Call ZBINU (Fortran lines 852-853)
    let nz = zbinu(zn, fnu, scaling, y, rl, fnul, tol, elim, alim)?;

    let status = if precision_warning {
        Accuracy::Reduced
    } else {
        Accuracy::Normal
    };

    // Apply phase factor (Fortran lines 855-878)
    let nl = n - nz;
    if nl == 0 {
        return Ok((nz, status));
    }

    let rtol = one / tol;
    let ascle = T::MACH_TINY * rtol * T::from_f64(1.0e3);

    for cy_item in y.iter_mut().take(nl) {
        let (scaled, atol) = if cy_item.re.abs().max(cy_item.im.abs()) <= ascle {
            (*cy_item * rtol, tol)
        } else {
            (*cy_item, one)
        };
        *cy_item = scaled * csgn * atol;

        // Advance CSGN: multiply by (0, cii) → CSGN *= i*CII
        // (Fortran lines 875-877)
        csgn = if cii > zero {
            mul_i(csgn)
        } else {
            mul_neg_i(csgn)
        };
    }

    Ok((nz, status))
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn besj_input_validation() {
        let z = Complex64::new(1.0, 0.0);
        assert!(zbesj(z, -1.0, Scaling::Unscaled, &mut [Complex64::new(0.0, 0.0)]).is_err());
        assert!(zbesj(z, 0.0, Scaling::Unscaled, &mut []).is_err());
    }
}
