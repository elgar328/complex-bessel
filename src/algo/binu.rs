//! I function master dispatcher.
//!
//! Translation of Fortran ZBINU from TOMS 644 (zbsubs.f lines 4378-4488).
//! Dispatches to ZSERI, ZASYI, ZMLRI, ZUOIK, ZWRSK based on |z| and fnu.

#![allow(clippy::too_many_arguments)]

use num_complex::Complex;

use crate::algo::asyi::zasyi;
use crate::algo::buni::zbuni;
use crate::algo::mlri::zmlri;
use crate::algo::seri::zseri;
use crate::algo::uoik::zuoik;
use crate::algo::wrsk::zwrsk;
use crate::machine::BesselFloat;
use crate::types::{Error, IkFlag, Scaling};
use crate::utils::zabs;

/// Compute I Bessel function in the right half z-plane.
///
/// Dispatches to the appropriate algorithm based on |z|, fnu, and n.
/// Writes results into `cy` and returns nz (underflow count) on success.
pub(crate) fn zbinu<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    cy: &mut [Complex<T>],
    rl: T,
    fnul: T,
    tol: T,
    elim: T,
    alim: T,
) -> Result<usize, Error> {
    let zero = T::zero();
    let one = T::one();
    let two = T::from_f64(2.0);
    let czero = Complex::new(zero, zero);

    let n = cy.len();
    cy.fill(czero);

    let mut nz: usize = 0;
    let az = zabs(z);
    let mut nn = n;
    let mut dfnu = fnu + T::from_f64((n - 1) as f64);

    // Dispatch: power series first (Fortran lines 4398-4411)
    if az <= two || az * az * T::from_f64(0.25) <= dfnu + one {
        // Label 10: power series (ZSERI)
        let nw = zseri(z, fnu, kode, &mut cy[..nn], tol, elim, alim);
        let inw = nw.unsigned_abs() as usize;
        nz += inw;
        nn -= inw;
        if nn == 0 {
            return Ok(nz);
        }
        if nw >= 0 {
            // Normal return from ZSERI
            return Ok(nz);
        }
        // nw < 0: need to continue with remaining terms
        dfnu = fnu + T::from_f64((nn - 1) as f64);
        // Fall through to label 20
    }

    // Label 20 onwards: dispatch to asymptotic, Miller, or Wronskian
    dispatch_20(
        z, fnu, kode, &mut nn, &mut nz, az, &mut dfnu, rl, fnul, tol, elim, alim, cy,
    )
}

/// Dispatch logic from Fortran label 20 onwards.
#[inline]
fn dispatch_20<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    nn: &mut usize,
    nz: &mut usize,
    az: T,
    dfnu: &mut T,
    rl: T,
    fnul: T,
    tol: T,
    elim: T,
    alim: T,
    cy: &mut [Complex<T>],
) -> Result<usize, Error> {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);

    // Fortran lines 4412-4422: check for asymptotic expansion
    if az >= rl {
        if *dfnu <= one || az + az >= *dfnu * *dfnu {
            // Label 30: use asymptotic expansion
            let nw = zasyi(z, fnu, kode, &mut cy[..*nn], rl, tol, elim, alim);
            if nw < 0 {
                return handle_error(nw);
            }
            return Ok(*nz);
        }
    } else if *dfnu <= one {
        // Label 70: use Miller algorithm (series normalization)
        let nw = zmlri(z, fnu, kode, &mut cy[..*nn], tol);
        if nw < 0 {
            return handle_error(nw);
        }
        return Ok(*nz);
    }

    // Label 50: overflow/underflow test on I sequence for Miller algorithm
    // (Fortran lines 4429-4437)
    let nw = zuoik(z, fnu, kode, IkFlag::I, &mut cy[..*nn], tol, elim, alim);
    if nw < 0 {
        return handle_error(nw);
    }
    *nz += nw as usize;
    *nn -= nw as usize;
    if *nn == 0 {
        return Ok(*nz);
    }
    *dfnu = fnu + T::from_f64((*nn - 1) as f64);

    if *dfnu > fnul || az > fnul {
        // Label 110: increment fnu+nn-1 up to fnul, compute and recur backward
        // (Fortran lines 4469-4481)
        // Safety: fnul and dfnu are finite f64-representable values
        let nui_f = (fnul - *dfnu).to_f64().unwrap() as i32 + 1;
        let nui = nui_f.max(0) as usize;
        let result = zbuni(z, fnu, kode, &mut cy[..*nn], nui, fnul, tol, elim, alim);
        if result.nz < 0 {
            return handle_error(result.nz);
        }
        *nz += result.nz as usize;
        if result.nlast == 0 {
            return Ok(*nz);
        }
        // NLAST != 0: retry from label 60 with reduced nn (Fortran GO TO 60)
        *nn = result.nlast;
    }

    // Label 60: check az vs rl (Fortran lines 4438-4439)
    if az > rl {
        // Label 80: Miller algorithm normalized by Wronskian
        // Overflow test on K functions used in Wronskian (Fortran lines 4454-4456)
        let mut cw_buf = [czero; 2];
        let nw = zuoik(z, fnu, kode, IkFlag::K, &mut cw_buf, tol, elim, alim);
        if nw < 0 {
            // All values underflow to zero (Fortran lines 4457-4462)
            *nz = *nn;
            for cy_item in cy.iter_mut().take(*nn) {
                *cy_item = czero;
            }
            return Ok(*nz);
        }
        if nw > 0 {
            return handle_error(-1);
        }
        // Call ZWRSK (Fortran lines 4465-4467)
        zwrsk(z, fnu, kode, &mut cy[..*nn], tol, elim, alim)?;
        return Ok(*nz);
    }

    // Label 70: Miller algorithm normalized by series
    let nw = zmlri(z, fnu, kode, &mut cy[..*nn], tol);
    if nw < 0 {
        return handle_error(nw);
    }
    Ok(*nz)
}

/// Map Fortran NW error codes to Error.
#[inline]
fn handle_error<T>(nw: i32) -> Result<T, Error> {
    if nw == -2 {
        Err(Error::ConvergenceFailure)
    } else {
        Err(Error::Overflow)
    }
}
