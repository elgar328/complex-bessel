//! I function master dispatcher.
//!
//! Translation of Fortran ZBINU from TOMS 644 (zbsubs.f lines 4378-4488).
//! Dispatches to ZSERI, ZASYI, ZMLRI, ZUOIK, ZWRSK based on |z| and fnu.

#![allow(clippy::too_many_arguments)]

use num_complex::Complex;
use num_traits::Float;

use crate::algo::asyi::zasyi;
use crate::algo::mlri::zmlri;
use crate::algo::seri::zseri;
use crate::algo::uoik::zuoik;
use crate::algo::wrsk::zwrsk;
use crate::machine::BesselFloat;
use crate::types::{BesselError, Scaling};
use crate::utils::zabs;

/// Compute I Bessel function in the right half z-plane.
///
/// Dispatches to the appropriate algorithm based on |z|, fnu, and n.
///
/// Returns (y, nz) on success, where nz is the underflow count.
pub(crate) fn zbinu<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    n: usize,
    rl: T,
    fnul: T,
    tol: T,
    elim: T,
    alim: T,
) -> Result<(Vec<Complex<T>>, usize), BesselError> {
    let zero = T::zero();
    let one = T::one();
    let two = T::from(2.0).unwrap();
    let czero = Complex::new(zero, zero);

    let mut nz: usize = 0;
    let az = zabs(z);
    let mut nn = n;
    let mut dfnu = fnu + T::from((n - 1) as f64).unwrap();

    // Dispatch: power series first (Fortran lines 4398-4411)
    if az <= two || az * az * T::from(0.25).unwrap() <= dfnu + one {
        // Label 10: power series (ZSERI)
        let (mut cy, nw) = zseri(z, fnu, kode, nn, tol, elim, alim);
        let inw = nw.unsigned_abs() as usize;
        nz += inw;
        nn -= inw;
        if nn == 0 {
            return Ok((cy, nz));
        }
        if nw >= 0 {
            // Normal return from ZSERI
            return Ok((cy, nz));
        }
        // nw < 0: need to continue with remaining terms
        dfnu = fnu + T::from((nn - 1) as f64).unwrap();
        // Fall through to label 20
        return dispatch_20(
            z, fnu, kode, nn, nz, az, dfnu, rl, fnul, tol, elim, alim, &mut cy,
        );
    }

    // Label 20: check asymptotic vs Miller
    let mut cy = vec![czero; n];
    dispatch_20(
        z, fnu, kode, nn, nz, az, dfnu, rl, fnul, tol, elim, alim, &mut cy,
    )
}

/// Dispatch logic from Fortran label 20 onwards.
fn dispatch_20<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    mut nn: usize,
    mut nz: usize,
    az: T,
    mut dfnu: T,
    rl: T,
    fnul: T,
    tol: T,
    elim: T,
    alim: T,
    cy: &mut [Complex<T>],
) -> Result<(Vec<Complex<T>>, usize), BesselError> {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);

    // Fortran lines 4412-4422: check for asymptotic expansion
    if az >= rl {
        if dfnu <= one {
            // Goto 30: use asymptotic expansion
            return call_zasyi(z, fnu, kode, nn, nz, rl, tol, elim, alim, cy);
        }
        if az + az >= dfnu * dfnu {
            // Goto 30: use asymptotic expansion
            return call_zasyi(z, fnu, kode, nn, nz, rl, tol, elim, alim, cy);
        }
    } else if dfnu <= one {
        // Goto 70: use Miller algorithm (series normalization)
        return call_zmlri(z, fnu, kode, nn, nz, tol, cy);
    }

    // Label 50: overflow/underflow test on I sequence for Miller algorithm
    // (Fortran lines 4429-4437)
    let (_uoik_y, nw) = zuoik(z, fnu, kode, 1, nn, tol, elim, alim);
    if nw < 0 {
        return handle_error(nw);
    }
    nz += nw as usize;
    nn -= nw as usize;
    if nn == 0 {
        return Ok((cy.to_vec(), nz));
    }
    dfnu = fnu + T::from((nn - 1) as f64).unwrap();

    if dfnu > fnul || az > fnul {
        // Label 110: ZBUNI — not implemented in Phase 4
        // (Fortran lines 4469-4481)
        return Err(BesselError::ConvergenceFailure);
    }

    // Label 60: check az vs rl (Fortran lines 4438-4439)
    if az > rl {
        // Label 80: Miller algorithm normalized by Wronskian
        return call_zwrsk(z, fnu, kode, nn, nz, tol, elim, alim, cy);
    }

    // Label 70: Miller algorithm normalized by series
    call_zmlri(z, fnu, kode, nn, nz, tol, cy)
}

/// Call ZASYI (asymptotic expansion for large z).
fn call_zasyi<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    nn: usize,
    nz: usize,
    rl: T,
    tol: T,
    elim: T,
    alim: T,
    cy: &mut [Complex<T>],
) -> Result<(Vec<Complex<T>>, usize), BesselError> {
    let (result, nw) = zasyi(z, fnu, kode, nn, rl, tol, elim, alim);
    if nw < 0 {
        return handle_error(nw);
    }
    // Copy results into cy
    cy[..nn].copy_from_slice(&result[..nn]);
    Ok((cy.to_vec(), nz))
}

/// Call ZMLRI (Miller algorithm, series normalization).
fn call_zmlri<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    nn: usize,
    nz: usize,
    tol: T,
    cy: &mut [Complex<T>],
) -> Result<(Vec<Complex<T>>, usize), BesselError> {
    let (result, nw) = zmlri(z, fnu, kode, nn, tol);
    if nw < 0 {
        return handle_error(nw);
    }
    cy[..nn].copy_from_slice(&result[..nn]);
    Ok((cy.to_vec(), nz))
}

/// Call ZWRSK (Miller algorithm, Wronskian normalization).
fn call_zwrsk<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    nn: usize,
    mut nz: usize,
    tol: T,
    elim: T,
    alim: T,
    cy: &mut [Complex<T>],
) -> Result<(Vec<Complex<T>>, usize), BesselError> {
    let zero = T::zero();
    let czero = Complex::new(zero, zero);

    // Overflow test on K functions used in Wronskian (Fortran lines 4454-4456)
    let (_cw_y, nw) = zuoik(z, fnu, kode, 2, 2, tol, elim, alim);
    if nw < 0 {
        // All values underflow to zero (Fortran lines 4457-4462)
        nz = nn;
        for cy_item in cy.iter_mut().take(nn) {
            *cy_item = czero;
        }
        return Ok((cy.to_vec(), nz));
    }
    if nw > 0 {
        return handle_error(-1);
    }

    // Call ZWRSK (Fortran lines 4465-4467)
    let (result, _nw_wrsk) = zwrsk(z, fnu, kode, nn, tol, elim, alim)?;
    cy[..nn].copy_from_slice(&result[..nn]);
    Ok((cy.to_vec(), nz))
}

/// Map Fortran NW error codes to BesselError.
fn handle_error<T>(nw: i32) -> Result<T, BesselError> {
    if nw == -2 {
        Err(BesselError::ConvergenceFailure)
    } else {
        Err(BesselError::Overflow)
    }
}
