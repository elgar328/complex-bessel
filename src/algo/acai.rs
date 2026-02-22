//! Airy-specific analytic continuation of K function to the left half z-plane.
//!
//! Translation of Fortran ZACAI from TOMS 644 (zbsubs.f lines 4678-4777).
//! Simplified version of ZACON for Airy functions where FNU is 1/3 or 2/3 and N=1.
//! Calls ZSERI/ZMLRI/ZASYI directly instead of ZBINU to avoid recursion.

#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]
#![allow(clippy::too_many_arguments)]

use num_complex::Complex;

use crate::algo::asyi::zasyi;
use crate::algo::bknu::zbknu;
use crate::algo::constants::PI;
use crate::algo::mlri::zmlri;
use crate::algo::s1s2::zs1s2;
use crate::algo::seri::zseri;
use crate::machine::BesselFloat;
use crate::types::{Error, Scaling};
use crate::utils::zabs;

/// Analytic continuation for Airy functions.
///
/// Applies: K(fnu, zn*exp(mp)) = K(fnu, zn)*exp(-mp*fnu) - mp*I(fnu, zn)
/// where mp = pi*mr*i.
///
/// Same formula as ZACON but with higher-order recurrence removed.
/// FNU is 1/3 or 2/3, N is always 1.
///
/// Returns nz where nz < 0 indicates error.
pub(crate) fn zacai<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    mr: i32,
    y: &mut [Complex<T>],
    rl: T,
    tol: T,
    elim: T,
    alim: T,
) -> Result<i32, Error> {
    let zero = T::zero();
    let one = T::one();
    let two = T::from_f64(2.0);
    let pi_t = T::from_f64(PI);
    let czero = Complex::new(zero, zero);

    let mut nz: i32 = 0;
    let n = y.len();
    y.fill(czero);

    // ZN = -Z (Fortran line 4703-4704)
    let zn = -z;
    let az = zabs(z);
    let dfnu = fnu + T::from_f64((n - 1) as f64);

    // I function dispatch (Fortran lines 4706-4731)
    // Direct calls to avoid ZBINU → ZBUNI → ZUNI2 → ZAIRY recursion
    // N is always 1 for acai, so use stack buffer
    let mut i_buf = [czero];
    if az <= two || az * az * T::from_f64(0.25) <= dfnu + one {
        // Label 10: power series (Fortran line 4710)
        zseri(zn, fnu, kode, &mut i_buf, tol, elim, alim);
    } else if az >= rl {
        // Label 20: asymptotic expansion for large z (Fortran line 4718)
        let nw = zasyi(zn, fnu, kode, &mut i_buf, rl, tol, elim, alim);
        if nw < 0 {
            return Err(if nw == -2 {
                Error::ConvergenceFailure
            } else {
                Error::Overflow
            });
        }
    } else {
        // Label 30: Miller algorithm (Fortran line 4726)
        let nw = zmlri(zn, fnu, kode, &mut i_buf, tol);
        if nw < 0 {
            return Err(if nw == -2 {
                Error::ConvergenceFailure
            } else {
                Error::Overflow
            });
        }
    }

    // Label 40: Analytic continuation (Fortran lines 4732-4777)
    // K function at -z (Fortran line 4736)
    let mut k_buf = [czero];
    let nw_k = zbknu(zn, fnu, kode, &mut k_buf, tol, elim, alim)?;
    if nw_k != 0 {
        // NW != 0 → error (Fortran line 4737 GO TO 80)
        return Err(Error::Overflow);
    }

    let fmr = T::from_f64(mr as f64);
    let sgn = -pi_t.copysign(fmr); // SGN = -DSIGN(PI, FMR) (Fortran line 4739)

    // CSGN = (0, sgn) (Fortran lines 4740-4741)
    let mut csgn = Complex::new(zero, sgn);

    if kode == Scaling::Exponential {
        // Fortran lines 4743-4745: KODE=2 adjustment
        let yy = -zn.im;
        csgn = csgn * Complex::new(yy.cos(), yy.sin());
    }

    // CSPN = exp(FNU*PI*I) with precision preservation (Fortran lines 4751-4758)
    // Safety: fnu is finite and < ~1e15 per upper-interface checks
    let inu = fnu.to_i32().unwrap();
    let arg = (fnu - T::from_f64(inu as f64)) * sgn;
    let mut cspn = Complex::new(arg.cos(), arg.sin());
    if inu % 2 != 0 {
        cspn = -cspn;
    }

    // C1 = K result, C2 = I result (Fortran lines 4760-4763)
    let mut c1 = k_buf[0];
    let mut c2 = i_buf[0];

    if kode == Scaling::Exponential {
        // Fortran lines 4765-4769: ZS1S2 call for KODE=2
        let iuf: i32 = 0;
        let ascle = T::from_f64(1.0e3) * T::MACH_TINY / tol;
        let s1s2_out = zs1s2(zn, c1, c2, ascle, alim, iuf);
        nz += s1s2_out.nz;
        c1 = s1s2_out.s1;
        c2 = s1s2_out.s2;
    }

    // Y(1) = CSPN*C1 + CSGN*C2 (Fortran lines 4771-4772)
    y[0] = cspn * c1 + csgn * c2;

    Ok(nz)
}
