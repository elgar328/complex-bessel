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
use crate::types::{BesselError, Scaling};
use crate::utils::zabs;

/// Analytic continuation for Airy functions.
///
/// Applies: K(fnu, zn*exp(mp)) = K(fnu, zn)*exp(-mp*fnu) - mp*I(fnu, zn)
/// where mp = pi*mr*i.
///
/// Same formula as ZACON but with higher-order recurrence removed.
/// FNU is 1/3 or 2/3, N is always 1.
///
/// Returns (y, nz) where nz < 0 indicates error.
pub(crate) fn zacai<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    mr: i32,
    n: usize,
    rl: T,
    tol: T,
    elim: T,
    alim: T,
) -> Result<(Vec<Complex<T>>, i32), BesselError> {
    let zero = T::zero();
    let one = T::one();
    let two = T::from(2.0).unwrap();
    let pi_t = T::from(PI).unwrap();

    let mut nz: i32 = 0;

    // ZN = -Z (Fortran line 4703-4704)
    let zn = Complex::new(-z.re, -z.im);
    let az = zabs(z);
    let nn = n;
    let dfnu = fnu + T::from((n - 1) as f64).unwrap();

    // I function dispatch (Fortran lines 4706-4731)
    // Direct calls to avoid ZBINU → ZBUNI → ZUNI2 → ZAIRY recursion
    let (mut y, _nw): (Vec<Complex<T>>, i32) =
        if az <= two || az * az * T::from(0.25).unwrap() <= dfnu + one {
            // Label 10: power series (Fortran line 4710)
            zseri(zn, fnu, kode, nn, tol, elim, alim)
        } else if az >= rl {
            // Label 20: asymptotic expansion for large z (Fortran line 4718)
            let result = zasyi(zn, fnu, kode, nn, rl, tol, elim, alim);
            if result.1 < 0 {
                // NW < 0 → error (Fortran line 4720 GO TO 80)
                let nz_err = if result.1 == -2 { -2 } else { -1 };
                return Err(if nz_err == -2 {
                    BesselError::ConvergenceFailure
                } else {
                    BesselError::Overflow
                });
            }
            result
        } else {
            // Label 30: Miller algorithm (Fortran line 4726)
            let result = zmlri(zn, fnu, kode, nn, tol);
            if result.1 < 0 {
                let nz_err = if result.1 == -2 { -2 } else { -1 };
                return Err(if nz_err == -2 {
                    BesselError::ConvergenceFailure
                } else {
                    BesselError::Overflow
                });
            }
            result
        };

    // Label 40: Analytic continuation (Fortran lines 4732-4777)
    // K function at -z (Fortran line 4736)
    let (cy, nw_k) = zbknu(zn, fnu, kode, 1, tol, elim, alim)?;
    if nw_k != 0 {
        // NW != 0 → error (Fortran line 4737 GO TO 80)
        return Err(BesselError::Overflow);
    }

    let fmr = T::from(mr as f64).unwrap();
    let sgn = -pi_t.copysign(fmr); // SGN = -DSIGN(PI, FMR) (Fortran line 4739)

    // CSGN = (0, sgn) (Fortran lines 4740-4741)
    let mut csgnr = zero;
    let mut csgni = sgn;

    if kode == Scaling::Exponential {
        // Fortran lines 4743-4745: KODE=2 adjustment
        let yy = -zn.im;
        let new_r = -csgni * yy.sin();
        let new_i = csgni * yy.cos();
        csgnr = new_r;
        csgni = new_i;
    }

    // CSPN = exp(FNU*PI*I) with precision preservation (Fortran lines 4751-4758)
    let inu = fnu.to_i32().unwrap();
    let arg = (fnu - T::from(inu as f64).unwrap()) * sgn;
    let mut cspnr = arg.cos();
    let mut cspni = arg.sin();
    if inu % 2 != 0 {
        cspnr = -cspnr;
        cspni = -cspni;
    }

    // C1 = K result, C2 = I result (Fortran lines 4760-4763)
    let mut c1r = cy[0].re;
    let mut c1i = cy[0].im;
    let mut c2r = y[0].re;
    let mut c2i = y[0].im;

    if kode == Scaling::Exponential {
        // Fortran lines 4765-4769: ZS1S2 call for KODE=2
        let iuf: i32 = 0;
        let ascle = T::from(1.0e3).unwrap() * T::MACH_TINY / tol;
        let s1s2_out = zs1s2(
            zn,
            Complex::new(c1r, c1i),
            Complex::new(c2r, c2i),
            ascle,
            alim,
            iuf,
        );
        nz += s1s2_out.nz;
        c1r = s1s2_out.s1.re;
        c1i = s1s2_out.s1.im;
        c2r = s1s2_out.s2.re;
        c2i = s1s2_out.s2.im;
    }

    // Y(1) = CSPN*C1 + CSGN*C2 (Fortran lines 4771-4772)
    y[0] = Complex::new(
        cspnr * c1r - cspni * c1i + csgnr * c2r - csgni * c2i,
        cspnr * c1i + cspni * c1r + csgnr * c2i + csgni * c2r,
    );

    Ok((y, nz))
}
