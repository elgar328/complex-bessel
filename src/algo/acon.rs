//! Analytic continuation of K function to the left half z-plane.
//!
//! Translation of Fortran ZACON from TOMS 644 (zbsubs.f lines 4174-4377).
//! Applies: K(fnu, zn*exp(mp)) = K(fnu, zn)*exp(-mp*fnu) - mp*I(fnu, zn)
//! where mp = pi*mr*i.

#![allow(clippy::too_many_arguments)]

use num_complex::Complex;

use crate::algo::binu::zbinu;
use crate::algo::bknu::zbknu;
use crate::algo::constants::PI;
use crate::algo::s1s2::zs1s2;
use crate::machine::BesselFloat;
use crate::types::{Error, Scaling};
use crate::utils::{reciprocal_z, zabs};

/// Analytic continuation of K function from right to left half-plane.
///
/// Writes results into `y` and returns nz (underflow count).
///
/// # Parameters
/// - `mr`: +1 or -1 (sign of Im(z) determines continuation direction)
pub(crate) fn zacon<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    mr: i32,
    y: &mut [Complex<T>],
    rl: T,
    fnul: T,
    tol: T,
    elim: T,
    alim: T,
) -> Result<usize, Error> {
    let zero = T::zero();
    let one = T::one();
    let pi_t = T::from_f64(PI);
    let czero = Complex::new(zero, zero);

    let n = y.len();
    let mut nz: usize = 0;

    // ZN = -Z (Fortran lines 4203-4204)
    let zn = -z;

    // Compute I function at -z via ZBINU, written directly into y (Fortran lines 4206-4208)
    zbinu(zn, fnu, kode, y, rl, fnul, tol, elim, alim)?;

    // Compute K function at -z via ZBKNU (Fortran lines 4212-4213)
    let nn_k = if n < 2 { n } else { 2 };
    let mut k_buf = [czero; 2];
    let nw_k = zbknu(zn, fnu, kode, &mut k_buf[..nn_k], tol, elim, alim)?;
    if nw_k != 0 {
        return Err(Error::Overflow);
    }

    // Analytic continuation formula (Fortran lines 4214-4276)
    let s1 = k_buf[0];
    let fmr = T::from_f64(mr as f64);
    let sgn = -pi_t.copysign(fmr); // -sign(pi, fmr)

    // CSGN = (0, sgn) (Fortran lines 4219-4220)
    let mut csgn = Complex::new(zero, sgn);

    if kode == Scaling::Exponential {
        // Multiply CSGN by exp(i*yy) (Fortran lines 4222-4225)
        let yy = -zn.im;
        csgn = csgn * Complex::new(yy.cos(), yy.sin());
    }

    // CSPN = exp(fnu*pi*i) with precision preservation (Fortran lines 4231-4239)
    // Safety: fnu is finite and < ~1e15 per upper-interface checks
    let inu = fnu.to_i32().unwrap();
    let arg = (fnu - T::from_f64(inu as f64)) * sgn;
    let mut cspn = Complex::new(arg.cos(), arg.sin());
    if inu % 2 != 0 {
        cspn = -cspn;
    }

    // First two terms (Fortran lines 4241-4276)
    let mut iuf: i32 = 0;
    let ascle = T::from_f64(1.0e3) * T::MACH_TINY / tol;

    let mut c1 = s1;
    let mut c2 = y[0]; // I value from zbinu
    // sc1/sc2 track S1S2 outputs across iterations for IUF=3 recovery
    // (Fortran lines 4339-4345). Declared here to match Fortran variable scope;
    // the "unused" first assignment mirrors the Fortran initialization pattern.
    #[allow(unused_assignments)]
    let mut sc1 = czero;
    #[allow(unused_assignments)]
    let mut sc2 = czero;

    if kode != Scaling::Unscaled {
        let s1s2_out = zs1s2(zn, c1, c2, ascle, alim, iuf);
        c1 = s1s2_out.s1;
        c2 = s1s2_out.s2;
        nz += s1s2_out.nz as usize;
        iuf = s1s2_out.iuf;
    }

    // Y(1) = CSPN*C1 + CSGN*C2 (Fortran lines 4253-4256)
    y[0] = cspn * c1 + csgn * c2;

    if n == 1 {
        return Ok(nz);
    }

    // Second term (Fortran lines 4258-4275)
    cspn = -cspn;
    let s2 = k_buf[1];
    c1 = s2;
    c2 = y[1]; // I value from zbinu

    if kode != Scaling::Unscaled {
        let s1s2_out = zs1s2(zn, c1, c2, ascle, alim, iuf);
        c1 = s1s2_out.s1;
        c2 = s1s2_out.s2;
        nz += s1s2_out.nz as usize;
        sc2 = c1;
        iuf = s1s2_out.iuf;
    }

    y[1] = cspn * c1 + csgn * c2;

    if n == 2 {
        return Ok(nz);
    }

    // Forward recurrence on K function for n > 2 (Fortran lines 4277-4371)
    cspn = -cspn;

    let rz = reciprocal_z(zn);

    let fn_val = fnu + one;
    let mut ck = rz * fn_val;

    // Scale near exponent extremes (Fortran lines 4291-4316)
    let cscl = one / tol;
    let cscr = tol;
    let cssr = [cscl, one, cscr];
    let csrr = [cscr, one, cscl];
    let bry = [ascle, one / ascle, T::MACH_HUGE];

    let as2 = zabs(s2);
    let mut kflag: usize = if as2 > bry[0] {
        if as2 < bry[1] { 1 } else { 2 }
    } else {
        0
    };

    let mut bscle = bry[kflag];
    let mut s1_k = s1 * cssr[kflag];
    let mut s2_k = s2 * cssr[kflag];
    let mut csr = csrr[kflag];

    for y_item in y[2..n].iter_mut() {
        let prev = s2_k;
        s2_k = ck * prev + s1_k;
        s1_k = prev;

        c1 = s2_k * csr;
        let mut saved_c1 = c1;
        c2 = *y_item; // I value from zbinu

        if kode != Scaling::Unscaled && iuf >= 0 {
            let s1s2_out = zs1s2(zn, c1, c2, ascle, alim, iuf);
            c1 = s1s2_out.s1;
            c2 = s1s2_out.s2;
            nz += s1s2_out.nz as usize;
            sc1 = sc2;
            sc2 = c1;
            iuf = s1s2_out.iuf;

            if iuf == 3 {
                // IUF=3 special handling (Fortran lines 4339-4345)
                iuf = -4;
                s1_k = sc1 * cssr[kflag];
                s2_k = sc2 * cssr[kflag];
                saved_c1 = sc2;
            }
        }

        // Y(I) = CSPN*C1 + CSGN*C2 (Fortran lines 4347-4350)
        *y_item = cspn * c1 + csgn * c2;

        ck = ck + rz;
        cspn = -cspn;

        // KFLAG scaling check (Fortran lines 4355-4370)
        if kflag < 2 {
            let c1m = c1.re.abs().max(c1.im.abs());
            if c1m > bscle {
                kflag += 1;
                bscle = bry[kflag];
                s1_k = s1_k * csr;
                s2_k = saved_c1;
                s1_k = s1_k * cssr[kflag];
                s2_k = s2_k * cssr[kflag];
                csr = csrr[kflag];
            }
        }
    }

    Ok(nz)
}
