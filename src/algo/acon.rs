//! Analytic continuation of K function to the left half z-plane.
//!
//! Translation of Fortran ZACON from TOMS 644 (zbsubs.f lines 4174-4377).
//! Applies: K(fnu, zn*exp(mp)) = K(fnu, zn)*exp(-mp*fnu) - mp*I(fnu, zn)
//! where mp = pi*mr*i.

#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]
#![allow(clippy::too_many_arguments)]
#![allow(unused_assignments)]

use num_complex::Complex;

use crate::algo::binu::zbinu;
use crate::algo::bknu::zbknu;
use crate::algo::constants::PI;
use crate::algo::s1s2::zs1s2;
use crate::machine::BesselFloat;
use crate::types::{BesselError, Scaling};
use crate::utils::zabs;

/// Analytic continuation of K function from right to left half-plane.
///
/// # Parameters
/// - `mr`: +1 or -1 (sign of Im(z) determines continuation direction)
pub(crate) fn zacon<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    mr: i32,
    n: usize,
    rl: T,
    fnul: T,
    tol: T,
    elim: T,
    alim: T,
) -> Result<(Vec<Complex<T>>, usize), BesselError> {
    let zero = T::zero();
    let one = T::one();
    let pi_t = T::from(PI).unwrap();

    let mut nz: usize = 0;

    // ZN = -Z (Fortran lines 4203-4204)
    let zn = Complex::new(-z.re, -z.im);
    let nn = n;

    // Compute I function at -z via ZBINU (Fortran lines 4206-4208)
    let (mut y, _nw) = zbinu(zn, fnu, kode, nn, rl, fnul, tol, elim, alim)?;
    // nw is returned as usize from zbinu; if it failed, we already got an Err

    // Compute K function at -z via ZBKNU (Fortran lines 4212-4213)
    let nn_k = if n < 2 { n } else { 2 };
    let (cy, nw_k) = zbknu(zn, fnu, kode, nn_k, tol, elim, alim)?;
    if nw_k != 0 {
        return Err(BesselError::Overflow);
    }

    // Analytic continuation formula (Fortran lines 4214-4276)
    let s1 = cy[0];
    let fmr = T::from(mr as f64).unwrap();
    let sgn = -pi_t.copysign(fmr); // -sign(pi, fmr)

    // CSGN = (0, sgn) (Fortran lines 4219-4220)
    let mut csgnr = zero;
    let mut csgni = sgn;

    if kode == Scaling::Exponential {
        // Multiply CSGN by exp(-i*zn.im) (Fortran lines 4222-4225)
        let yy = -zn.im;
        let cpn = yy.cos();
        let spn = yy.sin();
        let new_r = csgnr * cpn - csgni * spn;
        let new_i = csgnr * spn + csgni * cpn;
        csgnr = new_r;
        csgni = new_i;
    }

    // CSPN = exp(fnu*pi*i) with precision preservation (Fortran lines 4231-4239)
    let inu = fnu.to_i32().unwrap();
    let arg = (fnu - T::from(inu as f64).unwrap()) * sgn;
    let cpn = arg.cos();
    let spn = arg.sin();
    let mut cspnr = cpn;
    let mut cspni = spn;
    if inu % 2 != 0 {
        cspnr = -cspnr;
        cspni = -cspni;
    }

    // First two terms (Fortran lines 4241-4276)
    let mut iuf: i32 = 0;
    let ascle = T::from(1.0e3).unwrap() * T::MACH_TINY / tol;

    let mut c1r = s1.re;
    let mut c1i = s1.im;
    let mut c2r = y[0].re;
    let mut c2i = y[0].im;

    let mut sc1r = zero;
    let mut sc1i = zero;
    let mut sc2r = zero;
    let mut sc2i = zero;

    if kode != Scaling::Unscaled {
        let s1s2_out = zs1s2(
            zn,
            Complex::new(c1r, c1i),
            Complex::new(c2r, c2i),
            ascle,
            alim,
            iuf,
        );
        c1r = s1s2_out.s1.re;
        c1i = s1s2_out.s1.im;
        c2r = s1s2_out.s2.re;
        c2i = s1s2_out.s2.im;
        nz += s1s2_out.nz as usize;
        sc1r = c1r;
        sc1i = c1i;
        iuf = s1s2_out.iuf;
    }

    // Y(1) = CSPN*C1 + CSGN*C2 (Fortran lines 4253-4256)
    let str = cspnr * c1r - cspni * c1i;
    let sti = cspnr * c1i + cspni * c1r;
    let ptr = csgnr * c2r - csgni * c2i;
    let pti = csgnr * c2i + csgni * c2r;
    y[0] = Complex::new(str + ptr, sti + pti);

    if n == 1 {
        return Ok((y, nz));
    }

    // Second term (Fortran lines 4258-4275)
    cspnr = -cspnr;
    cspni = -cspni;
    let s2 = cy[1];
    c1r = s2.re;
    c1i = s2.im;
    c2r = y[1].re;
    c2i = y[1].im;

    if kode != Scaling::Unscaled {
        let s1s2_out = zs1s2(
            zn,
            Complex::new(c1r, c1i),
            Complex::new(c2r, c2i),
            ascle,
            alim,
            iuf,
        );
        c1r = s1s2_out.s1.re;
        c1i = s1s2_out.s1.im;
        c2r = s1s2_out.s2.re;
        c2i = s1s2_out.s2.im;
        nz += s1s2_out.nz as usize;
        sc2r = c1r;
        sc2i = c1i;
        iuf = s1s2_out.iuf;
    }

    let str2 = cspnr * c1r - cspni * c1i;
    let sti2 = cspnr * c1i + cspni * c1r;
    let ptr2 = csgnr * c2r - csgni * c2i;
    let pti2 = csgnr * c2i + csgni * c2r;
    y[1] = Complex::new(str2 + ptr2, sti2 + pti2);

    if n == 2 {
        return Ok((y, nz));
    }

    // Forward recurrence on K function for n > 2 (Fortran lines 4277-4371)
    cspnr = -cspnr;
    cspni = -cspni;

    let azn = zabs(zn);
    let razn = one / azn;
    let str_rz = zn.re * razn;
    let sti_rz = -zn.im * razn;
    let rzr = (str_rz + str_rz) * razn;
    let rzi = (sti_rz + sti_rz) * razn;

    let fn_val = fnu + one;
    let mut ckr = fn_val * rzr;
    let mut cki = fn_val * rzi;

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
    let mut s1r = s1.re * cssr[kflag];
    let mut s1i = s1.im * cssr[kflag];
    let mut s2r_k = s2.re * cssr[kflag];
    let mut s2i_k = s2.im * cssr[kflag];
    let mut csr = csrr[kflag];

    for y_item in y[2..n].iter_mut() {
        let str_k = s2r_k;
        let sti_k = s2i_k;
        s2r_k = ckr * str_k - cki * sti_k + s1r;
        s2i_k = ckr * sti_k + cki * str_k + s1i;
        s1r = str_k;
        s1i = sti_k;

        c1r = s2r_k * csr;
        c1i = s2i_k * csr;
        let mut str_c1 = c1r;
        let mut sti_c1 = c1i;
        c2r = y_item.re;
        c2i = y_item.im;

        if kode != Scaling::Unscaled && iuf >= 0 {
            let s1s2_out = zs1s2(
                zn,
                Complex::new(c1r, c1i),
                Complex::new(c2r, c2i),
                ascle,
                alim,
                iuf,
            );
            c1r = s1s2_out.s1.re;
            c1i = s1s2_out.s1.im;
            c2r = s1s2_out.s2.re;
            c2i = s1s2_out.s2.im;
            nz += s1s2_out.nz as usize;
            sc1r = sc2r;
            sc1i = sc2i;
            sc2r = c1r;
            sc2i = c1i;
            iuf = s1s2_out.iuf;

            if iuf == 3 {
                // IUF=3 special handling (Fortran lines 4339-4345)
                iuf = -4;
                s1r = sc1r * cssr[kflag];
                s1i = sc1i * cssr[kflag];
                s2r_k = sc2r * cssr[kflag];
                s2i_k = sc2i * cssr[kflag];
                str_c1 = sc2r;
                sti_c1 = sc2i;
            }
        }

        // Y(I) = CSPN*C1 + CSGN*C2 (Fortran lines 4347-4350)
        let ptr_k = cspnr * c1r - cspni * c1i;
        let pti_k = cspnr * c1i + cspni * c1r;
        *y_item = Complex::new(
            ptr_k + csgnr * c2r - csgni * c2i,
            pti_k + csgnr * c2i + csgni * c2r,
        );

        ckr = ckr + rzr;
        cki = cki + rzi;
        cspnr = -cspnr;
        cspni = -cspni;

        // KFLAG scaling check (Fortran lines 4355-4370)
        if kflag < 2 {
            let c1m = c1r.abs().max(c1i.abs());
            if c1m > bscle {
                kflag += 1;
                bscle = bry[kflag];
                s1r = s1r * csr;
                s1i = s1i * csr;
                s2r_k = str_c1;
                s2i_k = sti_c1;
                s1r = s1r * cssr[kflag];
                s1i = s1i * cssr[kflag];
                s2r_k = s2r_k * cssr[kflag];
                s2i_k = s2i_k * cssr[kflag];
                csr = csrr[kflag];
            }
        }
    }

    Ok((y, nz))
}
