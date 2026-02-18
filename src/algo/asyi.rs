//! Asymptotic expansion for I Bessel function.
//!
//! Translation of Fortran ZASYI from TOMS 644 (zbsubs.f lines 3813-3978).
//! Computes I(fnu, z) for |z| > max(RL, fnu²/2).

#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]
#![allow(clippy::too_many_arguments)]

use num_complex::Complex;

use crate::algo::constants::PI;
use crate::machine::BesselFloat;
use crate::types::Scaling;
use crate::utils::{zabs, zdiv};

/// Asymptotic expansion of I Bessel function for large |z|.
///
/// Returns (y, nz) where:
/// - nz = 0: normal return
/// - nz = -1: overflow on KODE=1
/// - nz = -2: convergence failure
pub(crate) fn zasyi<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    y: &mut [Complex<T>],
    rl: T,
    tol: T,
    elim: T,
    alim: T,
) -> i32 {
    let zero = T::zero();
    let one = T::one();
    let eight = T::from_f64(8.0);
    let czero = Complex::new(zero, zero);

    // Fortran DATA constants (line 3834)
    let pi = T::from_f64(PI);
    let rtpi = T::from_f64(0.159154943091895336); // 1/(2*pi)

    let n = y.len();
    for v in y.iter_mut() {
        *v = czero;
    }
    let nz: i32 = 0;

    let az = zabs(z);
    let arm = T::from_f64(1.0e3) * T::MACH_TINY;
    let rtr1 = arm.sqrt();
    let il = if n < 2 { n } else { 2 };
    let dfnu0 = fnu + T::from_f64((n - il) as f64);

    // Overflow test (Fortran lines 3846-3858)
    let raz = one / az;
    let str = z.re * raz;
    let sti = -z.im * raz;
    let mut ak1r = rtpi * str * raz;
    let mut ak1i = rtpi * sti * raz;
    // zsqrt(ak1) - use standard formula
    let ak1_abs = zabs(Complex::new(ak1r, ak1i));
    if ak1_abs > zero {
        let ak1_arg = ak1i.atan2(ak1r);
        let sqrt_abs = ak1_abs.sqrt();
        let half_arg = ak1_arg * T::from_f64(0.5);
        ak1r = sqrt_abs * half_arg.cos();
        ak1i = sqrt_abs * half_arg.sin();
    }

    let mut czr = z.re;
    let czi = z.im;
    if kode == Scaling::Exponential {
        czr = zero;
    }

    if czr.abs() > elim {
        // Overflow (Fortran label 100, line 3972)
        return -1;
    }

    let dnu2 = dfnu0 + dfnu0;
    let mut koded: i32 = 1;
    if czr.abs() <= alim || n <= 2 {
        koded = 0;
        // zexp(cz) and multiply ak1
        let exp_r = czr.exp();
        let str_e = exp_r * czi.cos();
        let sti_e = exp_r * czi.sin();
        // zmlt(ak1, exp(cz))
        let new_r = ak1r * str_e - ak1i * sti_e;
        let new_i = ak1r * sti_e + ak1i * str_e;
        ak1r = new_r;
        ak1i = new_i;
    }

    let mut fdn = if dnu2 > rtr1 { dnu2 * dnu2 } else { zero };
    let ezr = z.re * eight;
    let ezi = z.im * eight;

    // Error test parameters (Fortran lines 3875-3877)
    let aez = eight * az;
    let s = tol / aez;
    let jl = (rl + rl).to_i32().unwrap() + 2;

    // Phase factor P1 (Fortran lines 3878-3895)
    let mut p1r = zero;
    let mut p1i = zero;
    if z.im != zero {
        // Compute exp(pi*(0.5+fnu+n-il)*i) minimizing precision loss
        let inu = fnu.to_i32().unwrap();
        let arg = (fnu - T::from_f64(inu as f64)) * pi;
        let inu_adj = inu + n as i32 - il as i32;
        let ak_sign = -arg.sin();
        let mut bk_sign = arg.cos();
        if z.im < zero {
            bk_sign = -bk_sign;
        }
        p1r = ak_sign;
        p1i = bk_sign;
        if inu_adj % 2 != 0 {
            p1r = -p1r;
            p1i = -p1i;
        }
    }

    // Main asymptotic expansion loop (Fortran DO 70, lines 3897-3948)
    let mut dfnu_val = dfnu0;
    for k in 0..il {
        let sqk0 = fdn - one;
        let atol = s * sqk0.abs();
        let mut sgn = one;
        let mut cs1r = one;
        let mut cs1i = zero;
        let mut cs2r = one;
        let mut cs2i = zero;
        let mut ckr = one;
        let mut cki = zero;
        let mut ak_val = zero;
        let mut aa = one;
        let mut bb = aez;
        let mut dkr = ezr;
        let mut dki = ezi;
        let mut sqk = sqk0;

        let mut converged = false;
        for _j in 1..=jl {
            // zdiv(ck, dk)
            let st = zdiv(Complex::new(ckr, cki), Complex::new(dkr, dki));
            ckr = st.re * sqk;
            cki = st.im * sqk;
            cs2r = cs2r + ckr;
            cs2i = cs2i + cki;
            sgn = -sgn;
            cs1r = cs1r + ckr * sgn;
            cs1i = cs1i + cki * sgn;
            dkr = dkr + ezr;
            dki = dki + ezi;
            aa = aa * sqk.abs() / bb;
            bb = bb + aez;
            ak_val = ak_val + eight;
            sqk = sqk - ak_val;
            if aa <= atol {
                converged = true;
                break;
            }
        }

        if !converged {
            // Convergence failure (Fortran label 110, line 3975)
            return -2;
        }

        // Label 50: combine CS1 and CS2 (Fortran lines 3931-3947)
        let mut s2r = cs1r;
        let mut s2i = cs1i;
        if z.re + z.re < elim {
            let tzr = z.re + z.re;
            let tzi = z.im + z.im;
            // zexp(-2z)
            let exp_r = (-tzr).exp();
            let mut str_e = exp_r * (-tzi).cos();
            let mut sti_e = exp_r * (-tzi).sin();
            // zmlt(exp(-2z), p1)
            let new_r = str_e * p1r - sti_e * p1i;
            let new_i = str_e * p1i + sti_e * p1r;
            str_e = new_r;
            sti_e = new_i;
            // zmlt(str, cs2)
            let new_r2 = str_e * cs2r - sti_e * cs2i;
            let new_i2 = str_e * cs2i + sti_e * cs2r;
            s2r = s2r + new_r2;
            s2i = s2i + new_i2;
        }

        // Update fdn, negate p1, store result (Fortran lines 3942-3947)
        fdn = fdn + eight * dfnu_val + T::from_f64(4.0);
        p1r = -p1r;
        p1i = -p1i;
        // M = N - IL + K + 1 (Fortran 1-based) → 0-based: n - il + k
        let m = n - il + k;
        y[m] = Complex::new(s2r * ak1r - s2i * ak1i, s2r * ak1i + s2i * ak1r);

        if k == 0 {
            dfnu_val = dfnu_val + one; // For next iteration (il=2)
        }
    }

    // Forward recurrence for remaining terms (Fortran lines 3949-3970)
    if n <= 2 {
        return nz;
    }

    let nn = n;
    let mut k_idx: isize = nn as isize - 3;
    let mut ak_rec = T::from_f64((nn - 2) as f64);
    let str2 = z.re * raz;
    let sti2 = -z.im * raz;
    let rzr = (str2 + str2) * raz;
    let rzi = (sti2 + sti2) * raz;

    for _ in 2..nn {
        let ki = k_idx as usize;
        y[ki] = Complex::new(
            (ak_rec + fnu) * (rzr * y[ki + 1].re - rzi * y[ki + 1].im) + y[ki + 2].re,
            (ak_rec + fnu) * (rzr * y[ki + 1].im + rzi * y[ki + 1].re) + y[ki + 2].im,
        );
        ak_rec = ak_rec - one;
        k_idx -= 1;
    }

    if koded != 0 {
        // Multiply by exp(cz) (Fortran lines 3964-3970)
        let exp_r2 = czr.exp();
        let ckr_e = exp_r2 * czi.cos();
        let cki_e = exp_r2 * czi.sin();
        for y_item in y.iter_mut().take(nn) {
            let str_val = y_item.re * ckr_e - y_item.im * cki_e;
            let sti_val = y_item.re * cki_e + y_item.im * ckr_e;
            *y_item = Complex::new(str_val, sti_val);
        }
    }

    nz
}
