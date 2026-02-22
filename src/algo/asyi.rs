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
use crate::utils::{reciprocal_z, zabs, zdiv};

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
    // Note: caller (zbinu) already zeroes the output buffer.
    let nz: i32 = 0;

    let az = zabs(z);
    let arm = T::from_f64(1.0e3) * T::MACH_TINY;
    let rtr1 = arm.sqrt();
    let il = if n < 2 { n } else { 2 };
    let dfnu0 = fnu + T::from_f64((n - il) as f64);

    // Overflow test (Fortran lines 3846-3858)
    let raz = one / az;
    let mut ak1 = (z.conj() * (rtpi * raz * raz)).sqrt();

    let cz = if kode == Scaling::Exponential {
        Complex::new(zero, z.im)
    } else {
        z
    };

    if cz.re.abs() > elim {
        // Overflow (Fortran label 100, line 3972)
        return -1;
    }

    let dnu2 = dfnu0 + dfnu0;
    let deferred_exp = cz.re.abs() > alim && n > 2;
    if !deferred_exp {
        // zexp(cz) and multiply ak1
        ak1 = ak1 * cz.exp();
    }

    let mut fdn = if dnu2 > rtr1 { dnu2 * dnu2 } else { zero };
    let ezr = z.re * eight;
    let ezi = z.im * eight;

    // Error test parameters (Fortran lines 3875-3877)
    let aez = eight * az;
    let s = tol / aez;
    // Safety: rl (limit) is finite and bounded
    let jl = (rl + rl).to_i32().unwrap() + 2;

    // Phase factor P1 (Fortran lines 3878-3895)
    let mut p1 = czero;
    if z.im != zero {
        // Compute exp(pi*(0.5+fnu+n-il)*i) minimizing precision loss
        // Safety: fnu is finite and < ~1e15 per upper-interface checks
        let inu = fnu.to_i32().unwrap();
        let arg = (fnu - T::from_f64(inu as f64)) * pi;
        let inu_adj = inu + n as i32 - il as i32;
        let ak_sign = -arg.sin();
        let mut bk_sign = arg.cos();
        if z.im < zero {
            bk_sign = -bk_sign;
        }
        p1 = Complex::new(ak_sign, bk_sign);
        if inu_adj % 2 != 0 {
            p1 = -p1;
        }
    }

    // Main asymptotic expansion loop (Fortran DO 70, lines 3897-3948)
    let mut dfnu_val = dfnu0;
    for k in 0..il {
        let sqk0 = fdn - one;
        let atol = s * sqk0.abs();
        let mut sgn = one;
        let mut cs1 = Complex::from(one);
        let mut cs2 = Complex::from(one);
        let mut ck = Complex::from(one);
        let mut ak_val = zero;
        let mut aa = one;
        let mut bb = aez;
        let mut dk = Complex::new(ezr, ezi);
        let mut sqk = sqk0;

        let mut converged = false;
        for _j in 1..=jl {
            ck = zdiv(ck, dk) * sqk;
            cs2 = cs2 + ck;
            sgn = -sgn;
            cs1 = cs1 + ck * sgn;
            dk = dk + Complex::new(ezr, ezi);
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
        let mut s2 = cs1;
        if z.re + z.re < elim {
            let exp_2z = (-(z + z)).exp();
            s2 = s2 + exp_2z * p1 * cs2;
        }

        // Update fdn, negate p1, store result (Fortran lines 3942-3947)
        fdn = fdn + eight * dfnu_val + T::from_f64(4.0);
        p1 = -p1;
        // M = N - IL + K + 1 (Fortran 1-based) → 0-based: n - il + k
        let m = n - il + k;
        y[m] = s2 * ak1;

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
    let rz = reciprocal_z(z);

    for _ in 2..nn {
        let ki = k_idx as usize;
        y[ki] = rz * y[ki + 1] * (ak_rec + fnu) + y[ki + 2];
        ak_rec = ak_rec - one;
        k_idx -= 1;
    }

    if deferred_exp {
        // Multiply by exp(cz) (Fortran lines 3964-3970)
        let cz_exp = cz.exp();
        for y_item in y.iter_mut().take(nn) {
            *y_item = *y_item * cz_exp;
        }
    }

    nz
}
