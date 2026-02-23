//! I function dispatcher for large orders with order boosting.
//!
//! Translation of Fortran ZBUNI from TOMS 644 / SLATEC (zbsubs.f lines 6665-6839).
//! Dispatches to ZUNI1/ZUNI2 and performs backward recurrence to reduce the order.

#![allow(clippy::too_many_arguments)]

use num_complex::Complex;

use crate::algo::uni1::zuni1;
use crate::algo::uni2::zuni2;
use crate::machine::BesselFloat;
use crate::types::Scaling;
use crate::utils::{mul_add_scalar, reciprocal_z, zabs};

/// Output metadata of ZBUNI (no Vec — results written into caller-provided slice).
#[derive(Debug, Clone, Copy)]
pub(crate) struct BuniOutput {
    /// Underflow count. -1 or -2 indicates overflow/convergence failure.
    pub nz: i32,
    /// If nonzero, remaining items need a different method.
    pub nlast: usize,
}

/// Compute I(fnu,z) for large orders via uniform asymptotic expansion + backward recurrence.
///
/// Equivalent to Fortran ZBUNI in TOMS 644 (zbsubs.f lines 6665-6839).
///
/// # Parameters
/// - `z`: complex argument (Re(z) >= 0)
/// - `fnu`: starting order
/// - `kode`: scaling mode
/// - `y`: output slice for sequence members (length determines n)
/// - `nui`: number of extra orders to boost
/// - `fnul`: large-order threshold
/// - `tol`, `elim`, `alim`: machine-derived thresholds
pub(crate) fn zbuni<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    y: &mut [Complex<T>],
    nui: usize,
    fnul: T,
    tol: T,
    elim: T,
    alim: T,
) -> BuniOutput {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);

    let n = y.len();
    y.fill(czero);
    let nz: i32 = 0;

    // Region select (Fortran lines 6687-6690)
    let ax = z.re.abs() * T::from_f64(1.7321);
    let ay = z.im.abs();
    let iform: i32 = if ay > ax { 2 } else { 1 };

    // ── NUI == 0: direct call (Fortran lines 6691, 6815-6835) ──
    if nui == 0 {
        let (nw, nlast) = if iform == 1 {
            let result = zuni1(z, fnu, kode, y, fnul, tol, elim, alim);
            (result.nz, result.nlast as usize)
        } else {
            let result = zuni2(z, fnu, kode, y, fnul, tol, elim, alim);
            (result.nz, result.nlast as usize)
        };
        if nw < 0 {
            let nz_out = if nw == -2 { -2 } else { -1 };
            return BuniOutput {
                nz: nz_out,
                nlast: 0,
            };
        }
        return BuniOutput { nz: nw, nlast };
    }

    // ── NUI > 0: boost order and recur backward (Fortran lines 6692-6809) ──
    let fnui = T::from_f64(nui as f64);
    let dfnu = fnu + T::from_f64((n - 1) as f64);
    let gnu = dfnu + fnui;

    // Compute 2 values at boosted order (Fortran lines 6700-6711)
    let mut cy = [czero; 2];
    let (_nlast_val, nw) = if iform == 1 {
        let result = zuni1(z, gnu, kode, &mut cy, fnul, tol, elim, alim);
        (result.nlast as usize, result.nz)
    } else {
        let result = zuni2(z, gnu, kode, &mut cy, fnul, tol, elim, alim);
        (result.nlast as usize, result.nz)
    };

    if nw < 0 {
        let nz_out = if nw == -2 { -2 } else { -1 };
        return BuniOutput {
            nz: nz_out,
            nlast: 0,
        };
    }
    if nw != 0 {
        // NLAST = N (Fortran label 90)
        return BuniOutput { nz: 0, nlast: n };
    }

    // ── Scale backward recurrence (Fortran lines 6714-6772) ──
    let str_val = zabs(cy[0]);
    let bry0 = T::from_f64(1.0e3) * T::MACH_TINY / tol;
    let bry1 = one / bry0;

    let (mut iflag, mut ascle, mut csclr) = if str_val > bry0 {
        if str_val < bry1 {
            (2, bry1, one)
        } else {
            (3, bry1, tol) // BRY(3)=BRY(2) in Fortran
        }
    } else {
        (1, bry0, one / tol)
    };

    let mut cscrr = one / csclr;
    let mut s1 = cy[1] * csclr;
    let mut s2 = cy[0] * csclr;

    let rz = reciprocal_z(z);

    // Backward recurrence for NUI steps (Fortran DO 30)
    let mut fnui_val = fnui;
    for _i in 0..nui {
        let prev = s2;
        let cfn = dfnu + fnui_val;
        s2 = mul_add_scalar(rz * prev, cfn, s1);
        s1 = prev;
        fnui_val = fnui_val - one;

        if iflag >= 3 {
            continue;
        }
        let c2_scaled = s2 * cscrr;
        let c1m = c2_scaled.re.abs().max(c2_scaled.im.abs());
        if c1m <= ascle {
            continue;
        }
        iflag += 1;
        ascle = bry1; // BRY(2)=BRY(3) in Fortran
        s1 = s1 * cscrr;
        s2 = c2_scaled;
        csclr = csclr * tol;
        cscrr = one / csclr;
        s1 = s1 * csclr;
        s2 = s2 * csclr;
    }

    // Store Y(N) (Fortran line 6773)
    y[n - 1] = s2 * cscrr;
    if n == 1 {
        return BuniOutput { nz, nlast: 0 };
    }

    // Backward recurrence for remaining N-1 orders (Fortran DO 40, lines 6779-6809)
    let nl = n - 1;
    let mut fnui_val = T::from_f64(nl as f64);
    let mut k = nl; // 1-based index counting down
    for _i in 0..nl {
        let prev = s2;
        let cfn = fnu + fnui_val;
        s2 = mul_add_scalar(rz * prev, cfn, s1);
        s1 = prev;
        let c2_scaled = s2 * cscrr;
        y[k - 1] = c2_scaled;
        fnui_val = fnui_val - one;
        k -= 1;

        if iflag >= 3 {
            continue;
        }
        let c1m = c2_scaled.re.abs().max(c2_scaled.im.abs());
        if c1m <= ascle {
            continue;
        }
        iflag += 1;
        ascle = bry1; // BRY(2)=BRY(3) in Fortran
        s1 = s1 * cscrr;
        s2 = c2_scaled * csclr;
        csclr = csclr * tol;
        cscrr = one / csclr;
        s1 = s1 * csclr;
    }

    BuniOutput { nz, nlast: 0 }
}
