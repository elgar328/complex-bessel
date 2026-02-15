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
use crate::utils::zabs;

/// Output of ZBUNI.
pub(crate) struct BuniOutput<T> {
    /// Computed I function values.
    pub y: Vec<Complex<T>>,
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
/// - `n`: number of sequence members
/// - `nui`: number of extra orders to boost
/// - `fnul`: large-order threshold
/// - `tol`, `elim`, `alim`: machine-derived thresholds
pub(crate) fn zbuni<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    n: usize,
    nui: usize,
    fnul: T,
    tol: T,
    elim: T,
    alim: T,
) -> BuniOutput<T> {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);

    let nz: i32 = 0;

    // Region select (Fortran lines 6687-6690)
    let ax = z.re.abs() * T::from(1.7321).unwrap();
    let ay = z.im.abs();
    let iform: i32 = if ay > ax { 2 } else { 1 };

    // ── NUI == 0: direct call (Fortran lines 6691, 6815-6835) ──
    if nui == 0 {
        let (y_out, nw, nlast) = if iform == 1 {
            let result = zuni1(z, fnu, kode, n, fnul, tol, elim, alim);
            (result.y, result.nz, result.nlast as usize)
        } else {
            let result = zuni2(z, fnu, kode, n, fnul, tol, elim, alim);
            (result.y, result.nz, result.nlast as usize)
        };
        if nw < 0 {
            let nz_out = if nw == -2 { -2 } else { -1 };
            return BuniOutput {
                y: y_out,
                nz: nz_out,
                nlast: 0,
            };
        }
        return BuniOutput {
            y: y_out,
            nz: nw,
            nlast,
        };
    }

    // ── NUI > 0: boost order and recur backward (Fortran lines 6692-6809) ──
    let fnui = T::from(nui as f64).unwrap();
    let dfnu = fnu + T::from((n - 1) as f64).unwrap();
    let gnu = dfnu + fnui;

    // Compute 2 values at boosted order (Fortran lines 6700-6711)
    let mut cy = [czero; 2];
    let (_nlast_val, nw) = if iform == 1 {
        let result = zuni1(z, gnu, kode, 2, fnul, tol, elim, alim);
        cy[0] = result.y[0];
        cy[1] = result.y[1];
        (result.nlast as usize, result.nz)
    } else {
        let result = zuni2(z, gnu, kode, 2, fnul, tol, elim, alim);
        cy[0] = result.y[0];
        cy[1] = result.y[1];
        (result.nlast as usize, result.nz)
    };

    if nw < 0 {
        let nz_out = if nw == -2 { -2 } else { -1 };
        return BuniOutput {
            y: vec![czero; n],
            nz: nz_out,
            nlast: 0,
        };
    }
    if nw != 0 {
        // NLAST = N (Fortran label 90)
        return BuniOutput {
            y: vec![czero; n],
            nz: 0,
            nlast: n,
        };
    }

    let mut y = vec![czero; n];

    // ── Scale backward recurrence (Fortran lines 6714-6772) ──
    let str_val = zabs(Complex::new(cy[0].re, cy[0].im));
    let bry0 = T::from(1.0e3).unwrap() * T::MACH_TINY / tol;
    let bry1 = one / bry0;

    let mut iflag: usize;
    let mut ascle: T;
    let mut csclr: T;

    if str_val > bry0 {
        if str_val < bry1 {
            iflag = 2;
            ascle = bry1;
            csclr = one;
        } else {
            iflag = 3;
            ascle = bry1; // BRY(3)=BRY(2) in Fortran
            csclr = tol;
        }
    } else {
        iflag = 1;
        ascle = bry0;
        csclr = one / tol;
    }

    let mut cscrr = one / csclr;
    let mut s1 = Complex::new(cy[1].re * csclr, cy[1].im * csclr);
    let mut s2 = Complex::new(cy[0].re * csclr, cy[0].im * csclr);

    let raz = one / zabs(z);
    let str2 = z.re * raz;
    let sti2 = -z.im * raz;
    let rzr = (str2 + str2) * raz;
    let rzi = (sti2 + sti2) * raz;

    // Backward recurrence for NUI steps (Fortran DO 30)
    let mut fnui_val = fnui;
    for _i in 0..nui {
        let str = s2.re;
        let sti = s2.im;
        let cfn = dfnu + fnui_val;
        s2 = Complex::new(
            cfn * (rzr * str - rzi * sti) + s1.re,
            cfn * (rzr * sti + rzi * str) + s1.im,
        );
        s1 = Complex::new(str, sti);
        fnui_val = fnui_val - one;

        if iflag >= 3 {
            continue;
        }
        let str_s = s2.re * cscrr;
        let sti_s = s2.im * cscrr;
        let c1m = str_s.abs().max(sti_s.abs());
        if c1m <= ascle {
            continue;
        }
        iflag += 1;
        ascle = bry1; // BRY(2)=BRY(3) in Fortran
        s1 = Complex::new(s1.re * cscrr, s1.im * cscrr);
        s2 = Complex::new(str_s, sti_s);
        csclr = csclr * tol;
        cscrr = one / csclr;
        s1 = Complex::new(s1.re * csclr, s1.im * csclr);
        s2 = Complex::new(s2.re * csclr, s2.im * csclr);
    }

    // Store Y(N) (Fortran line 6773)
    y[n - 1] = Complex::new(s2.re * cscrr, s2.im * cscrr);
    if n == 1 {
        return BuniOutput { y, nz, nlast: 0 };
    }

    // Backward recurrence for remaining N-1 orders (Fortran DO 40, lines 6779-6809)
    let nl = n - 1;
    let mut fnui_val = T::from(nl as f64).unwrap();
    let mut k = nl; // 1-based index counting down
    for _i in 0..nl {
        let str = s2.re;
        let sti = s2.im;
        let cfn = fnu + fnui_val;
        s2 = Complex::new(
            cfn * (rzr * str - rzi * sti) + s1.re,
            cfn * (rzr * sti + rzi * str) + s1.im,
        );
        s1 = Complex::new(str, sti);
        let str_s = s2.re * cscrr;
        let sti_s = s2.im * cscrr;
        y[k - 1] = Complex::new(str_s, sti_s);
        fnui_val = fnui_val - one;
        k -= 1;

        if iflag >= 3 {
            continue;
        }
        let c1m = str_s.abs().max(sti_s.abs());
        if c1m <= ascle {
            continue;
        }
        iflag += 1;
        ascle = bry1; // BRY(2)=BRY(3) in Fortran
        s1 = Complex::new(s1.re * cscrr, s1.im * cscrr);
        s2 = Complex::new(str_s * csclr, sti_s * csclr); // Note: s2 = scaled * csclr
        csclr = csclr * tol;
        cscrr = one / csclr;
        s1 = Complex::new(s1.re * csclr, s1.im * csclr);
        s2 = Complex::new(s2.re, s2.im); // already rescaled above
    }

    BuniOutput { y, nz, nlast: 0 }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    const TOL: f64 = 2.220446049250313e-16;
    const ELIM: f64 = 700.9217936944459;
    const ALIM: f64 = 664.8716455337102;
    const FNUL: f64 = 85.92135864716212;

    #[test]
    fn zbuni_nui_zero() {
        // nui=0: direct passthrough to ZUNI1/ZUNI2
        let z = Complex64::new(3.0, 1.0);
        let result = zbuni(z, 90.0, Scaling::Unscaled, 1, 0, FNUL, TOL, ELIM, ALIM);
        assert!(result.nz >= 0);
    }

    #[test]
    fn zbuni_nui_positive() {
        // nui > 0: order boosting + backward recurrence
        let z = Complex64::new(50.0, 30.0);
        let result = zbuni(z, 50.0, Scaling::Unscaled, 2, 40, FNUL, TOL, ELIM, ALIM);
        assert!(result.nz >= 0, "nz = {}", result.nz);
    }

    #[test]
    fn zbuni_region2() {
        // |Im(z)| > |Re(z)|*sqrt(3) → iform=2 → ZUNI2
        let z = Complex64::new(1.0, 10.0);
        let result = zbuni(z, 90.0, Scaling::Unscaled, 1, 0, FNUL, TOL, ELIM, ALIM);
        assert!(result.nz >= 0);
    }

    #[test]
    fn zbuni_region1() {
        // |Im(z)| <= |Re(z)|*sqrt(3) → iform=1 → ZUNI1
        let z = Complex64::new(5.0, 1.0);
        let result = zbuni(z, 90.0, Scaling::Unscaled, 1, 0, FNUL, TOL, ELIM, ALIM);
        assert!(result.nz >= 0);
    }
}
