//! Power series for I Bessel function.
//!
//! Translation of Fortran ZSERI from TOMS 644 (zbsubs.f lines 3622-3812).
//! Computes I(fnu, z) for |z| <= 2*sqrt(fnu+1) region.

#![allow(clippy::excessive_precision)]

use num_complex::Complex;

use crate::algo::gamln::gamln;
use crate::algo::uchk::zuchk;
use crate::machine::BesselFloat;
use crate::types::Scaling;
use crate::utils::{zabs, zdiv};

/// Power series computation of I Bessel function.
///
/// Returns (y, nz) where:
/// - nz > 0: last nz components set to zero due to underflow
/// - nz < 0: |nz| components underflowed but |z²/4| > fnu+n-nz-1,
///   must complete in zbinu with n = n - |nz|
pub(crate) fn zseri<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    n: usize,
    tol: T,
    elim: T,
    alim: T,
) -> (Vec<Complex<T>>, i32) {
    let zero = T::zero();
    let one = T::one();
    let two = T::from(2.0).unwrap();
    let half = T::from(0.5).unwrap();
    let czero = Complex::new(zero, zero);
    let cone = Complex::new(one, zero);

    let mut y = vec![czero; n];
    let mut nz: i32 = 0;

    let az = zabs(z);

    // z = 0 special case (Fortran label 160, line 3792)
    if az == zero {
        if fnu == zero {
            y[0] = cone;
        }
        nz = if fnu == zero { n as i32 - 1 } else { n as i32 };
        return (y, nz);
    }

    // Fortran lines 3651-3654
    let arm = T::from(1.0e3).unwrap() * T::MACH_TINY;
    let rtr1 = arm.sqrt();
    let mut crscr = one;
    let mut iflag: i32 = 0;

    // |z| very tiny (Fortran label 150, line 3789)
    if az < arm {
        nz = if fnu == zero { n as i32 - 1 } else { n as i32 };
        if fnu == zero {
            y[0] = cone;
        }
        return (y, nz);
    }

    // hz = z/2, cz = hz² (Fortran lines 3656-3663)
    let hz = Complex::new(z.re * half, z.im * half);
    let cz = if az > rtr1 {
        // zmlt(hz, hz)
        Complex::new(hz.re * hz.re - hz.im * hz.im, two * hz.re * hz.im)
    } else {
        czero
    };
    let acz = zabs(cz);
    let mut nn = n;

    // ck = log(hz) (Fortran line 3665)
    let hz_r = zabs(hz);
    let ck = Complex::new(hz_r.ln(), hz.im.atan2(hz.re));

    // Outer loop: underflow scan (Fortran labels 20-30)
    'outer: loop {
        let mut dfnu = fnu + T::from((nn - 1) as f64).unwrap();
        let fnup = dfnu + one;

        // Underflow test (Fortran lines 3672-3677)
        let mut ak1r = ck.re * dfnu;
        let ak1i = ck.im * dfnu;
        let ak = gamln(fnup).unwrap();
        ak1r = ak1r - ak;
        if kode == Scaling::Exponential {
            ak1r = ak1r - z.re;
        }

        if ak1r <= -elim {
            // Label 30: underflow (Fortran lines 3678-3685)
            nz += 1;
            y[nn - 1] = czero;
            if acz > dfnu {
                // Label 190: nz = -nz
                nz = -nz;
                return (y, nz);
            }
            nn -= 1;
            if nn == 0 {
                return (y, nz);
            }
            continue 'outer;
        }

        // Label 40: check near-underflow scaling (Fortran lines 3687-3691)
        let mut ss = one;
        let mut ascle = zero;
        if ak1r <= -alim {
            iflag = 1;
            ss = one / tol;
            crscr = tol;
            ascle = arm * ss;
        }

        // Label 50: compute coefficient (Fortran lines 3693-3698)
        let mut aa = ak1r.exp();
        if iflag == 1 {
            aa = aa * ss;
        }
        let mut coefr = aa * ak1i.cos();
        let mut coefi = aa * ak1i.sin();
        let atol = tol * acz / fnup;

        let il = if nn < 2 { nn } else { 2 };

        // Working array for up to 2 scaled values (Fortran W array)
        let mut w = [czero; 2];

        // Series computation loop (Fortran DO 90, lines 3699-3738)
        let mut went_to_30 = false;
        for (i, w_item) in w.iter_mut().enumerate().take(il) {
            dfnu = fnu + T::from((nn - 1 - i) as f64).unwrap();
            let fnup_i = dfnu + one;
            let mut s1 = cone;

            if acz >= tol * fnup_i {
                // Power series sum (Fortran lines 3705-3721)
                let mut ak1 = cone;
                let mut ak_val = fnup_i + two;
                let mut s = fnup_i;
                let mut aa_val = two;

                loop {
                    let rs = one / s;
                    let str = ak1.re * cz.re - ak1.im * cz.im;
                    let sti = ak1.re * cz.im + ak1.im * cz.re;
                    ak1 = Complex::new(str * rs, sti * rs);
                    s1 = Complex::new(s1.re + ak1.re, s1.im + ak1.im);
                    s = s + ak_val;
                    ak_val = ak_val + two;
                    aa_val = aa_val * acz * rs;
                    if aa_val <= atol {
                        break;
                    }
                }
            }

            // s2 = s1 * coef (Fortran lines 3723-3724)
            let s2 = Complex::new(s1.re * coefr - s1.im * coefi, s1.re * coefi + s1.im * coefr);
            *w_item = s2;

            if iflag != 0 {
                // Underflow check (Fortran lines 3728-3729)
                let underflowed = zuchk(s2, ascle, tol);
                if underflowed {
                    // Goto label 30
                    nz += 1;
                    y[nn - 1] = czero;
                    if acz > dfnu {
                        nz = -nz;
                        return (y, nz);
                    }
                    nn -= 1;
                    if nn == 0 {
                        return (y, nz);
                    }
                    went_to_30 = true;
                    break;
                }
            }

            // M = NN - I + 1 (Fortran 1-based) → 0-based: nn - 1 - i
            let m = nn - 1 - i;
            y[m] = Complex::new(s2.re * crscr, s2.im * crscr);

            if i < il - 1 {
                // Update coefficient (Fortran lines 3735-3737)
                let st = zdiv(Complex::new(coefr, coefi), hz);
                coefr = st.re * dfnu;
                coefi = st.im * dfnu;
            }
        }

        if went_to_30 {
            continue 'outer;
        }

        // Forward recurrence for remaining terms (Fortran lines 3739-3788)
        if nn <= 2 {
            return (y, nz);
        }

        let mut k: isize = nn as isize - 3; // 0-based (Fortran K=NN-2 → k=NN-3)
        let mut ak_val = T::from((nn - 2) as f64).unwrap();
        let raz = one / az;
        let str = z.re * raz;
        let sti = -z.im * raz;
        let rzr = (str + str) * raz;
        let rzi = (sti + sti) * raz;

        if iflag == 1 {
            // Scaled recurrence (Fortran label 120, lines 3760-3788)
            let mut s1 = w[0];
            let mut s2 = w[1];

            let mut l = 2usize; // iteration index (Fortran L=3..NN, we use 2..nn-1)
            while l < nn {
                let ck_val = s2;
                s2 = Complex::new(
                    s1.re + (ak_val + fnu) * (rzr * ck_val.re - rzi * ck_val.im),
                    s1.im + (ak_val + fnu) * (rzr * ck_val.im + rzi * ck_val.re),
                );
                s1 = ck_val;
                let ck_scaled = Complex::new(s2.re * crscr, s2.im * crscr);
                y[k as usize] = ck_scaled;
                ak_val = ak_val - one;
                k -= 1;
                l += 1;

                if zabs(ck_scaled) > ascle {
                    // Label 140: switch to unscaled recurrence (Fortran lines 3785-3788)
                    // IB = L + 1 in Fortran; remaining iterations continue unscaled
                    while l < nn {
                        let ki = k as usize;
                        y[ki] = Complex::new(
                            (ak_val + fnu) * (rzr * y[ki + 1].re - rzi * y[ki + 1].im)
                                + y[ki + 2].re,
                            (ak_val + fnu) * (rzr * y[ki + 1].im + rzi * y[ki + 1].re)
                                + y[ki + 2].im,
                        );
                        ak_val = ak_val - one;
                        k -= 1;
                        l += 1;
                    }
                    return (y, nz);
                }
            }
            return (y, nz);
        }

        // Unscaled recurrence (Fortran label 100, lines 3749-3756)
        for _ in 2..nn {
            let ki = k as usize;
            y[ki] = Complex::new(
                (ak_val + fnu) * (rzr * y[ki + 1].re - rzi * y[ki + 1].im) + y[ki + 2].re,
                (ak_val + fnu) * (rzr * y[ki + 1].im + rzi * y[ki + 1].re) + y[ki + 2].im,
            );
            ak_val = ak_val - one;
            k -= 1;
        }
        return (y, nz);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn seri_i0_at_origin() {
        // I_0(0) = 1
        let z = Complex64::new(0.0, 0.0);
        let (y, nz) = zseri(
            z,
            0.0,
            Scaling::Unscaled,
            1,
            f64::tol(),
            f64::elim(),
            f64::alim(),
        );
        assert_eq!(y[0].re, 1.0);
        assert_eq!(y[0].im, 0.0);
        assert_eq!(nz, 0);
    }
}
