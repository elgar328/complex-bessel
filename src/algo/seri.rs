//! Power series for I Bessel function.
//!
//! Translation of Fortran ZSERI from TOMS 644 (zbsubs.f lines 3622-3812).
//! Computes I(fnu, z) for |z| <= 2*sqrt(fnu+1) region.

use num_complex::Complex;

use crate::algo::gamln::gamln;
use crate::algo::uchk::zuchk;
use crate::machine::BesselFloat;
use crate::types::Scaling;
use crate::utils::{mul_add_scalar, reciprocal_z, zabs, zdiv};

/// Power series computation of I Bessel function.
///
/// Writes results into `y` and returns `nz` where:
/// - nz > 0: last nz components set to zero due to underflow
/// - nz < 0: |nz| components underflowed but |z²/4| > fnu+n-nz-1,
///   must complete in zbinu with n = n - |nz|
pub(crate) fn zseri<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    y: &mut [Complex<T>],
    tol: T,
    elim: T,
    alim: T,
) -> i32 {
    let zero = T::zero();
    let one = T::one();
    let two = T::from_f64(2.0);
    let half = T::from_f64(0.5);
    let czero = Complex::new(zero, zero);
    let cone = Complex::from(one);

    let n = y.len();
    // Note: caller (zbinu) already zeroes the output buffer.
    let mut nz: i32 = 0;

    let az = zabs(z);

    // z = 0 special case (Fortran label 160, line 3792)
    if az == zero {
        if fnu == zero {
            y[0] = cone;
            nz = n as i32 - 1;
        } else {
            nz = n as i32;
        }
        return nz;
    }

    // Fortran lines 3651-3654
    let arm = T::from_f64(1.0e3) * T::MACH_TINY;
    let rtr1 = arm.sqrt();
    let mut crscr = one;
    let mut iflag: i32 = 0;

    // |z| very tiny (Fortran label 150, line 3789)
    if az < arm {
        nz = if fnu == zero { n as i32 - 1 } else { n as i32 };
        if fnu == zero {
            y[0] = cone;
        }
        return nz;
    }

    // hz = z/2, cz = hz² (Fortran lines 3656-3663)
    let hz = z * half;
    let cz = if az > rtr1 { hz * hz } else { czero };
    let acz = zabs(cz);
    let mut nn = n;

    // ck = log(hz) (Fortran line 3665)
    let ck = hz.ln();

    // Outer loop: underflow scan (Fortran labels 20-30)
    'outer: loop {
        let mut dfnu = fnu + T::from_f64((nn - 1) as f64);
        let fnup = dfnu + one;

        // Underflow test (Fortran lines 3672-3677)
        let mut ak1r = ck.re * dfnu;
        let ak1i = ck.im * dfnu;
        // Safety: fnup > 0 guaranteed by algorithm invariant
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
                return nz;
            }
            nn -= 1;
            if nn == 0 {
                return nz;
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
        let mut coef = Complex::new(ak1r, ak1i).exp();
        if iflag == 1 {
            coef = coef * ss;
        }
        let atol = tol * acz / fnup;

        let il = if nn < 2 { nn } else { 2 };

        // Working array for up to 2 scaled values (Fortran W array)
        let mut w = [czero; 2];

        // Series computation loop (Fortran DO 90, lines 3699-3738)
        let mut went_to_30 = false;
        for (i, w_item) in w.iter_mut().enumerate().take(il) {
            dfnu = fnu + T::from_f64((nn - 1 - i) as f64);
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
                    ak1 = ak1 * cz * rs;
                    s1 = s1 + ak1;
                    s = s + ak_val;
                    ak_val = ak_val + two;
                    aa_val = aa_val * acz * rs;
                    if aa_val <= atol {
                        break;
                    }
                }
            }

            // s2 = s1 * coef (Fortran lines 3723-3724)
            let s2 = s1 * coef;
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
                        return nz;
                    }
                    nn -= 1;
                    if nn == 0 {
                        return nz;
                    }
                    went_to_30 = true;
                    break;
                }
            }

            // M = NN - I + 1 (Fortran 1-based) → 0-based: nn - 1 - i
            let m = nn - 1 - i;
            y[m] = s2 * crscr;

            if i < il - 1 {
                // Update coefficient (Fortran lines 3735-3737)
                coef = zdiv(coef, hz) * dfnu;
            }
        }

        if went_to_30 {
            continue 'outer;
        }

        // Forward recurrence for remaining terms (Fortran lines 3739-3788)
        if nn <= 2 {
            return nz;
        }

        let mut k: isize = nn as isize - 3; // 0-based (Fortran K=NN-2 → k=NN-3)
        let mut ak_val = T::from_f64((nn - 2) as f64);
        let rz = reciprocal_z(z);

        if iflag == 1 {
            // Scaled recurrence (Fortran label 120, lines 3760-3788)
            let mut s1 = w[0];
            let mut s2 = w[1];

            let mut l = 2usize; // iteration index (Fortran L=3..NN, we use 2..nn-1)
            while l < nn {
                let ck_val = s2;
                s2 = mul_add_scalar(rz * ck_val, ak_val + fnu, s1);
                s1 = ck_val;
                let ck_scaled = s2 * crscr;
                y[k as usize] = ck_scaled;
                ak_val = ak_val - one;
                k -= 1;
                l += 1;

                if zabs(ck_scaled) > ascle {
                    // Label 140: switch to unscaled recurrence (Fortran lines 3785-3788)
                    // IB = L + 1 in Fortran; remaining iterations continue unscaled
                    while l < nn {
                        let ki = k as usize;
                        y[ki] = mul_add_scalar(rz * y[ki + 1], ak_val + fnu, y[ki + 2]);
                        ak_val = ak_val - one;
                        k -= 1;
                        l += 1;
                    }
                    return nz;
                }
            }
            return nz;
        }

        // Unscaled recurrence (Fortran label 100, lines 3749-3756)
        for _ in 2..nn {
            let ki = k as usize;
            y[ki] = mul_add_scalar(rz * y[ki + 1], ak_val + fnu, y[ki + 2]);
            ak_val = ak_val - one;
            k -= 1;
        }
        return nz;
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
        let mut y = [Complex64::new(0.0, 0.0)];
        let nz = zseri(
            z,
            0.0,
            Scaling::Unscaled,
            &mut y,
            f64::tol(),
            f64::elim(),
            f64::alim(),
        );
        assert_eq!(y[0].re, 1.0);
        assert_eq!(y[0].im, 0.0);
        assert_eq!(nz, 0);
    }
}
