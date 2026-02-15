//! Region 2 uniform asymptotic expansion for the I function.
//!
//! Translation of Fortran ZUNI2 from TOMS 644 / SLATEC (zbsubs.f lines 7045-7312).
//! Computes I(fnu,z) in the region |arg(z)| > pi/3 using the uniform
//! asymptotic expansion for J(fnu, z*exp(m*hpi*i)).

#![allow(clippy::too_many_arguments)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]

use num_complex::Complex;

use crate::airy::zairy;
use crate::algo::uchk::zuchk;
use crate::algo::unhj::zunhj;
use crate::algo::uoik::zuoik;
use crate::machine::BesselFloat;
use crate::types::{AiryDerivative, Scaling};
use crate::utils::zabs;

/// HPI = pi/2 (Fortran line 7079)
const HPI: f64 = 1.57079632679489662e+00;

/// AIC = ln(sqrt(pi/2)) (Fortran line 7080)
const AIC: f64 = 1.265512123484645396e+00;

/// CIP rotation table (Fortran lines 7077-7078)
/// CIPR + i*CIPI: {1, i, -1, -i}
const CIPR: [f64; 4] = [1.0, 0.0, -1.0, 0.0];
const CIPI: [f64; 4] = [0.0, 1.0, 0.0, -1.0];

/// Output of ZUNI2.
pub(crate) struct Uni2Output<T> {
    /// Computed I function values.
    pub y: Vec<Complex<T>>,
    /// Underflow count (number of zeroed trailing members).
    /// -1 indicates overflow.
    pub nz: i32,
    /// If nonzero, remaining items need a different method.
    pub nlast: i32,
}

/// Compute I(fnu,z) via Region 2 uniform asymptotic expansion.
///
/// Equivalent to Fortran ZUNI2 in TOMS 644 (zbsubs.f lines 7045-7312).
pub(crate) fn zuni2<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    n: usize,
    fnul: T,
    tol: T,
    elim: T,
    alim: T,
) -> Uni2Output<T> {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);

    let mut nz: i32 = 0;
    let mut nd = n;
    let mut nlast: i32 = 0;

    let mut y = vec![czero; n];

    // ── 3-level scaling (Fortran lines 7090-7097) ──
    let cscl = one / tol;
    let crsc = tol;
    let cssr = [cscl, one, crsc];
    let csrr = [crsc, one, cscl];
    let bry0 = T::from(1.0e3).unwrap() * T::MACH_TINY / tol;

    // ── Rotate z to right half plane (Fortran lines 7102-7106) ──
    let mut znr = z.im; // ZNR = ZI
    let zni = -z.re; // ZNI = -ZR
    let zbr = z.re;
    let mut zbi = z.im;
    let mut cidi = -one; // CIDI = -1

    let inu = fnu.to_i32().unwrap() as usize;
    let ang = T::from(HPI).unwrap() * (fnu - T::from(inu as f64).unwrap());
    let c2r_init = ang.cos();
    let c2i_init = ang.sin();
    let car = c2r_init;
    let sar = c2i_init;

    // IN = (INU + N - 1) mod 4 (0-based)
    let in_idx = (inu + n - 1) % 4;
    let mut c2r =
        c2r_init * T::from(CIPR[in_idx]).unwrap() - c2i_init * T::from(CIPI[in_idx]).unwrap();
    let mut c2i =
        c2r_init * T::from(CIPI[in_idx]).unwrap() + c2i_init * T::from(CIPR[in_idx]).unwrap();

    if z.im <= zero {
        // Fortran lines 7119-7122
        znr = -znr;
        zbi = -zbi;
        cidi = -cidi;
        c2i = -c2i;
    }

    // ── Check first member for overflow/underflow (Fortran lines 7127-7144) ──
    let fn_val = fnu.max(one);
    let result0 = zunhj(Complex::new(znr, zni), fn_val, 1, tol);

    let rs1 = if kode == Scaling::Exponential {
        let str = zbr + result0.zeta2.re;
        let sti = zbi + result0.zeta2.im;
        let rast = fn_val / zabs(Complex::new(str, sti));
        let str2 = str * rast * rast;
        let _sti2 = -sti * rast * rast;
        -result0.zeta1.re + str2
    } else {
        -result0.zeta1.re + result0.zeta2.re
    };

    if rs1.abs() > elim {
        if rs1 > zero {
            return Uni2Output {
                y,
                nz: -1,
                nlast: 0,
            };
        }
        // All underflow
        nz = n as i32;
        for item in y.iter_mut() {
            *item = czero;
        }
        return Uni2Output { y, nz, nlast: 0 };
    }

    // ── Label 40: main computation loop (Fortran lines 7145-7263) ──
    'label40: loop {
        let nn = nd.min(2);
        let mut iflag: usize = 2;
        let mut cy = [czero; 2];

        let mut computed_ok = true;

        #[allow(clippy::needless_range_loop)]
        for i in 0..nn {
            let fn_val = fnu + T::from((nd - 1 - i) as f64).unwrap();
            let result = zunhj(Complex::new(znr, zni), fn_val, 0, tol);

            let (s1r, s1i) = if kode == Scaling::Exponential {
                // KODE=2 path (Fortran lines 7151-7158)
                // NOTE: s1i += |zi| (absolute value!), NOT zi like ZUNI1
                let str = zbr + result.zeta2.re;
                let sti = zbi + result.zeta2.im;
                let rast = fn_val / zabs(Complex::new(str, sti));
                let str2 = str * rast * rast;
                let sti2 = -sti * rast * rast;
                (
                    -result.zeta1.re + str2,
                    -result.zeta1.im + sti2 + z.im.abs(),
                )
            } else {
                (
                    -result.zeta1.re + result.zeta2.re,
                    -result.zeta1.im + result.zeta2.im,
                )
            };

            // ── Overflow/underflow test (Fortran lines 7167-7203) ──
            let mut rs1 = s1r;
            if rs1.abs() > elim {
                if rs1 > zero {
                    return Uni2Output {
                        y,
                        nz: -1,
                        nlast: 0,
                    };
                }
                // Underflow (Fortran label 120)
                y[nd - 1] = czero;
                nz += 1;
                nd -= 1;
                if nd == 0 {
                    return Uni2Output { y, nz, nlast: 0 };
                }
                let (_y_oik, nuf) = zuoik(z, fnu, kode, 1, nd, tol, elim, alim);
                if nuf < 0 {
                    return Uni2Output {
                        y,
                        nz: -1,
                        nlast: 0,
                    };
                }
                nd -= nuf as usize;
                nz += nuf;
                if nd == 0 {
                    return Uni2Output { y, nz, nlast: 0 };
                }
                let fn_check = fnu + T::from((nd - 1) as f64).unwrap();
                if fn_check >= fnul {
                    // Recalculate C2 (Fortran lines 7292-7296)
                    let in_idx2 = (inu + nd - 1) % 4;
                    c2r = car * T::from(CIPR[in_idx2]).unwrap()
                        - sar * T::from(CIPI[in_idx2]).unwrap();
                    c2i = car * T::from(CIPI[in_idx2]).unwrap()
                        + sar * T::from(CIPR[in_idx2]).unwrap();
                    if z.im <= zero {
                        c2i = -c2i;
                    }
                    computed_ok = false;
                    break;
                }
                nlast = nd as i32;
                return Uni2Output { y, nz, nlast };
            }

            if i == 0 {
                iflag = 2;
            }
            if rs1.abs() >= alim {
                // Refine test (Fortran lines 7175-7181)
                let aphi = zabs(result.phi);
                let aarg = zabs(result.arg);
                rs1 = rs1 + aphi.ln() - T::from(0.25).unwrap() * aarg.ln() - T::from(AIC).unwrap();
                if rs1.abs() > elim {
                    if rs1 > zero {
                        return Uni2Output {
                            y,
                            nz: -1,
                            nlast: 0,
                        };
                    }
                    y[nd - 1] = czero;
                    nz += 1;
                    nd -= 1;
                    if nd == 0 {
                        return Uni2Output { y, nz, nlast: 0 };
                    }
                    let (_y_oik, nuf) = zuoik(z, fnu, kode, 1, nd, tol, elim, alim);
                    if nuf < 0 {
                        return Uni2Output {
                            y,
                            nz: -1,
                            nlast: 0,
                        };
                    }
                    nd -= nuf as usize;
                    nz += nuf;
                    if nd == 0 {
                        return Uni2Output { y, nz, nlast: 0 };
                    }
                    let fn_check = fnu + T::from((nd - 1) as f64).unwrap();
                    if fn_check >= fnul {
                        let in_idx2 = (inu + nd - 1) % 4;
                        c2r = car * T::from(CIPR[in_idx2]).unwrap()
                            - sar * T::from(CIPI[in_idx2]).unwrap();
                        c2i = car * T::from(CIPI[in_idx2]).unwrap()
                            + sar * T::from(CIPR[in_idx2]).unwrap();
                        if z.im <= zero {
                            c2i = -c2i;
                        }
                        computed_ok = false;
                        break;
                    }
                    nlast = nd as i32;
                    return Uni2Output { y, nz, nlast };
                }
                if i == 0 {
                    iflag = 1;
                }
                if rs1 < zero {
                    // iflag stays
                } else if i == 0 {
                    iflag = 3;
                }
            }

            // ── Scale S1 and compute S2 (Fortran lines 7187-7216) ──
            // S2 = PHI * (Ai*ASUM + Ai'*BSUM)
            let (ai, _nai) = zairy(result.arg, AiryDerivative::Value, Scaling::Exponential)
                .unwrap_or((czero, 0));
            let (dai, _ndai) = zairy(result.arg, AiryDerivative::Derivative, Scaling::Exponential)
                .unwrap_or((czero, 0));

            let str_d = dai.re * result.bsum.re - dai.im * result.bsum.im;
            let sti_d = dai.re * result.bsum.im + dai.im * result.bsum.re;
            let str_a = str_d + (ai.re * result.asum.re - ai.im * result.asum.im);
            let sti_a = sti_d + (ai.re * result.asum.im + ai.im * result.asum.re);
            let s2r = result.phi.re * str_a - result.phi.im * sti_a;
            let s2i = result.phi.re * sti_a + result.phi.im * str_a;

            let str_exp = s1r.exp() * cssr[iflag - 1];
            let s1_scaled_re = str_exp * s1i.cos();
            let s1_scaled_im = str_exp * s1i.sin();
            let mut s2_re = s2r * s1_scaled_re - s2i * s1_scaled_im;
            let mut s2_im = s2r * s1_scaled_im + s2i * s1_scaled_re;

            if iflag == 1 && zuchk(Complex::new(s2_re, s2_im), bry0, tol) {
                // Underflow (label 120)
                y[nd - 1] = czero;
                nz += 1;
                nd -= 1;
                if nd == 0 {
                    return Uni2Output { y, nz, nlast: 0 };
                }
                let (_y_oik, nuf) = zuoik(z, fnu, kode, 1, nd, tol, elim, alim);
                if nuf < 0 {
                    return Uni2Output {
                        y,
                        nz: -1,
                        nlast: 0,
                    };
                }
                nd -= nuf as usize;
                nz += nuf;
                if nd == 0 {
                    return Uni2Output { y, nz, nlast: 0 };
                }
                let fn_check = fnu + T::from((nd - 1) as f64).unwrap();
                if fn_check >= fnul {
                    let in_idx2 = (inu + nd - 1) % 4;
                    c2r = car * T::from(CIPR[in_idx2]).unwrap()
                        - sar * T::from(CIPI[in_idx2]).unwrap();
                    c2i = car * T::from(CIPI[in_idx2]).unwrap()
                        + sar * T::from(CIPR[in_idx2]).unwrap();
                    if z.im <= zero {
                        c2i = -c2i;
                    }
                    computed_ok = false;
                    break;
                }
                nlast = nd as i32;
                return Uni2Output { y, nz, nlast };
            }

            // Conjugation for ZI <= 0 (Fortran line 7205)
            if z.im <= zero {
                s2_im = -s2_im;
            }

            // Phase rotation by C2 (Fortran lines 7206-7208)
            let str_c = s2_re * c2r - s2_im * c2i;
            s2_im = s2_re * c2i + s2_im * c2r;
            s2_re = str_c;

            // Store (Fortran lines 7209-7213)
            cy[i] = Complex::new(s2_re, s2_im);
            let j_idx = nd - 1 - i;
            y[j_idx] = Complex::new(s2_re * csrr[iflag - 1], s2_im * csrr[iflag - 1]);

            // C2 update: C2 *= (0, CIDI) → C2R = -C2I*CIDI, C2I = C2R*CIDI
            // (Fortran lines 7214-7216)
            let str2 = -c2i * cidi;
            c2i = c2r * cidi;
            c2r = str2;
        }

        if !computed_ok {
            continue 'label40;
        }

        // ── Forward recurrence for remaining terms (Fortran lines 7218-7263) ──
        if nd <= 2 {
            return Uni2Output { y, nz, nlast };
        }

        let raz = one / zabs(z);
        let str = z.re * raz;
        let sti = -z.im * raz;
        let rzr = (str + str) * raz;
        let rzi = (sti + sti) * raz;

        let bry1 = one / bry0;
        let bry2 = T::MACH_HUGE;

        let mut s1 = cy[0];
        let mut s2 = cy[1];
        let mut c1r = csrr[iflag - 1];
        let mut ascle = if iflag == 1 {
            bry0
        } else if iflag == 2 {
            bry1
        } else {
            bry2
        };

        let mut k = nd - 2;
        let mut fn_rec = T::from(k as f64).unwrap();

        for _i in 2..nd {
            let c2 = s2;
            let cfn = fnu + fn_rec;
            s2 = s1
                + Complex::new(
                    cfn * (rzr * c2.re - rzi * c2.im),
                    cfn * (rzr * c2.im + rzi * c2.re),
                );
            s1 = c2;
            let c2_scaled = Complex::new(s2.re * c1r, s2.im * c1r);
            y[k] = c2_scaled;
            k -= 1;
            fn_rec = fn_rec - one;

            if iflag >= 3 {
                continue;
            }
            let c2m = c2_scaled.re.abs().max(c2_scaled.im.abs());
            if c2m <= ascle {
                continue;
            }
            iflag += 1;
            ascle = if iflag == 2 { bry1 } else { bry2 };
            s1 = Complex::new(s1.re * c1r, s1.im * c1r);
            s2 = c2_scaled;
            s1 = Complex::new(s1.re * cssr[iflag - 1], s1.im * cssr[iflag - 1]);
            s2 = Complex::new(s2.re * cssr[iflag - 1], s2.im * cssr[iflag - 1]);
            c1r = csrr[iflag - 1];
        }

        return Uni2Output { y, nz, nlast };
    }
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
    fn zuni2_basic_single() {
        // Region 2 argument: |Im(z)| > |Re(z)|*sqrt(3)
        let z = Complex64::new(1.0, 10.0);
        let result = zuni2(z, 90.0, Scaling::Unscaled, 1, FNUL, TOL, ELIM, ALIM);
        assert!(result.nz >= 0, "nz={}", result.nz);
        assert_eq!(result.nlast, 0);
        assert!(result.y[0].re.is_finite());
    }

    #[test]
    fn zuni2_sequence() {
        let z = Complex64::new(2.0, 15.0);
        let result = zuni2(z, 90.0, Scaling::Unscaled, 3, FNUL, TOL, ELIM, ALIM);
        assert!(result.nz >= 0);
    }

    #[test]
    fn zuni2_scaled() {
        let z = Complex64::new(1.0, 10.0);
        let result = zuni2(z, 90.0, Scaling::Exponential, 1, FNUL, TOL, ELIM, ALIM);
        assert!(result.nz >= 0);
        assert!(result.y[0].re.is_finite());
    }

    #[test]
    fn zuni2_negative_imaginary() {
        // z.im < 0 triggers conjugation path
        let z = Complex64::new(1.0, -10.0);
        let result = zuni2(z, 90.0, Scaling::Unscaled, 1, FNUL, TOL, ELIM, ALIM);
        assert!(result.nz >= 0);
    }
}
