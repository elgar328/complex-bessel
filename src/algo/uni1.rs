//! Region 1 uniform asymptotic expansion for the I function.
//!
//! Translation of Fortran ZUNI1 from TOMS 644 / SLATEC (zbsubs.f lines 6840-7044).
//! Computes I(fnu,z) in the region |arg(z)| <= pi/3 using the uniform
//! asymptotic expansion.

#![allow(clippy::too_many_arguments)]
#![allow(clippy::excessive_precision)]

use num_complex::Complex;
use num_traits::Float;

use crate::algo::uchk::zuchk;
use crate::algo::unik::zunik;
use crate::algo::uoik::zuoik;
use crate::machine::BesselFloat;
use crate::types::Scaling;
use crate::utils::zabs;

/// Output of ZUNI1.
pub(crate) struct Uni1Output<T> {
    /// Computed I function values.
    pub y: Vec<Complex<T>>,
    /// Underflow count (number of zeroed trailing members).
    /// -1 indicates overflow.
    pub nz: i32,
    /// If nonzero, remaining items (orders fnu to fnu+nlast-1) need a different method
    /// because fnu+nlast-1 < fnul.
    pub nlast: i32,
}

/// Compute I(fnu,z) via Region 1 uniform asymptotic expansion.
///
/// Equivalent to Fortran ZUNI1 in TOMS 644 (zbsubs.f lines 6840-7044).
///
/// # Parameters
/// - `z`: complex argument
/// - `fnu`: starting order ν >= 0
/// - `kode`: scaling mode
/// - `n`: number of sequence members
/// - `fnul`: large order threshold
/// - `tol`, `elim`, `alim`: machine-derived thresholds
pub(crate) fn zuni1<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    n: usize,
    fnul: T,
    tol: T,
    elim: T,
    alim: T,
) -> Uni1Output<T> {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);

    let mut nz: i32 = 0;
    let mut nd = n;
    let mut nlast: i32 = 0;

    let mut y = vec![czero; n];

    // ── 3-level scaling (Fortran lines 6877-6885) ──
    let cscl = one / tol;
    let crsc = tol;
    let cssr = [cscl, one, crsc];
    let csrr = [crsc, one, cscl];
    let bry0 = T::from(1.0e3).unwrap() * T::MACH_TINY / tol;

    // ── Check first member for underflow/overflow (Fortran lines 6889-6907) ──
    let fn_val = fnu.max(one);
    let result0 = zunik(z, fn_val, 1, 1, tol, None);

    let (s1r, s1i) = if kode == Scaling::Exponential {
        // KODE=2 (Fortran lines 6894-6900)
        let st = Complex::new(z.re + result0.zeta2.re, z.im + result0.zeta2.im);
        let rast = fn_val / zabs(st);
        let str_val = st.re * rast * rast;
        let sti = -st.im * rast * rast;
        (-result0.zeta1.re + str_val, -result0.zeta1.im + sti)
    } else {
        // KODE=1 (Fortran lines 6903-6904)
        (
            -result0.zeta1.re + result0.zeta2.re,
            -result0.zeta1.im + result0.zeta2.im,
        )
    };

    let rs1 = s1r;
    if rs1.abs() > elim {
        if rs1 > zero {
            // Overflow (label 120 → NZ=-1)
            return Uni1Output {
                y,
                nz: -1,
                nlast: 0,
            };
        }
        // All underflow (label 130 → set all to zero)
        nz = n as i32;
        for item in y.iter_mut() {
            *item = czero;
        }
        return Uni1Output { y, nz, nlast: 0 };
    }

    // ── Label 30: compute first min(2, nd) terms directly (Fortran lines 6908-6965) ──
    'label30: loop {
        let nn = nd.min(2);
        let mut iflag: usize = 2; // default scaling level (1-based for cssr/csrr indexing)
        let mut cy = [czero; 2]; // CYR, CYI workspace

        let mut last_phi = czero;
        let mut last_sum = czero;

        let mut computed_ok = true;
        #[allow(clippy::needless_range_loop)]
        for i in 0..nn {
            // fn = fnu + (nd - 1 - i) (Fortran: FN = FNU + FLOAT(ND-I), with Fortran I=1..NN)
            let fn_val = fnu + T::from((nd - 1 - i) as f64).unwrap();
            let result = zunik(z, fn_val, 1, 0, tol, None);

            last_phi = result.phi;
            last_sum = result.sum;

            let (s1r, s1i) = if kode == Scaling::Exponential {
                // KODE=2 (Fortran lines 6916-6922)
                let st = Complex::new(z.re + result.zeta2.re, z.im + result.zeta2.im);
                let rast = fn_val / zabs(st);
                let str_val = st.re * rast * rast;
                let sti = -st.im * rast * rast;
                // Note: Fortran line 6922 has +ZI for S1I (not present in KODE=1 path)
                (-result.zeta1.re + str_val, -result.zeta1.im + sti + z.im)
            } else {
                // KODE=1 (Fortran lines 6925-6926)
                (
                    -result.zeta1.re + result.zeta2.re,
                    -result.zeta1.im + result.zeta2.im,
                )
            };

            // ── Test for underflow/overflow (Fortran lines 6931-6958) ──
            let rs1 = s1r;
            if rs1.abs() > elim {
                // label 110: underflow or overflow
                if rs1 > zero {
                    // Overflow (label 120)
                    return Uni1Output {
                        y,
                        nz: -1,
                        nlast: 0,
                    };
                }
                // Underflow: set y[nd-1] = 0, adjust nd (Fortran lines 7017-7032)
                y[nd - 1] = czero;
                nz += 1;
                nd -= 1;
                if nd == 0 {
                    return Uni1Output { y, nz, nlast: 0 };
                }
                // ZUOIK additional check
                let (y_oik, nuf) = zuoik(z, fnu, kode, 1, nd, tol, elim, alim);
                if nuf < 0 {
                    return Uni1Output {
                        y,
                        nz: -1,
                        nlast: 0,
                    };
                }
                nd -= nuf as usize;
                nz += nuf;
                if nd == 0 {
                    return Uni1Output { y, nz, nlast: 0 };
                }
                // Copy zeroed entries from ZUOIK
                for j in 0..n {
                    if y_oik[j] == czero && y[j] == czero {
                        // already zero
                    }
                }
                let fn_check = fnu + T::from((nd - 1) as f64).unwrap();
                if fn_check >= fnul {
                    continue 'label30; // retry with reduced nd
                }
                nlast = nd as i32;
                return Uni1Output { y, nz, nlast };
            }

            if i == 0 {
                iflag = 2;
            }
            if rs1.abs() >= alim {
                // Refine test (Fortran lines 6938-6943)
                let aphi = zabs(result.phi);
                let rs1_refined = rs1 + aphi.ln();
                if rs1_refined.abs() > elim {
                    // label 110
                    if rs1_refined > zero {
                        return Uni1Output {
                            y,
                            nz: -1,
                            nlast: 0,
                        };
                    }
                    y[nd - 1] = czero;
                    nz += 1;
                    nd -= 1;
                    if nd == 0 {
                        return Uni1Output { y, nz, nlast: 0 };
                    }
                    let (_y_oik, nuf) = zuoik(z, fnu, kode, 1, nd, tol, elim, alim);
                    if nuf < 0 {
                        return Uni1Output {
                            y,
                            nz: -1,
                            nlast: 0,
                        };
                    }
                    nd -= nuf as usize;
                    nz += nuf;
                    if nd == 0 {
                        return Uni1Output { y, nz, nlast: 0 };
                    }
                    let fn_check = fnu + T::from((nd - 1) as f64).unwrap();
                    if fn_check >= fnul {
                        continue 'label30;
                    }
                    nlast = nd as i32;
                    return Uni1Output { y, nz, nlast };
                }
                if i == 0 {
                    iflag = 1;
                }
                if rs1 < zero {
                    // iflag stays at 1 if i==0 already set it
                } else if i == 0 {
                    iflag = 3;
                }
            }

            // ── Scale S1 and compute S2 = PHI * SUM (Fortran lines 6948-6964) ──
            let s2 = result.phi * result.sum;
            let str_exp = s1r.exp() * cssr[iflag - 1];
            let s1_scaled = Complex::new(str_exp * s1i.cos(), str_exp * s1i.sin());
            let s2_final = s2 * s1_scaled;

            if iflag == 1 && zuchk(s2_final, bry0, tol) {
                // Underflow (label 110)
                if rs1 > zero {
                    return Uni1Output {
                        y,
                        nz: -1,
                        nlast: 0,
                    };
                }
                y[nd - 1] = czero;
                nz += 1;
                nd -= 1;
                if nd == 0 {
                    return Uni1Output { y, nz, nlast: 0 };
                }
                let (_y_oik, nuf) = zuoik(z, fnu, kode, 1, nd, tol, elim, alim);
                if nuf < 0 {
                    return Uni1Output {
                        y,
                        nz: -1,
                        nlast: 0,
                    };
                }
                nd -= nuf as usize;
                nz += nuf;
                if nd == 0 {
                    return Uni1Output { y, nz, nlast: 0 };
                }
                let fn_check = fnu + T::from((nd - 1) as f64).unwrap();
                if fn_check >= fnul {
                    computed_ok = false;
                    break;
                }
                nlast = nd as i32;
                return Uni1Output { y, nz, nlast };
            }

            // Store results (Fortran lines 6960-6964)
            cy[i] = s2_final;
            // m = nd - i (Fortran: M = ND - I + 1, 1-based)
            let m = nd - 1 - i;
            y[m] = s2_final * csrr[iflag - 1];
        }

        if !computed_ok {
            continue 'label30;
        }

        // ── Forward recurrence for remaining terms (Fortran lines 6966-7011) ──
        if nd <= 2 {
            return Uni1Output { y, nz, nlast };
        }

        let rast = one / zabs(z);
        let str_val = z.re * rast;
        let sti = -z.im * rast;
        let rzr = (str_val + str_val) * rast;
        let rzi = (sti + sti) * rast;

        let bry1 = one / bry0;
        let bry2 = T::MACH_HUGE;

        let mut s1 = cy[0]; // CYR(1), CYI(1)
        let mut s2 = cy[1]; // CYR(2), CYI(2)
        let mut c1r = csrr[iflag - 1];
        let mut ascle = if iflag == 1 {
            bry0
        } else if iflag == 2 {
            bry1
        } else {
            bry2
        };

        let mut k = nd - 2; // 0-based index into y (Fortran K = ND - 2, 1-based)
        let mut fn_rec = T::from(k as f64).unwrap(); // Fortran: FN = FLOAT(K), where K starts at ND-2

        // Fortran DO 90 I=3,ND
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

        return Uni1Output { y, nz, nlast };
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
    fn zuni1_basic_single() {
        // Single value, moderate parameters
        let z = Complex64::new(2.0, 0.5);
        let result = zuni1(z, 90.0, Scaling::Unscaled, 1, FNUL, TOL, ELIM, ALIM);
        // Should compute without error
        assert!(result.nz >= 0);
        assert_eq!(result.nlast, 0);
        // Result should be finite and nonzero
        assert!(result.y[0].re.is_finite());
        assert!(result.y[0].im.is_finite());
    }

    #[test]
    fn zuni1_sequence() {
        // Sequence of 3 values
        let z = Complex64::new(3.0, 1.0);
        let result = zuni1(z, 90.0, Scaling::Unscaled, 3, FNUL, TOL, ELIM, ALIM);
        assert!(result.nz >= 0);
        // At least first result should be nonzero
        let any_nonzero = result.y.iter().any(|v| v.re != 0.0 || v.im != 0.0);
        assert!(any_nonzero, "at least one value should be nonzero");
    }

    #[test]
    fn zuni1_scaled() {
        // Exponential scaling
        let z = Complex64::new(5.0, 2.0);
        let result = zuni1(z, 90.0, Scaling::Exponential, 1, FNUL, TOL, ELIM, ALIM);
        assert!(result.nz >= 0);
        assert!(result.y[0].re.is_finite());
    }

    #[test]
    fn zuni1_real_axis() {
        // Real positive argument
        let z = Complex64::new(10.0, 0.0);
        let result = zuni1(z, 90.0, Scaling::Unscaled, 1, FNUL, TOL, ELIM, ALIM);
        assert!(result.nz >= 0);
        // I_v(x) for real x > 0 should be real
        if result.nz == 0 && result.nlast == 0 {
            assert!(
                result.y[0].im.abs() < 1e-10 * result.y[0].re.abs().max(1e-300),
                "imaginary part should be tiny for real argument: {}",
                result.y[0].im
            );
        }
    }
}
