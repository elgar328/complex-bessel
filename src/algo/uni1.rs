//! Region 1 uniform asymptotic expansion for the I function.
//!
//! Translation of Fortran ZUNI1 from TOMS 644 / SLATEC (zbsubs.f lines 6840-7044).
//! Computes I(fnu,z) in the region |arg(z)| <= pi/3 using the uniform
//! asymptotic expansion.

#![allow(clippy::too_many_arguments)]

use num_complex::Complex;

use crate::algo::uchk::zuchk;
use crate::algo::unik::zunik;
use crate::algo::uoik::zuoik;
use crate::machine::BesselFloat;
use crate::types::{IkFlag, Scaling, SumOption};
use crate::utils::{reciprocal_z, zabs};

/// Output of ZUNI1.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Uni1Output {
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
/// - `y`: output slice for sequence members (length determines n)
/// - `fnul`: large order threshold
/// - `tol`, `elim`, `alim`: machine-derived thresholds
pub(crate) fn zuni1<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    y: &mut [Complex<T>],
    fnul: T,
    tol: T,
    elim: T,
    alim: T,
) -> Uni1Output {
    let n = y.len();

    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);

    let mut nz: i32 = 0;
    let mut nd = n;
    let mut nlast: i32 = 0;

    // Initialize output slice to zero
    y.fill(czero);

    // ── 3-level scaling (Fortran lines 6877-6885) ──
    let cscl = one / tol;
    let crsc = tol;
    let cssr = [cscl, one, crsc];
    let csrr = [crsc, one, cscl];
    let bry0 = T::from_f64(1.0e3) * T::MACH_TINY / tol;

    // ── Check first member for underflow/overflow (Fortran lines 6889-6907) ──
    let fn_val = fnu.max(one);
    let result0 = zunik(z, fn_val, IkFlag::I, SumOption::SkipSum, tol, None);

    let s1_check = if kode == Scaling::Exponential {
        // KODE=2 (Fortran lines 6894-6900)
        let st = z + result0.zeta2;
        let rast = fn_val / zabs(st);
        st.conj() * (rast * rast) - result0.zeta1
    } else {
        // KODE=1 (Fortran lines 6903-6904)
        result0.zeta2 - result0.zeta1
    };

    let rs1 = s1_check.re;
    if rs1.abs() > elim {
        if rs1 > zero {
            // Overflow (label 120 → NZ=-1)
            return Uni1Output { nz: -1, nlast: 0 };
        }
        // All underflow (label 130 → set all to zero)
        nz = n as i32;
        y.fill(czero);
        return Uni1Output { nz, nlast: 0 };
    }

    // ── Label 30: compute first min(2, nd) terms directly (Fortran lines 6908-6965) ──
    'label30: loop {
        let nn = nd.min(2);
        let mut iflag: usize = 2; // default scaling level (1-based for cssr/csrr indexing)
        let mut cy = [czero; 2]; // CYR, CYI workspace

        let mut computed_ok = true;
        #[allow(clippy::needless_range_loop)]
        for i in 0..nn {
            // fn = fnu + (nd - 1 - i) (Fortran: FN = FNU + FLOAT(ND-I), with Fortran I=1..NN)
            let fn_val = fnu + T::from_f64((nd - 1 - i) as f64);
            let result = zunik(z, fn_val, IkFlag::I, SumOption::Full, tol, None);

            let s1_exp = if kode == Scaling::Exponential {
                // KODE=2 (Fortran lines 6916-6922)
                let st = z + result.zeta2;
                let rast = fn_val / zabs(st);
                // Note: Fortran line 6922 has +ZI for S1I (not present in KODE=1 path)
                let mut s = st.conj() * (rast * rast) - result.zeta1;
                s.im = s.im + z.im;
                s
            } else {
                // KODE=1 (Fortran lines 6925-6926)
                result.zeta2 - result.zeta1
            };

            // ── Test for underflow/overflow (Fortran lines 6931-6958) ──
            let rs1 = s1_exp.re;
            if rs1.abs() > elim {
                // label 110: underflow or overflow
                if rs1 > zero {
                    // Overflow (label 120)
                    return Uni1Output { nz: -1, nlast: 0 };
                }
                // Underflow: set y[nd-1] = 0, adjust nd (Fortran lines 7017-7032)
                y[nd - 1] = czero;
                nz += 1;
                nd -= 1;
                if nd == 0 {
                    return Uni1Output { nz, nlast: 0 };
                }
                // ZUOIK additional check
                let nuf = zuoik(z, fnu, kode, IkFlag::I, &mut y[..nd], tol, elim, alim);
                if nuf < 0 {
                    return Uni1Output { nz: -1, nlast: 0 };
                }
                nd -= nuf as usize;
                nz += nuf;
                if nd == 0 {
                    return Uni1Output { nz, nlast: 0 };
                }
                let fn_check = fnu + T::from_f64((nd - 1) as f64);
                if fn_check >= fnul {
                    continue 'label30; // retry with reduced nd
                }
                nlast = nd as i32;
                return Uni1Output { nz, nlast };
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
                        return Uni1Output { nz: -1, nlast: 0 };
                    }
                    y[nd - 1] = czero;
                    nz += 1;
                    nd -= 1;
                    if nd == 0 {
                        return Uni1Output { nz, nlast: 0 };
                    }
                    let nuf = zuoik(z, fnu, kode, IkFlag::I, &mut y[..nd], tol, elim, alim);
                    if nuf < 0 {
                        return Uni1Output { nz: -1, nlast: 0 };
                    }
                    nd -= nuf as usize;
                    nz += nuf;
                    if nd == 0 {
                        return Uni1Output { nz, nlast: 0 };
                    }
                    let fn_check = fnu + T::from_f64((nd - 1) as f64);
                    if fn_check >= fnul {
                        continue 'label30;
                    }
                    nlast = nd as i32;
                    return Uni1Output { nz, nlast };
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
            let s1_scaled = s1_exp.exp() * cssr[iflag - 1];
            let s2_final = s2 * s1_scaled;

            if iflag == 1 && zuchk(s2_final, bry0, tol) {
                // Underflow (label 110)
                if rs1 > zero {
                    return Uni1Output { nz: -1, nlast: 0 };
                }
                y[nd - 1] = czero;
                nz += 1;
                nd -= 1;
                if nd == 0 {
                    return Uni1Output { nz, nlast: 0 };
                }
                let nuf = zuoik(z, fnu, kode, IkFlag::I, &mut y[..nd], tol, elim, alim);
                if nuf < 0 {
                    return Uni1Output { nz: -1, nlast: 0 };
                }
                nd -= nuf as usize;
                nz += nuf;
                if nd == 0 {
                    return Uni1Output { nz, nlast: 0 };
                }
                let fn_check = fnu + T::from_f64((nd - 1) as f64);
                if fn_check >= fnul {
                    computed_ok = false;
                    break;
                }
                nlast = nd as i32;
                return Uni1Output { nz, nlast };
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
            return Uni1Output { nz, nlast };
        }

        let rz = reciprocal_z(z);

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
        let mut fn_rec = T::from_f64(k as f64); // Fortran: FN = FLOAT(K), where K starts at ND-2

        // Fortran DO 90 I=3,ND
        for _i in 2..nd {
            let c2 = s2;
            let cfn = fnu + fn_rec;
            s2 = s1 + rz * c2 * cfn;
            s1 = c2;
            let c2_scaled = s2 * c1r;
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
            s1 = s1 * c1r;
            s2 = c2_scaled;
            s1 = s1 * cssr[iflag - 1];
            s2 = s2 * cssr[iflag - 1];
            c1r = csrr[iflag - 1];
        }

        return Uni1Output { nz, nlast };
    }
}
