//! Region 2 uniform asymptotic expansion for the K function + analytic continuation.
//!
//! Translation of Fortran ZUNK2 from TOMS 644 / SLATEC (zbsubs.f lines 6159-6664).
//! Computes K(fnu,z) and its analytic continuation from the right half plane
//! to the left half plane using the uniform asymptotic expansion for H(kind,fnu,zn)
//! and J(fnu,zn).

#![allow(clippy::too_many_arguments)]
#![allow(clippy::excessive_precision)]

use num_complex::Complex;

use crate::airy::zairy;
use crate::algo::constants::{AIC, HPI, PI};
use crate::algo::s1s2::zs1s2;
use crate::algo::uchk::zuchk;
use crate::algo::unhj::zunhj;
use crate::machine::BesselFloat;
use crate::types::{Accuracy, AiryDerivative, Scaling, SumOption};
use crate::utils::{mul_neg_i, reciprocal_z, zabs};

/// CR1 = (1, sqrt(3)) (Fortran line 6196-6197)
const CR1: [f64; 2] = [1.0, 1.73205080756887729];

/// CR2 = (-0.5, -sqrt(3)/2) (Fortran line 6197)
const CR2: [f64; 2] = [-0.5, -8.66025403784438647e-01];

/// CIP rotation table (Fortran lines 6201-6203)
const CIPR: [f64; 4] = [1.0, 0.0, -1.0, 0.0];
const CIPI: [f64; 4] = [0.0, -1.0, 0.0, 1.0];

/// Compute K(fnu,z) via Region 2 uniform asymptotic expansion + analytic continuation.
///
/// Equivalent to Fortran ZUNK2 in TOMS 644 (zbsubs.f lines 6159-6664).
///
/// # Parameters
/// - `z`: complex argument
/// - `fnu`: starting order ν >= 0
/// - `kode`: scaling mode
/// - `mr`: analytic continuation direction (0 = none, ±1 = left half plane)
/// - `y`: output slice for sequence members (length determines n)
/// - `tol`, `elim`, `alim`: machine-derived thresholds
///
/// # Returns
/// `nz` where nz = -1 indicates overflow.
pub(crate) fn zunk2<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    mr: i32,
    y: &mut [Complex<T>],
    tol: T,
    elim: T,
    alim: T,
) -> i32 {
    let n = y.len();
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);
    let pi_t = T::from_f64(PI);
    let hpi_t = T::from_f64(HPI);
    let aic_t = T::from_f64(AIC);

    let cr1 = Complex::new(T::from_f64(CR1[0]), T::from_f64(CR1[1]));
    let cr2 = Complex::new(T::from_f64(CR2[0]), T::from_f64(CR2[1]));

    y.fill(czero);
    let mut kdflg: usize = 1;
    let mut nz: i32 = 0;

    // ── 3-level scaling (Fortran lines 6211-6221) ──
    let cscl = one / tol;
    let crsc = tol;
    let cssr = [cscl, one, crsc];
    let csrr = [crsc, one, cscl];
    let bry = [
        T::from_f64(1.0e3) * T::MACH_TINY / tol,
        one / (T::from_f64(1.0e3) * T::MACH_TINY / tol),
        T::MACH_HUGE,
    ];

    // ── Reflect to right half plane (Fortran lines 6222-6226) ──
    let zr = if z.re >= zero { z } else { -z };

    // (Fortran lines 6228-6232)
    let yy = zr.im;
    let zn = Complex::new(zr.im.abs(), -zr.re);
    let zb = Complex::new(zr.re, zr.im.abs());

    // Safety: fnu is finite and < ~1e15 per upper-interface checks
    let inu = fnu.to_i32().unwrap() as usize;
    let fnf = fnu - T::from_f64(inu as f64);
    let ang = -hpi_t * fnf;
    let car = ang.cos();
    let sar = ang.sin();
    let c2r_init = hpi_t * sar;
    let c2i_init = -hpi_t * car;

    // KK = MOD(INU,4) + 1 → 0-based: kk = inu % 4
    let kk = inu % 4;
    let c2_init = Complex::new(c2r_init, c2i_init);
    let cip_kk = Complex::new(T::from_f64(CIPR[kk]), T::from_f64(CIPI[kk]));
    let mut cs = cr1 * c2_init * cip_kk;

    // ── Phase 1: K function computation (Fortran lines 6254-6352) ──
    let mut j: usize = 2; // Fortran J starts at 2
    let mut phi_arr = [czero; 2];
    let mut arg_arr = [czero; 2];
    let mut zeta1_arr = [czero; 2];
    let mut zeta2_arr = [czero; 2];
    let mut asum_arr = [czero; 2];
    let mut bsum_arr = [czero; 2];
    let mut cy = [czero; 2];
    let mut kflag: usize = 2;

    let mut i_exit = n;

    for i in 0..n {
        // J flip-flop (Fortran: J = 3 - J)
        j = 3 - j;
        let jj = j - 1; // 0-based index
        let fn_val = fnu + T::from_f64(i as f64);

        let result = zunhj(zn, fn_val, SumOption::Full, tol);
        phi_arr[jj] = result.phi;
        arg_arr[jj] = result.arg;
        zeta1_arr[jj] = result.zeta1;
        zeta2_arr[jj] = result.zeta2;
        asum_arr[jj] = result.asum;
        bsum_arr[jj] = result.bsum;

        // S1 exponent (Fortran lines 6264-6275)
        let s1_exp = if kode == Scaling::Exponential {
            let st = zb + result.zeta2;
            let rast = fn_val / zabs(st);
            result.zeta1 - st.conj() * (rast * rast)
        } else {
            result.zeta1 - result.zeta2
        };

        let mut rs1 = s1_exp.re;

        // ── Overflow/underflow test (Fortran lines 6280-6351) ──
        if rs1.abs() > elim {
            // label 70
            if rs1 > zero {
                return -1;
            }
            if z.re < zero {
                return -1;
            }
            kdflg = 1;
            y[i] = czero;
            nz += 1;
            cs = mul_neg_i(cs);
            if i > 0 && (y[i - 1] != czero) {
                y[i - 1] = czero;
                nz += 1;
            }
            continue;
        }

        if kdflg == 1 {
            kflag = 2;
        }
        if rs1.abs() >= alim {
            // Refine test (Fortran lines 6287-6293)
            let aphi = zabs(result.phi);
            let aarg = zabs(result.arg);
            rs1 = rs1 + aphi.ln() - T::from_f64(0.25) * aarg.ln() - aic_t;
            if rs1.abs() > elim {
                if rs1 > zero {
                    return -1;
                }
                if z.re < zero {
                    return -1;
                }
                kdflg = 1;
                y[i] = czero;
                nz += 1;
                cs = mul_neg_i(cs);
                if i > 0 && (y[i - 1] != czero) {
                    y[i - 1] = czero;
                    nz += 1;
                }
                continue;
            }
            if kdflg == 1 {
                kflag = 1;
            }
            if rs1 >= zero && kdflg == 1 {
                kflag = 3;
            }
        }

        // ── Compute S2 = PHI*(Ai*ASUM + CR2*(Ai'*BSUM)) * CS (Fortran lines 6299-6318) ──
        // ARG is multiplied by CR2 before ZAIRY call
        let c2_arg = arg_arr[jj] * cr2;

        let (ai, _nai, _) = zairy(c2_arg, AiryDerivative::Value, Scaling::Exponential).unwrap_or((
            czero,
            0,
            Accuracy::Normal,
        ));
        let (dai, _ndai, _) = zairy(c2_arg, AiryDerivative::Derivative, Scaling::Exponential)
            .unwrap_or((czero, 0, Accuracy::Normal));

        let mut s2 = phi_arr[jj] * (ai * asum_arr[jj] + cr2 * (dai * bsum_arr[jj])) * cs;

        // Scale by exp(S1) (Fortran lines 6313-6318)
        let s1_scaled = s1_exp.exp() * cssr[kflag - 1];
        s2 = s2 * s1_scaled;

        if kflag == 1 && zuchk(s2, bry[0], tol) {
            // label 70: underflow
            if rs1 > zero {
                return -1;
            }
            if z.re < zero {
                return -1;
            }
            kdflg = 1;
            y[i] = czero;
            nz += 1;
            cs = mul_neg_i(cs);
            if i > 0 && (y[i - 1] != czero) {
                y[i - 1] = czero;
                nz += 1;
            }
            continue;
        }

        // label 60: store result
        if yy <= zero {
            s2 = s2.conj();
        }
        cy[kdflg - 1] = s2;
        y[i] = s2 * csrr[kflag - 1];

        // CS rotation: CS *= -i (Fortran lines 6328-6330)
        cs = mul_neg_i(cs);

        if kdflg == 2 {
            i_exit = i + 1;
            break;
        }
        kdflg = 2;
    }

    // ── Phase 2: forward recurrence fill (Fortran lines 6354-6443) ──
    if i_exit == n {
        // Loop completed without early exit from kdflg==2
    }

    let fn_at_exit = fnu + T::from_f64((i_exit - 1) as f64);
    let rz = reciprocal_z(zr);
    let mut ck = rz * fn_at_exit;

    let ib = i_exit + 1;

    if ib <= n {
        // Test last member (Fortran lines 6368-6406)
        let fn_last = fnu + T::from_f64((n - 1) as f64);
        let ipard = if mr != 0 {
            SumOption::Full
        } else {
            SumOption::SkipSum
        };
        let result_last = zunhj(zn, fn_last, ipard, tol);

        let s1r_last = if kode == Scaling::Exponential {
            let st = zb + result_last.zeta2;
            let rast = fn_last / zabs(st);
            result_last.zeta1.re - st.re * rast * rast
        } else {
            result_last.zeta1.re - result_last.zeta2.re
        };

        let mut rs1 = s1r_last;
        if rs1.abs() > elim {
            if rs1 > zero {
                return -1;
            }
            if z.re < zero {
                return -1;
            }
            nz = n as i32;
            y.fill(czero);
            return nz;
        }
        if rs1.abs() >= alim {
            let aphi = zabs(result_last.phi);
            rs1 = rs1 + aphi.ln();
            if rs1.abs() >= elim {
                if rs1 > zero {
                    return -1;
                }
                if z.re < zero {
                    return -1;
                }
                nz = n as i32;
                y.fill(czero);
                return nz;
            }
        }

        // Label 120: forward recurrence (Fortran lines 6407-6443)
        let mut s1 = cy[0];
        let mut s2 = cy[1];
        let mut c1r_rec = csrr[kflag - 1];
        let mut ascle = bry[kflag - 1];

        for y_item in y.iter_mut().take(n).skip(ib - 1) {
            let prev = s2;
            s2 = ck * prev + s1;
            s1 = prev;
            ck = ck + rz;
            let c2_scaled = s2 * c1r_rec;
            *y_item = c2_scaled;

            if kflag >= 3 {
                continue;
            }
            let c2m = c2_scaled.re.abs().max(c2_scaled.im.abs());
            if c2m <= ascle {
                continue;
            }
            kflag += 1;
            ascle = bry[kflag - 1];
            s1 = s1 * c1r_rec;
            s2 = c2_scaled;
            s1 = s1 * cssr[kflag - 1];
            s2 = s2 * cssr[kflag - 1];
            c1r_rec = csrr[kflag - 1];
        }
    }

    // ── Phase 3: Analytic continuation (Fortran lines 6444-6660) ──
    if mr == 0 {
        return nz;
    }

    nz = 0;
    let fmr = T::from_f64(mr as f64);
    let sgn = -pi_t.copysign(fmr); // -DSIGN(PI, FMR)

    // CSGN and CSPN (Fortran lines 6455-6463)
    let mut csgni = sgn;
    if yy <= zero {
        csgni = -csgni;
    }
    let ifn = inu + n - 1;
    let ang2 = fnf * sgn;
    let mut cspn = Complex::new(ang2.cos(), ang2.sin());
    if ifn % 2 != 0 {
        cspn = -cspn;
    }

    // CS for I function (Fortran lines 6471-6478)
    let in_idx = ifn % 4;
    let cip_ac = Complex::new(T::from_f64(CIPR[in_idx]), T::from_f64(CIPI[in_idx]));
    let mut cs_ac = Complex::new(sar * csgni, car * csgni) * cip_ac.conj();

    let asc = bry[0];
    let mut iuf: i32 = 0;
    let mut kk_idx = n; // 1-based counting down
    kdflg = 1;
    let ib_ac = ib.saturating_sub(1);
    let ic = ib_ac.saturating_sub(1);
    let mut iflag: usize = 2;

    // Phase 3 loop (Fortran DO 290 K=1,N)
    let mut k_exit = n;
    for k in 0..n {
        let fn_val = fnu + T::from_f64((kk_idx - 1) as f64);

        // ── Logic to reuse phase 1 params or recompute (Fortran labels 172-175) ──
        let (phid, argd, zet1d, zet2d, asumd, bsumd);

        if n <= 2 {
            // label 172: use cached J-slot params
            let jj = j - 1; // current j slot (0-based)
            phid = phi_arr[jj];
            argd = arg_arr[jj];
            zet1d = zeta1_arr[jj];
            zet2d = zeta2_arr[jj];
            asumd = asum_arr[jj];
            bsumd = bsum_arr[jj];
            j = 3 - j; // flip j
        } else if kk_idx == n && ib_ac < n {
            // label 210: fresh computation
            let result_ac = zunhj(zn, fn_val, SumOption::Full, tol);
            phid = result_ac.phi;
            argd = result_ac.arg;
            zet1d = result_ac.zeta1;
            zet2d = result_ac.zeta2;
            asumd = result_ac.asum;
            bsumd = result_ac.bsum;
        } else if kk_idx == ib_ac || kk_idx == ic {
            // label 172: reuse cached
            let jj = j - 1;
            phid = phi_arr[jj];
            argd = arg_arr[jj];
            zet1d = zeta1_arr[jj];
            zet2d = zeta2_arr[jj];
            asumd = asum_arr[jj];
            bsumd = bsum_arr[jj];
            j = 3 - j;
        } else {
            // Fresh ZUNHJ
            let result_ac = zunhj(zn, fn_val, SumOption::Full, tol);
            phid = result_ac.phi;
            argd = result_ac.arg;
            zet1d = result_ac.zeta1;
            zet2d = result_ac.zeta2;
            asumd = result_ac.asum;
            bsumd = result_ac.bsum;
        }

        // S1 exponent for I function (Fortran lines 6514-6525)
        let s1_exp = if kode == Scaling::Exponential {
            let st = zb + zet2d;
            let rast = fn_val / zabs(st);
            st.conj() * (rast * rast) - zet1d
        } else {
            zet2d - zet1d
        };

        let mut rs1 = s1_exp.re;

        // ── Overflow/underflow test (Fortran lines 6530-6602) ──
        if rs1.abs() > elim {
            // label 280
            if rs1 > zero {
                return -1;
            }
            // S2 = 0
            let c2_save = czero;
            let s2_scaled = czero;

            // label 250 continuation
            if yy <= zero {
                // s2i negation not needed since s2 is zero
            }
            cy[kdflg - 1] = czero;
            let s2_k = s2_scaled;

            // Add I and K functions (Fortran lines 6577-6584)
            let mut s1_k = y[kk_idx - 1];
            let mut s2_k_final = s2_k;
            if kode == Scaling::Exponential {
                let s1s2_result = zs1s2(zr, s1_k, s2_k_final, asc, alim, iuf);
                s1_k = s1s2_result.s1;
                s2_k_final = s1s2_result.s2;
                nz += s1s2_result.nz;
                iuf = s1s2_result.iuf;
            }
            y[kk_idx - 1] = cspn * s1_k + s2_k_final;
            kk_idx -= 1;
            cspn = -cspn;
            cs_ac = mul_neg_i(cs_ac);

            if c2_save == czero {
                kdflg = 1;
            } else if kdflg == 2 {
                k_exit = k + 1;
                break;
            } else {
                kdflg = 2;
            }
            continue;
        }

        if kdflg == 1 {
            iflag = 2;
        }
        if rs1.abs() >= alim {
            let aphi = zabs(phid);
            let aarg = zabs(argd);
            rs1 = rs1 + aphi.ln() - T::from_f64(0.25) * aarg.ln() - aic_t;
            if rs1.abs() > elim {
                if rs1 > zero {
                    return -1;
                }
                // S2 = 0 path (goto 250)
                let c2_save = czero;
                cy[kdflg - 1] = czero;
                let s2_scaled = czero;

                let mut s1_k = y[kk_idx - 1];
                let mut s2_k_final = s2_scaled;
                if kode == Scaling::Exponential {
                    let s1s2_result = zs1s2(zr, s1_k, s2_k_final, asc, alim, iuf);
                    s1_k = s1s2_result.s1;
                    s2_k_final = s1s2_result.s2;
                    nz += s1s2_result.nz;
                    iuf = s1s2_result.iuf;
                }
                y[kk_idx - 1] = cspn * s1_k + s2_k_final;
                kk_idx -= 1;
                cspn = -cspn;
                cs_ac = mul_neg_i(cs_ac);

                if c2_save == czero {
                    kdflg = 1;
                } else if kdflg == 2 {
                    k_exit = k + 1;
                    break;
                } else {
                    kdflg = 2;
                }
                continue;
            }
            if kdflg == 1 {
                iflag = 1;
            }
            if rs1 >= zero && kdflg == 1 {
                iflag = 3;
            }
        }

        // ── Compute I function term: S2 = PHI*(Ai*ASUM + Ai'*BSUM)*CS (Fortran lines 6544-6560) ──
        let (ai, _nai, _) = zairy(argd, AiryDerivative::Value, Scaling::Exponential).unwrap_or((
            czero,
            0,
            Accuracy::Normal,
        ));
        let (dai, _ndai, _) = zairy(argd, AiryDerivative::Derivative, Scaling::Exponential)
            .unwrap_or((czero, 0, Accuracy::Normal));

        let mut s2 = phid * (ai * asumd + dai * bsumd) * cs_ac;

        // Scale by exp(S1) (Fortran lines 6555-6560)
        let s1_scaled = s1_exp.exp() * cssr[iflag - 1];
        s2 = s2 * s1_scaled;

        if iflag == 1 && zuchk(s2, bry[0], tol) {
            s2 = czero;
        }

        // label 250 (Fortran lines 6566-6597)
        if yy <= zero {
            s2 = s2.conj();
        }
        cy[kdflg - 1] = s2;
        let c2_save = s2;
        let s2_scaled = s2 * csrr[iflag - 1];

        // Add I and K functions (Fortran lines 6577-6590)
        let mut s1_k = y[kk_idx - 1];
        let mut s2_k = s2_scaled;
        if kode == Scaling::Exponential {
            let s1s2_result = zs1s2(zr, s1_k, s2_k, asc, alim, iuf);
            s1_k = s1s2_result.s1;
            s2_k = s1s2_result.s2;
            nz += s1s2_result.nz;
            iuf = s1s2_result.iuf;
        }
        y[kk_idx - 1] = cspn * s1_k + s2_k;
        kk_idx -= 1;
        cspn = -cspn;
        cs_ac = mul_neg_i(cs_ac);

        // KDFLG state machine (Fortran lines 6591-6597)
        if c2_save == czero {
            kdflg = 1;
        } else if kdflg == 2 {
            k_exit = k + 1;
            break;
        } else {
            kdflg = 2;
        }
    }

    // ── Backward recurrence for remaining I terms (Fortran lines 6604-6660) ──
    let il = n - k_exit;
    if il == 0 {
        return nz;
    }

    let mut s1 = cy[0];
    let mut s2 = cy[1];
    let mut csr_rec = csrr[iflag - 1];
    let mut ascle = bry[iflag - 1];
    let mut fn_rec = T::from_f64((inu + il) as f64);

    for _i in 0..il {
        let prev = s2;
        let cfn = fn_rec + fnf;
        s2 = rz * prev * cfn + s1;
        s1 = prev;
        fn_rec = fn_rec - one;
        let ck_val = s2 * csr_rec;

        // Add K function (Fortran lines 6632-6639)
        let mut c1_k = y[kk_idx - 1];
        let mut c2_k = ck_val;
        if kode == Scaling::Exponential {
            let s1s2_result = zs1s2(zr, c1_k, c2_k, asc, alim, iuf);
            c1_k = s1s2_result.s1;
            c2_k = s1s2_result.s2;
            nz += s1s2_result.nz;
            iuf = s1s2_result.iuf;
        }
        y[kk_idx - 1] = cspn * c1_k + c2_k;
        kk_idx -= 1;
        cspn = -cspn;

        if iflag >= 3 {
            continue;
        }
        let c2m = ck_val.re.abs().max(ck_val.im.abs());
        if c2m <= ascle {
            continue;
        }
        iflag += 1;
        ascle = bry[iflag - 1];
        s1 = s1 * csr_rec;
        s2 = ck_val;
        s1 = s1 * cssr[iflag - 1];
        s2 = s2 * cssr[iflag - 1];
        csr_rec = csrr[iflag - 1];
    }

    nz
}
