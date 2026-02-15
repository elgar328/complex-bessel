//! Region 2 uniform asymptotic expansion for the K function + analytic continuation.
//!
//! Translation of Fortran ZUNK2 from TOMS 644 / SLATEC (zbsubs.f lines 6159-6664).
//! Computes K(fnu,z) and its analytic continuation from the right half plane
//! to the left half plane using the uniform asymptotic expansion for H(kind,fnu,zn)
//! and J(fnu,zn).

#![allow(clippy::too_many_arguments)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]

use num_complex::Complex;

use crate::airy::zairy;
use crate::algo::s1s2::zs1s2;
use crate::algo::uchk::zuchk;
use crate::algo::unhj::zunhj;
use crate::machine::BesselFloat;
use crate::types::{AiryDerivative, Scaling};
use crate::utils::zabs;

const HPI: f64 = 1.57079632679489662e+00;
const PI: f64 = 3.14159265358979324e+00;
const AIC: f64 = 1.26551212348464539e+00; // ln(sqrt(pi/2))

/// CR1 = (1, sqrt(3)) (Fortran line 6196-6197)
const CR1R: f64 = 1.0;
const CR1I: f64 = 1.73205080756887729;

/// CR2 = (-0.5, -sqrt(3)/2) (Fortran line 6197)
const CR2R: f64 = -0.5;
const CR2I: f64 = -8.66025403784438647e-01;

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
/// - `n`: number of sequence members
/// - `tol`, `elim`, `alim`: machine-derived thresholds
///
/// # Returns
/// `(y, nz)` where nz = -1 indicates overflow.
pub(crate) fn zunk2<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    mr: i32,
    n: usize,
    tol: T,
    elim: T,
    alim: T,
) -> (Vec<Complex<T>>, i32) {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);
    let pi_t = T::from(PI).unwrap();
    let hpi_t = T::from(HPI).unwrap();
    let aic_t = T::from(AIC).unwrap();

    let cr1r = T::from(CR1R).unwrap();
    let cr1i = T::from(CR1I).unwrap();
    let cr2r = T::from(CR2R).unwrap();
    let cr2i = T::from(CR2I).unwrap();

    let mut y = vec![czero; n];
    let mut kdflg: usize = 1;
    let mut nz: i32 = 0;

    // ── 3-level scaling (Fortran lines 6211-6221) ──
    let cscl = one / tol;
    let crsc = tol;
    let cssr = [cscl, one, crsc];
    let csrr = [crsc, one, cscl];
    let bry = [
        T::from(1.0e3).unwrap() * T::MACH_TINY / tol,
        one / (T::from(1.0e3).unwrap() * T::MACH_TINY / tol),
        T::MACH_HUGE,
    ];

    // ── Reflect to right half plane (Fortran lines 6222-6226) ──
    let (zrr, zri) = if z.re >= zero {
        (z.re, z.im)
    } else {
        (-z.re, -z.im)
    };

    // (Fortran lines 6228-6232)
    let yy = zri;
    let mut znr = zri;
    let zni = -zrr;
    let zbr = zrr;
    let mut zbi = zri;

    let inu = fnu.to_i32().unwrap() as usize;
    let fnf = fnu - T::from(inu as f64).unwrap();
    let ang = -hpi_t * fnf;
    let car = ang.cos();
    let sar = ang.sin();
    let c2r_init = hpi_t * sar;
    let c2i_init = -hpi_t * car;

    // KK = MOD(INU,4) + 1 → 0-based: kk = inu % 4
    let kk = inu % 4;
    let str0 = c2r_init * T::from(CIPR[kk]).unwrap() - c2i_init * T::from(CIPI[kk]).unwrap();
    let sti0 = c2r_init * T::from(CIPI[kk]).unwrap() + c2i_init * T::from(CIPR[kk]).unwrap();
    let mut csr = cr1r * str0 - cr1i * sti0;
    let mut csi = cr1r * sti0 + cr1i * str0;

    if yy <= zero {
        znr = -znr;
        zbi = -zbi;
    }

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
        let fn_val = fnu + T::from(i as f64).unwrap();

        let result = zunhj(Complex::new(znr, zni), fn_val, 0, tol);
        phi_arr[jj] = result.phi;
        arg_arr[jj] = result.arg;
        zeta1_arr[jj] = result.zeta1;
        zeta2_arr[jj] = result.zeta2;
        asum_arr[jj] = result.asum;
        bsum_arr[jj] = result.bsum;

        // S1 exponent (Fortran lines 6264-6275)
        let (s1r, s1i) = if kode == Scaling::Exponential {
            let str = zbr + result.zeta2.re;
            let sti = zbi + result.zeta2.im;
            let rast = fn_val / zabs(Complex::new(str, sti));
            let str2 = str * rast * rast;
            let sti2 = -sti * rast * rast;
            (result.zeta1.re - str2, result.zeta1.im - sti2)
        } else {
            (
                result.zeta1.re - result.zeta2.re,
                result.zeta1.im - result.zeta2.im,
            )
        };

        let mut rs1 = s1r;

        // ── Overflow/underflow test (Fortran lines 6280-6351) ──
        if rs1.abs() > elim {
            // label 70
            if rs1 > zero {
                return (y, -1);
            }
            if z.re < zero {
                return (y, -1);
            }
            kdflg = 1;
            y[i] = czero;
            nz += 1;
            let str_cs = csi;
            csi = -csr;
            csr = str_cs;
            if i > 0 && (y[i - 1].re != zero || y[i - 1].im != zero) {
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
            rs1 = rs1 + aphi.ln() - T::from(0.25).unwrap() * aarg.ln() - aic_t;
            if rs1.abs() > elim {
                if rs1 > zero {
                    return (y, -1);
                }
                if z.re < zero {
                    return (y, -1);
                }
                kdflg = 1;
                y[i] = czero;
                nz += 1;
                let str_cs = csi;
                csi = -csr;
                csr = str_cs;
                if i > 0 && (y[i - 1].re != zero || y[i - 1].im != zero) {
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

        // ── Compute S2 = PHI*(Ai*ASUM + CR2*(Ai'*BSUM)) (Fortran lines 6299-6318) ──
        // ARG is multiplied by CR2 before ZAIRY call
        let c2_arg = Complex::new(
            arg_arr[jj].re * cr2r - arg_arr[jj].im * cr2i,
            arg_arr[jj].re * cr2i + arg_arr[jj].im * cr2r,
        );

        let (ai, _nai) =
            zairy(c2_arg, AiryDerivative::Value, Scaling::Exponential).unwrap_or((czero, 0));
        let (dai, _ndai) =
            zairy(c2_arg, AiryDerivative::Derivative, Scaling::Exponential).unwrap_or((czero, 0));

        // STR = DAIR*BSUM - DAII*BSUMI (Fortran lines 6303-6304)
        let str_d = dai.re * bsum_arr[jj].re - dai.im * bsum_arr[jj].im;
        let sti_d = dai.re * bsum_arr[jj].im + dai.im * bsum_arr[jj].re;
        // PTR = STR*CR2R - STI*CR2I (multiply by CR2)
        let ptr = str_d * cr2r - sti_d * cr2i;
        let pti = str_d * cr2i + sti_d * cr2r;
        // STR = PTR + (AIR*ASUMR-AII*ASUMI)
        let str_a = ptr + (ai.re * asum_arr[jj].re - ai.im * asum_arr[jj].im);
        let sti_a = pti + (ai.re * asum_arr[jj].im + ai.im * asum_arr[jj].re);
        // PTR = STR*PHIR - STI*PHII
        let ptr2 = str_a * phi_arr[jj].re - sti_a * phi_arr[jj].im;
        let pti2 = str_a * phi_arr[jj].im + sti_a * phi_arr[jj].re;
        // S2 = PTR*CS
        let mut s2r = ptr2 * csr - pti2 * csi;
        let mut s2i = ptr2 * csi + pti2 * csr;

        // Scale by exp(S1) (Fortran lines 6313-6318)
        let str_exp = s1r.exp() * cssr[kflag - 1];
        let s1_re = str_exp * s1i.cos();
        let s1_im = str_exp * s1i.sin();
        let str_final = s2r * s1_re - s2i * s1_im;
        s2i = s1_re * s2i + s2r * s1_im;
        s2r = str_final;

        if kflag == 1 && zuchk(Complex::new(s2r, s2i), bry[0], tol) {
            // label 70: underflow
            if rs1 > zero {
                return (y, -1);
            }
            if z.re < zero {
                return (y, -1);
            }
            kdflg = 1;
            y[i] = czero;
            nz += 1;
            let str_cs = csi;
            csi = -csr;
            csr = str_cs;
            if i > 0 && (y[i - 1].re != zero || y[i - 1].im != zero) {
                y[i - 1] = czero;
                nz += 1;
            }
            continue;
        }

        // label 60: store result
        if yy <= zero {
            s2i = -s2i;
        }
        cy[kdflg - 1] = Complex::new(s2r, s2i);
        y[i] = Complex::new(s2r * csrr[kflag - 1], s2i * csrr[kflag - 1]);

        // CS rotation: (CSR,CSI) → (CSI,-CSR) (Fortran lines 6328-6330)
        let str_cs = csi;
        csi = -csr;
        csr = str_cs;

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

    let fn_at_exit = fnu + T::from((i_exit - 1) as f64).unwrap();
    let razr = one / zabs(Complex::new(zrr, zri));
    let str_rz = zrr * razr;
    let sti_rz = -zri * razr;
    let rzr = (str_rz + str_rz) * razr;
    let rzi = (sti_rz + sti_rz) * razr;
    let mut ckr = fn_at_exit * rzr;
    let mut cki = fn_at_exit * rzi;

    let ib = i_exit + 1;

    if ib <= n {
        // Test last member (Fortran lines 6368-6406)
        let fn_last = fnu + T::from((n - 1) as f64).unwrap();
        let ipard = if mr != 0 { 0 } else { 1 };
        let result_last = zunhj(Complex::new(znr, zni), fn_last, ipard, tol);

        let s1r_last = if kode == Scaling::Exponential {
            let str = zbr + result_last.zeta2.re;
            let sti = zbi + result_last.zeta2.im;
            let rast = fn_last / zabs(Complex::new(str, sti));
            let str2 = str * rast * rast;
            let _sti2 = -sti * rast * rast;
            result_last.zeta1.re - str2
        } else {
            result_last.zeta1.re - result_last.zeta2.re
        };

        let mut rs1 = s1r_last;
        if rs1.abs() > elim {
            if rs1 > zero {
                return (y, -1);
            }
            if z.re < zero {
                return (y, -1);
            }
            nz = n as i32;
            for item in y.iter_mut() {
                *item = czero;
            }
            return (y, nz);
        }
        if rs1.abs() >= alim {
            let aphi = zabs(result_last.phi);
            rs1 = rs1 + aphi.ln();
            if rs1.abs() >= elim {
                if rs1 > zero {
                    return (y, -1);
                }
                if z.re < zero {
                    return (y, -1);
                }
                nz = n as i32;
                for item in y.iter_mut() {
                    *item = czero;
                }
                return (y, nz);
            }
        }

        // Label 120: forward recurrence (Fortran lines 6407-6443)
        let mut s1 = cy[0];
        let mut s2 = cy[1];
        let mut c1r_rec = csrr[kflag - 1];
        let mut ascle = bry[kflag - 1];

        for y_item in y.iter_mut().take(n).skip(ib - 1) {
            let c2 = s2;
            s2 = Complex::new(
                ckr * c2.re - cki * c2.im + s1.re,
                ckr * c2.im + cki * c2.re + s1.im,
            );
            s1 = c2;
            ckr = ckr + rzr;
            cki = cki + rzi;
            let c2_scaled = Complex::new(s2.re * c1r_rec, s2.im * c1r_rec);
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
            s1 = Complex::new(s1.re * c1r_rec, s1.im * c1r_rec);
            s2 = c2_scaled;
            s1 = Complex::new(s1.re * cssr[kflag - 1], s1.im * cssr[kflag - 1]);
            s2 = Complex::new(s2.re * cssr[kflag - 1], s2.im * cssr[kflag - 1]);
            c1r_rec = csrr[kflag - 1];
        }
    }

    // ── Phase 3: Analytic continuation (Fortran lines 6444-6660) ──
    if mr == 0 {
        return (y, nz);
    }

    nz = 0;
    let fmr = T::from(mr as f64).unwrap();
    let sgn = -pi_t.copysign(fmr); // -DSIGN(PI, FMR)

    // CSGN and CSPN (Fortran lines 6455-6463)
    let mut csgni = sgn;
    if yy <= zero {
        csgni = -csgni;
    }
    let ifn = inu + n - 1;
    let ang2 = fnf * sgn;
    let mut cspnr = ang2.cos();
    let mut cspni = ang2.sin();
    if ifn % 2 != 0 {
        cspnr = -cspnr;
        cspni = -cspni;
    }

    // CS for I function (Fortran lines 6471-6478)
    let mut csr_ac = sar * csgni;
    let mut csi_ac = car * csgni;
    let in_idx = ifn % 4;
    let c2_ac_r = T::from(CIPR[in_idx]).unwrap();
    let c2_ac_i = T::from(CIPI[in_idx]).unwrap();
    let str_ac = csr_ac * c2_ac_r + csi_ac * c2_ac_i;
    csi_ac = -csr_ac * c2_ac_i + csi_ac * c2_ac_r;
    csr_ac = str_ac;

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
        let fn_val = fnu + T::from((kk_idx - 1) as f64).unwrap();

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
            let result_ac = zunhj(Complex::new(znr, zni), fn_val, 0, tol);
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
            let result_ac = zunhj(Complex::new(znr, zni), fn_val, 0, tol);
            phid = result_ac.phi;
            argd = result_ac.arg;
            zet1d = result_ac.zeta1;
            zet2d = result_ac.zeta2;
            asumd = result_ac.asum;
            bsumd = result_ac.bsum;
        }

        // S1 exponent for I function (Fortran lines 6514-6525)
        let (s1r, s1i) = if kode == Scaling::Exponential {
            let str = zbr + zet2d.re;
            let sti = zbi + zet2d.im;
            let rast = fn_val / zabs(Complex::new(str, sti));
            let str2 = str * rast * rast;
            let sti2 = -sti * rast * rast;
            (-zet1d.re + str2, -zet1d.im + sti2)
        } else {
            (-zet1d.re + zet2d.re, -zet1d.im + zet2d.im)
        };

        let mut rs1 = s1r;

        // ── Overflow/underflow test (Fortran lines 6530-6602) ──
        if rs1.abs() > elim {
            // label 280
            if rs1 > zero {
                return (y, -1);
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
                let s1s2_result = zs1s2(Complex::new(zrr, zri), s1_k, s2_k_final, asc, alim, iuf);
                s1_k = s1s2_result.s1;
                s2_k_final = s1s2_result.s2;
                nz += s1s2_result.nz;
                iuf = s1s2_result.iuf;
            }
            y[kk_idx - 1] = Complex::new(
                s1_k.re * cspnr - s1_k.im * cspni + s2_k_final.re,
                cspnr * s1_k.im + cspni * s1_k.re + s2_k_final.im,
            );
            kk_idx -= 1;
            cspnr = -cspnr;
            cspni = -cspni;
            let str_cs = csi_ac;
            csi_ac = -csr_ac;
            csr_ac = str_cs;

            if c2_save.re == zero && c2_save.im == zero {
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
            rs1 = rs1 + aphi.ln() - T::from(0.25).unwrap() * aarg.ln() - aic_t;
            if rs1.abs() > elim {
                if rs1 > zero {
                    return (y, -1);
                }
                // S2 = 0 path (goto 250)
                let c2_save = czero;
                cy[kdflg - 1] = czero;
                let s2_scaled = czero;

                if yy <= zero {
                    // no-op for zero
                }
                let s2_k = s2_scaled;
                let mut s1_k = y[kk_idx - 1];
                let mut s2_k_final = s2_k;
                if kode == Scaling::Exponential {
                    let s1s2_result =
                        zs1s2(Complex::new(zrr, zri), s1_k, s2_k_final, asc, alim, iuf);
                    s1_k = s1s2_result.s1;
                    s2_k_final = s1s2_result.s2;
                    nz += s1s2_result.nz;
                    iuf = s1s2_result.iuf;
                }
                y[kk_idx - 1] = Complex::new(
                    s1_k.re * cspnr - s1_k.im * cspni + s2_k_final.re,
                    cspnr * s1_k.im + cspni * s1_k.re + s2_k_final.im,
                );
                kk_idx -= 1;
                cspnr = -cspnr;
                cspni = -cspni;
                let str_cs = csi_ac;
                csi_ac = -csr_ac;
                csr_ac = str_cs;

                if c2_save.re == zero && c2_save.im == zero {
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
        let (ai, _nai) =
            zairy(argd, AiryDerivative::Value, Scaling::Exponential).unwrap_or((czero, 0));
        let (dai, _ndai) =
            zairy(argd, AiryDerivative::Derivative, Scaling::Exponential).unwrap_or((czero, 0));

        let str_d = dai.re * bsumd.re - dai.im * bsumd.im;
        let sti_d = dai.re * bsumd.im + dai.im * bsumd.re;
        let str_a = str_d + (ai.re * asumd.re - ai.im * asumd.im);
        let sti_a = sti_d + (ai.re * asumd.im + ai.im * asumd.re);
        let ptr = str_a * phid.re - sti_a * phid.im;
        let pti = str_a * phid.im + sti_a * phid.re;
        let mut s2r = ptr * csr_ac - pti * csi_ac;
        let mut s2i = ptr * csi_ac + pti * csr_ac;

        // Scale by exp(S1) (Fortran lines 6555-6560)
        let str_exp = s1r.exp() * cssr[iflag - 1];
        let s1_re = str_exp * s1i.cos();
        let s1_im = str_exp * s1i.sin();
        let str_final = s2r * s1_re - s2i * s1_im;
        s2i = s2r * s1_im + s2i * s1_re;
        s2r = str_final;

        if iflag == 1 && zuchk(Complex::new(s2r, s2i), bry[0], tol) {
            s2r = zero;
            s2i = zero;
        }

        // label 250 (Fortran lines 6566-6597)
        if yy <= zero {
            s2i = -s2i;
        }
        cy[kdflg - 1] = Complex::new(s2r, s2i);
        let c2_save = Complex::new(s2r, s2i);
        let s2_scaled = Complex::new(s2r * csrr[iflag - 1], s2i * csrr[iflag - 1]);

        // Add I and K functions (Fortran lines 6577-6590)
        let mut s1_k = y[kk_idx - 1];
        let mut s2_k = s2_scaled;
        if kode == Scaling::Exponential {
            let s1s2_result = zs1s2(Complex::new(zrr, zri), s1_k, s2_k, asc, alim, iuf);
            s1_k = s1s2_result.s1;
            s2_k = s1s2_result.s2;
            nz += s1s2_result.nz;
            iuf = s1s2_result.iuf;
        }
        y[kk_idx - 1] = Complex::new(
            s1_k.re * cspnr - s1_k.im * cspni + s2_k.re,
            cspnr * s1_k.im + cspni * s1_k.re + s2_k.im,
        );
        kk_idx -= 1;
        cspnr = -cspnr;
        cspni = -cspni;
        let str_cs = csi_ac;
        csi_ac = -csr_ac;
        csr_ac = str_cs;

        // KDFLG state machine (Fortran lines 6591-6597)
        if c2_save.re == zero && c2_save.im == zero {
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
        return (y, nz);
    }

    let mut s1 = cy[0];
    let mut s2 = cy[1];
    let mut csr_rec = csrr[iflag - 1];
    let mut ascle = bry[iflag - 1];
    let mut fn_rec = T::from((inu + il) as f64).unwrap();

    for _i in 0..il {
        let c2 = s2;
        let cfn = fn_rec + fnf;
        s2 = Complex::new(
            s1.re + cfn * (rzr * c2.re - rzi * c2.im),
            s1.im + cfn * (rzr * c2.im + rzi * c2.re),
        );
        s1 = c2;
        fn_rec = fn_rec - one;
        let ck_val = Complex::new(s2.re * csr_rec, s2.im * csr_rec);

        // Add K function (Fortran lines 6632-6639)
        let mut c1_k = y[kk_idx - 1];
        let mut c2_k = ck_val;
        if kode == Scaling::Exponential {
            let s1s2_result = zs1s2(Complex::new(zrr, zri), c1_k, c2_k, asc, alim, iuf);
            c1_k = s1s2_result.s1;
            c2_k = s1s2_result.s2;
            nz += s1s2_result.nz;
            iuf = s1s2_result.iuf;
        }
        y[kk_idx - 1] = Complex::new(
            c1_k.re * cspnr - c1_k.im * cspni + c2_k.re,
            c1_k.re * cspni + c1_k.im * cspnr + c2_k.im,
        );
        kk_idx -= 1;
        cspnr = -cspnr;
        cspni = -cspni;

        if iflag >= 3 {
            continue;
        }
        let c2m = ck_val.re.abs().max(ck_val.im.abs());
        if c2m <= ascle {
            continue;
        }
        iflag += 1;
        ascle = bry[iflag - 1];
        s1 = Complex::new(s1.re * csr_rec, s1.im * csr_rec);
        s2 = ck_val;
        s1 = Complex::new(s1.re * cssr[iflag - 1], s1.im * cssr[iflag - 1]);
        s2 = Complex::new(s2.re * cssr[iflag - 1], s2.im * cssr[iflag - 1]);
        csr_rec = csrr[iflag - 1];
    }

    (y, nz)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    const TOL: f64 = 2.220446049250313e-16;
    const ELIM: f64 = 700.9217936944459;
    const ALIM: f64 = 664.8716455337102;

    #[test]
    fn zunk2_basic_k_only() {
        // MR=0: K function only (region 2: |Im(z)| > |Re(z)|*sqrt(3))
        let z = Complex64::new(1.0, 10.0);
        let (y, nz) = zunk2(z, 90.0, Scaling::Unscaled, 0, 1, TOL, ELIM, ALIM);
        assert!(nz >= 0, "nz = {}", nz);
        assert!(y[0].re.is_finite());
    }

    #[test]
    fn zunk2_k_sequence() {
        let z = Complex64::new(1.0, 10.0);
        let (_y, nz) = zunk2(z, 90.0, Scaling::Unscaled, 0, 3, TOL, ELIM, ALIM);
        assert!(nz >= 0);
    }

    #[test]
    fn zunk2_analytic_continuation() {
        // MR=1: analytic continuation
        let z = Complex64::new(-1.0, 10.0);
        let (y, nz) = zunk2(z, 90.0, Scaling::Unscaled, 1, 1, TOL, ELIM, ALIM);
        assert!(nz >= 0, "nz = {}", nz);
        assert!(y[0].re.is_finite());
    }

    #[test]
    fn zunk2_scaled() {
        let z = Complex64::new(1.0, 10.0);
        let (y, nz) = zunk2(z, 90.0, Scaling::Exponential, 0, 1, TOL, ELIM, ALIM);
        assert!(nz >= 0);
        assert!(y[0].re.is_finite());
    }

    #[test]
    fn zunk2_real_argument() {
        let z = Complex64::new(10.0, 0.0);
        let (_y, nz) = zunk2(z, 90.0, Scaling::Unscaled, 0, 1, TOL, ELIM, ALIM);
        assert!(nz >= 0);
    }
}
