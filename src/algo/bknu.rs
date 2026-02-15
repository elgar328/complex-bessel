//! Core computation of the K Bessel function in the right half z-plane.
//!
//! Translation of Fortran ZBKNU from TOMS 644 / SLATEC.
//! This is the largest and most complex routine in the Amos algorithm.

// Exact Fortran constants — preserve verbatim.
#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]
#![allow(unused_assignments)]

use num_complex::Complex;

use crate::algo::gamln::gamln;
use crate::algo::kscl::zkscl;
use crate::algo::shch::zshch;
use crate::algo::uchk::zuchk;
use crate::machine::BesselFloat;
use crate::types::{BesselError, Scaling};
use crate::utils::{zabs, zdiv};

// ── Constants from Fortran DATA statements ──

const KMAX: i32 = 30;
const R1: f64 = 2.0;
const RTHPI: f64 = 1.25331413731550025; // sqrt(pi/2)
const SPI: f64 = 1.90985931710274403; // sqrt(6/pi)? Actually 6/(pi^2) related
const FPI: f64 = 1.89769999331517738; // 2^(1.75)/sqrt(pi)

use crate::algo::constants::{HPI, PI, R1M5, TTH};

/// Chebyshev coefficients for the f_0 series (small |DNU| case).
/// CC(1) = γ (Euler-Mascheroni constant).
#[rustfmt::skip]
const CC: [f64; 8] = [
    5.77215664901532861e-01,
   -4.20026350340952355e-02,
   -4.21977345555443367e-02,
    7.21894324666309954e-03,
   -2.15241674114950973e-04,
   -2.01348547807882387e-05,
    1.13302723198169588e-06,
    6.11609510448141582e-09,
];

/// Compute K Bessel function in the right half z-plane.
///
/// Computes K_{ν+j}(z) for j = 0, 1, ..., n-1, where Re(z) > 0.
///
/// # Parameters
/// - `z`: complex argument with Re(z) ≥ 0
/// - `fnu`: starting order ν ≥ 0
/// - `kode`: scaling mode (Unscaled = K_ν(z), Exponential = e^z · K_ν(z))
/// - `n`: number of sequence members
/// - `tol`: tolerance
/// - `elim`: underflow elimination threshold
/// - `alim`: overflow elimination threshold
///
/// # Returns
/// `(y, nz)` where `y` is the vector of K values and `nz` is the number
/// of underflowed (zeroed) leading components.
/// Returns `Err(BesselError::ConvergenceFailure)` if the forward recurrence
/// loop in the Miller algorithm does not converge.
pub(crate) fn zbknu<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    n: usize,
    tol: T,
    elim: T,
    alim: T,
) -> Result<(Vec<Complex<T>>, usize), BesselError> {
    // Convenience conversions
    let zero = T::zero();
    let one = T::one();
    let two = T::from(2.0).unwrap();
    let half = T::from(0.5).unwrap();
    let czero = Complex::new(zero, zero);

    let caz = zabs(z);
    let csclr = one / tol;
    let crscr = tol;

    // Scaling/unscaling arrays (1-indexed in Fortran, 0-indexed here)
    let cssr = [csclr, one, crscr];
    let csrr = [crscr, one, csclr];

    // D1MACH(1) = MACH_TINY, D1MACH(2) = MACH_HUGE
    let bry = [
        T::from(1.0e3).unwrap() * T::MACH_TINY / tol,
        one / (T::from(1.0e3).unwrap() * T::MACH_TINY / tol),
        T::MACH_HUGE,
    ];

    let mut y = vec![czero; n];
    let mut nz: usize = 0;
    let mut iflag = 0_i32;
    let mut koded = kode;
    let rcaz = one / caz;
    let str_val = z.re * rcaz;
    let sti = -z.im * rcaz;
    let rzr = (str_val + str_val) * rcaz;
    let rzi = (sti + sti) * rcaz;
    let rz = Complex::new(rzr, rzi);

    let inu = (fnu + half).floor().to_i32().unwrap();
    let dnu = fnu - T::from(inu).unwrap();
    let dnu2 = if dnu.abs() > tol { dnu * dnu } else { zero };

    // ── Branch: series (|z| ≤ R1) vs Miller/asymptotic (|z| > R1) ──
    if dnu.abs() != half && caz <= T::from(R1).unwrap() {
        // ══════════════════════════════════════════════════════════════
        // SERIES FOR |Z| ≤ R1  (Fortran lines 82-216)
        // ══════════════════════════════════════════════════════════════
        let mut fc = one;
        let rz_log = Complex::new(rzr, rzi).ln();
        let fmu = Complex::new(rz_log.re * dnu, rz_log.im * dnu);
        let (csh, cch) = zshch(fmu);

        let smu;
        if dnu != zero {
            fc = dnu * T::from(PI).unwrap();
            fc = fc / fc.sin();
            smu = Complex::new(csh.re / dnu, csh.im / dnu);
        } else {
            // DNU=0: sinh(x·DNU)/DNU → x as DNU→0 by L'Hôpital.
            // Fortran reuses SMUR/SMUI from ZLOG(RZ) here (variable reuse).
            smu = rz_log;
        }

        let a2 = one + dnu;
        // T2 = 1/Γ(1+DNU), T1 = 1/Γ(1-DNU)
        let t2 = (-gamln(a2)?).exp();
        let t1 = one / (t2 * fc);

        // Compute G1 and G2
        let g1;
        if dnu.abs() > T::from(0.1).unwrap() {
            // Large |DNU|: direct formula
            g1 = (t1 - t2) / (dnu + dnu);
        } else {
            // Small |DNU|: Chebyshev series for f_0
            let mut ak = one;
            let mut s = T::from(CC[0]).unwrap();
            for cc_k in &CC[1..] {
                ak = ak * dnu2;
                let tm = T::from(*cc_k).unwrap() * ak;
                s = s + tm;
                if tm.abs() < tol {
                    break;
                }
            }
            g1 = -s;
        }
        let g2 = (t1 + t2) * half;

        // F = FC * (CCH*G1 + SMU*G2)
        let mut fr = fc * (cch.re * g1 + smu.re * g2);
        let mut fi = fc * (cch.im * g1 + smu.im * g2);

        // P = 0.5 * exp(FMU) / T2
        let efmu = fmu.exp();
        let mut pr = half * efmu.re / t2;
        let mut pi_val = half * efmu.im / t2;

        // Q = 0.5 / (exp(FMU) * T1)
        let q_tmp = zdiv(Complex::new(half, zero), efmu);
        let mut qr = q_tmp.re / t1;
        let mut qi = q_tmp.im / t1;

        let mut s1 = Complex::new(fr, fi);
        let mut s2 = Complex::new(pr, pi_val);
        let mut ak = one;
        let mut a1 = one;
        let mut ck = Complex::new(one, zero);
        let mut bk = one - dnu2;

        let mut kflag: usize; // 0-based index into cssr/csrr arrays

        if inu > 0 || n > 1 {
            // ── Generate K(DNU,Z) and K(DNU+1,Z) for forward recurrence ──
            // (Fortran label 80-100)
            if caz >= tol {
                let cz = z * z * T::from(0.25).unwrap();
                let t1_sq = T::from(0.25).unwrap() * caz * caz;
                loop {
                    fr = (fr * ak + pr + qr) / bk;
                    fi = (fi * ak + pi_val + qi) / bk;
                    let str_ak = one / (ak - dnu);
                    pr = pr * str_ak;
                    pi_val = pi_val * str_ak;
                    let str_ak2 = one / (ak + dnu);
                    qr = qr * str_ak2;
                    qi = qi * str_ak2;
                    let rak = one / ak;
                    let ck_new = Complex::new(
                        (ck.re * cz.re - ck.im * cz.im) * rak,
                        (ck.re * cz.im + ck.im * cz.re) * rak,
                    );
                    ck = ck_new;
                    s1 = s1 + Complex::new(ck.re * fr - ck.im * fi, ck.re * fi + ck.im * fr);
                    let str_s2 = pr - fr * ak;
                    let sti_s2 = pi_val - fi * ak;
                    s2 = s2
                        + Complex::new(
                            ck.re * str_s2 - ck.im * sti_s2,
                            ck.re * sti_s2 + ck.im * str_s2,
                        );
                    a1 = a1 * t1_sq * rak;
                    bk = bk + ak + ak + one;
                    ak = ak + one;
                    if a1 <= tol {
                        break;
                    }
                }
            }

            // Scale and prepare for forward recurrence
            kflag = 1; // KFLAG=2 in Fortran (0-based: 1)
            let a1_check = (fnu + one) * smu.re.abs();
            if a1_check > alim {
                kflag = 2; // KFLAG=3
            }
            let scale = cssr[kflag];
            let p2 = s2 * scale;
            s2 = Complex::new(p2.re * rzr - p2.im * rzi, p2.re * rzi + p2.im * rzr);
            s1 = s1 * scale;

            if matches!(koded, Scaling::Exponential) {
                let ez = z.exp();
                s1 = s1 * ez;
                s2 = s2 * ez;
            }
        } else {
            // ── Generate K(FNU,Z), 0 ≤ FNU < 0.5 and N=1 ──
            // (Fortran label 60-70)
            if caz >= tol {
                let cz = z * z * T::from(0.25).unwrap();
                let t1_sq = T::from(0.25).unwrap() * caz * caz;
                loop {
                    fr = (fr * ak + pr + qr) / bk;
                    fi = (fi * ak + pi_val + qi) / bk;
                    let str_ak = one / (ak - dnu);
                    pr = pr * str_ak;
                    pi_val = pi_val * str_ak;
                    let str_ak2 = one / (ak + dnu);
                    qr = qr * str_ak2;
                    qi = qi * str_ak2;
                    let rak = one / ak;
                    let ck_new = Complex::new(
                        (ck.re * cz.re - ck.im * cz.im) * rak,
                        (ck.re * cz.im + ck.im * cz.re) * rak,
                    );
                    ck = ck_new;
                    s1 = s1 + Complex::new(ck.re * fr - ck.im * fi, ck.re * fi + ck.im * fr);
                    a1 = a1 * t1_sq * rak;
                    bk = bk + ak + ak + one;
                    ak = ak + one;
                    if a1 <= tol {
                        break;
                    }
                }
            }

            y[0] = s1;
            if matches!(koded, Scaling::Exponential) {
                let ez = z.exp();
                y[0] = s1 * ez;
            }
            return Ok((y, nz));
        }

        // ── Forward recurrence: label 210 ──
        return forward_recurrence(
            z, fnu, &koded, n, &mut y, &mut nz, tol, elim, alim, iflag, s1, s2, rz, inu, kflag,
            &cssr, &csrr, &bry,
        );
    }

    // ══════════════════════════════════════════════════════════════════
    // MILLER ALGORITHM / ASYMPTOTIC PATH (Fortran label 110+)
    // |z| > R1 or |DNU| == 0.5
    // ══════════════════════════════════════════════════════════════════

    // COEF = RTHPI / sqrt(z) = sqrt(π/2) / sqrt(z) = sqrt(π/(2z))
    let sz = z.sqrt();
    let coef_base = zdiv(Complex::new(T::from(RTHPI).unwrap(), zero), sz);

    let mut kflag: usize = 1; // KFLAG=2 in Fortran

    let coef;
    if matches!(koded, Scaling::Unscaled) {
        if z.re > alim {
            // Label 290: scale by exp(z), iflag=1
            koded = Scaling::Exponential;
            iflag = 1;
            kflag = 1; // stays KFLAG=2
            coef = coef_base; // will apply exp later
        } else {
            // Multiply coef by exp(-z) * cssr[kflag]
            let exp_neg_z_scaled = {
                let e = (-z.re).exp() * cssr[kflag];
                Complex::new(e * z.im.cos(), -e * z.im.sin())
            };
            coef = Complex::new(
                coef_base.re * exp_neg_z_scaled.re - coef_base.im * exp_neg_z_scaled.im,
                coef_base.re * exp_neg_z_scaled.im + coef_base.im * exp_neg_z_scaled.re,
            );
        }
    } else {
        coef = coef_base;
    }

    // Check for DNU = ±0.5 special case (label 300)
    if dnu.abs() == half {
        let s1 = coef;
        let s2 = coef;
        return forward_recurrence(
            z, fnu, &koded, n, &mut y, &mut nz, tol, elim, alim, iflag, s1, s2, rz, inu, kflag,
            &cssr, &csrr, &bry,
        );
    }

    // Check if cos(π·DNU) == 0 or FHS == 0
    let ak_cos = T::from(PI).unwrap() * dnu;
    let ak_val = ak_cos.cos().abs();
    let fhs = (T::from(0.25).unwrap() - dnu2).abs();
    if ak_val == zero || fhs == zero {
        // DNU = ±0.5 effective (label 300)
        let s1 = coef;
        let s2 = coef;
        return forward_recurrence(
            z, fnu, &koded, n, &mut y, &mut nz, tol, elim, alim, iflag, s1, s2, rz, inu, kflag,
            &cssr, &csrr, &bry,
        );
    }

    // ── Compute R2 = F(E) for determining backward index K ──
    // T1 = (I1MACH(14) - 1) * D1MACH(5) * 3.321928094
    // D1MACH(5) = log10(2) for binary
    let r1m5 = T::from(R1M5).unwrap();
    let t1_raw = T::from(T::MACH_DIGITS - 1).unwrap() * r1m5 * T::from(3.321928094).unwrap();
    let t1_clamped = t1_raw
        .max(T::from(12.0).unwrap())
        .min(T::from(60.0).unwrap());
    let t2_miller = T::from(TTH).unwrap() * t1_clamped - T::from(6.0).unwrap();

    let t1_angle = if z.re != zero {
        (z.im / z.re).abs().atan()
    } else {
        T::from(HPI).unwrap()
    };

    let fk: T;
    let mut fhs_miller = fhs;

    if t2_miller <= caz {
        // ── Forward recurrence to find backward index K (label 150) ──
        let etest = ak_val / (T::from(PI).unwrap() * caz * tol);
        let mut fk_val = one;
        if etest >= one {
            let mut fks = two;
            let mut ckr_fwd = caz + caz + two;
            let mut p1r = zero;
            let mut p2r = one;
            let mut converged = false;

            for _i in 0..KMAX {
                let ak_inner = fhs_miller / fks;
                let cbr = ckr_fwd / (fk_val + one);
                let ptr = p2r;
                p2r = cbr * p2r - p1r * ak_inner;
                p1r = ptr;
                ckr_fwd = ckr_fwd + two;
                fks = fks + fk_val + fk_val + two;
                fhs_miller = fhs_miller + fk_val + fk_val;
                fk_val = fk_val + one;
                let str_test = p2r.abs() * fk_val;
                if etest < str_test {
                    converged = true;
                    break;
                }
            }
            if !converged {
                // Label 310: convergence failure
                return Err(BesselError::ConvergenceFailure);
            }
            // Label 160
            fk_val = fk_val + T::from(SPI).unwrap() * t1_angle * (t2_miller / caz).sqrt();
            fhs_miller = (T::from(0.25).unwrap() - dnu2).abs();
        }
        fk = fk_val;
    } else {
        // ── Compute backward index K for |z| < R2 (label 170) ──
        let a2 = caz.sqrt();
        let ak_inner = T::from(FPI).unwrap() * ak_val / (tol * a2.sqrt());
        let three = T::from(3.0).unwrap();
        let aa = three * t1_angle / (one + caz);
        let bb = T::from(14.7).unwrap() * t1_angle / (T::from(28.0).unwrap() + caz);
        let ak_log =
            (ak_inner.ln() + caz * aa.cos() / (one + T::from(0.008).unwrap() * caz)) / bb.cos();
        fk = T::from(0.12125).unwrap() * ak_log * ak_log / caz + T::from(1.5).unwrap();
    }

    // ── Backward recurrence loop (Miller algorithm, label 180-190) ──
    let k = fk.floor().to_i32().unwrap();
    let fk_int = T::from(k).unwrap();
    let mut fks = fk_int * fk_int;
    let mut p1 = czero;
    let mut p2 = Complex::new(tol, zero);
    let mut cs = p2;
    let mut fk_cur = fk_int;

    for _i in 0..k {
        let a1 = fks - fk_cur;
        let ak_inner = (fks + fk_cur) / (a1 + fhs_miller);
        let rak = two / (fk_cur + one);
        let cb = Complex::new((fk_cur + z.re) * rak, z.im * rak);
        let pt = p2;
        p2 = Complex::new(
            (pt.re * cb.re - pt.im * cb.im - p1.re) * ak_inner,
            (pt.im * cb.re + pt.re * cb.im - p1.im) * ak_inner,
        );
        p1 = pt;
        cs = cs + p2;
        fks = a1 - fk_cur + one;
        fk_cur = fk_cur - one;
    }

    // ── Compute (P2/CS) with better scaling (label after 190) ──
    let tm = zabs(cs);
    let ptr = one / tm;
    let s1_raw = p2 * ptr;
    let cs_conj = Complex::new(cs.re * ptr, -cs.im * ptr);
    // S1 = COEF * (P2/|CS|) * (conj(CS)/|CS|)
    let tmp = Complex::new(
        coef.re * s1_raw.re - coef.im * s1_raw.im,
        coef.re * s1_raw.im + coef.im * s1_raw.re,
    );
    let s1 = Complex::new(
        tmp.re * cs_conj.re - tmp.im * cs_conj.im,
        tmp.re * cs_conj.im + tmp.im * cs_conj.re,
    );

    if inu == 0 && n == 1 {
        // Direct output path
        let zd = z;
        if iflag == 1 {
            return handle_iflag1_final(
                zd, fnu, n, &mut y, &mut nz, tol, elim, alim, s1, s1, rz, &cssr, &csrr, &bry,
            );
        }
        // Label 240: store and return
        let str_scale = csrr[kflag];
        y[0] = s1 * str_scale;
        return Ok((y, nz));
    }

    // ── Compute P1/P2 ratio for S2 (label 200) ──
    let tm2 = zabs(p2);
    let ptr2 = one / tm2;
    let p1_scaled = p1 * ptr2;
    let p2_conj = Complex::new(p2.re * ptr2, -p2.im * ptr2);
    let pt_ratio = Complex::new(
        p1_scaled.re * p2_conj.re - p1_scaled.im * p2_conj.im,
        p1_scaled.re * p2_conj.im + p1_scaled.im * p2_conj.re,
    );

    // S2 = S1 * (DNU + 0.5 - P1/P2·conj(P2)/|P2|) / z + S1... wait
    // Fortran: STR = DNU + 0.5 - PTR; STI = -PTI
    //          CALL ZDIV(STR, STI, ZR, ZI, STR, STI)
    //          STR = STR + 1.0
    //          CALL ZMLT(STR, STI, S1R, S1I, S2R, S2I)
    let ratio = Complex::new(dnu + half - pt_ratio.re, -pt_ratio.im);
    let ratio_div_z = zdiv(ratio, z) + Complex::new(one, zero);
    let s2 = Complex::new(
        ratio_div_z.re * s1.re - ratio_div_z.im * s1.im,
        ratio_div_z.re * s1.im + ratio_div_z.im * s1.re,
    );

    // ── Forward recurrence (label 210) ──
    forward_recurrence(
        z, fnu, &koded, n, &mut y, &mut nz, tol, elim, alim, iflag, s1, s2, rz, inu, kflag, &cssr,
        &csrr, &bry,
    )
}

/// Forward recurrence on the three-term recursion relation.
///
/// Computes K_{ν+j}(z) for j = 0, 1, ..., n-1 starting from S1 = K(DNU, z)
/// and S2 = K(DNU+1, z), using the recurrence:
///   K_{ν+1}(z) = (2(ν+1)/z) K_ν(z) + K_{ν-1}(z)
///
/// Handles scaling (KFLAG) and underflow (IFLAG=1) paths.
#[allow(clippy::too_many_arguments)]
fn forward_recurrence<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    _koded: &Scaling,
    n: usize,
    y: &mut [Complex<T>],
    nz: &mut usize,
    tol: T,
    elim: T,
    alim: T,
    iflag: i32,
    mut s1: Complex<T>,
    mut s2: Complex<T>,
    rz: Complex<T>,
    inu: i32,
    mut kflag: usize,
    cssr: &[T; 3],
    csrr: &[T; 3],
    bry: &[T; 3],
) -> Result<(Vec<Complex<T>>, usize), BesselError> {
    let zero = T::zero();
    let one = T::one();
    let _czero = Complex::new(zero, zero);

    // CK = (DNU + 1) * RZ
    let str_init = fnu - T::from(inu).unwrap() + one;
    let mut ck = Complex::new(str_init * rz.re, str_init * rz.im);

    let mut inu_adj = inu;
    if n == 1 {
        inu_adj -= 1;
    }

    if inu_adj > 0 {
        if iflag == 1 {
            // ── IFLAG=1: forward recurrence on scaled values (label 261) ──
            return iflag1_recurrence(
                z, fnu, n, y, nz, tol, elim, alim, s1, s2, rz, ck, inu_adj, cssr, csrr, bry,
            );
        }

        // ── Normal forward recurrence (labels 225-230) ──
        let inub = 0_i32; // 0-based start
        let mut p1r = csrr[kflag];
        let mut ascle = bry[kflag];

        for _i in inub..inu_adj {
            let st = s2;
            s2 = Complex::new(
                ck.re * st.re - ck.im * st.im + s1.re,
                ck.im * st.re + ck.re * st.im + s1.im,
            );
            s1 = st;
            ck = ck + rz;

            if kflag < 2 {
                let p2_test = s2 * p1r;
                let p2m = p2_test.re.abs().max(p2_test.im.abs());
                if p2m > ascle {
                    kflag += 1;
                    ascle = bry[kflag];
                    s1 = s1 * p1r;
                    s2 = p2_test;
                    let str_scale = cssr[kflag];
                    s1 = s1 * str_scale;
                    s2 = s2 * str_scale;
                    p1r = csrr[kflag];
                }
            }
        }

        if n == 1 {
            s1 = s2;
        }
    } else {
        // inu_adj <= 0
        if n == 1 {
            s1 = s2;
        }
        // Fortran label 215: IF(IFLAG.EQ.1) GO TO 270 (ZKSCL path)
        if iflag == 1 {
            return handle_iflag1_final(
                z, fnu, n, y, nz, tol, elim, alim, s1, s2, rz, cssr, csrr, bry,
            );
        }
    }

    // ── Label 240: store results and continue recurrence ──
    let str_scale = csrr[kflag];
    y[0] = s1 * str_scale;
    if n == 1 {
        return Ok((y.to_vec(), *nz));
    }
    y[1] = s2 * str_scale;
    if n == 2 {
        return Ok((y.to_vec(), *nz));
    }

    // ── Label 250-260: additional values via recurrence ──
    let mut kk = 2_usize; // next index to fill (0-based)
    loop {
        if kk >= n {
            return Ok((y.to_vec(), *nz));
        }

        let mut p1r = csrr[kflag];
        let mut ascle = bry[kflag];

        for i in kk..n {
            let p2_val = s2;
            s2 = Complex::new(
                ck.re * p2_val.re - ck.im * p2_val.im + s1.re,
                ck.im * p2_val.re + ck.re * p2_val.im + s1.im,
            );
            s1 = p2_val;
            ck = ck + rz;
            let p2_scaled = s2 * p1r;
            y[i] = p2_scaled;

            if kflag < 2 {
                let p2m = p2_scaled.re.abs().max(p2_scaled.im.abs());
                if p2m > ascle {
                    kflag += 1;
                    ascle = bry[kflag];
                    s1 = s1 * p1r;
                    s2 = p2_scaled;
                    let str_s = cssr[kflag];
                    s1 = s1 * str_s;
                    s2 = s2 * str_s;
                    p1r = csrr[kflag];
                    kk = i + 1;
                    // Fortran: GO TO 250 (restart outer loop)
                    break;
                }
            }

            if i == n - 1 {
                return Ok((y.to_vec(), *nz));
            }
        }
        // If we didn't break early, we've filled all values
        if kk >= n {
            return Ok((y.to_vec(), *nz));
        }
    }
}

/// Handle IFLAG=1 forward recurrence on scaled values (Fortran labels 261-280).
#[allow(clippy::too_many_arguments)]
fn iflag1_recurrence<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    n: usize,
    y: &mut [Complex<T>],
    nz: &mut usize,
    tol: T,
    elim: T,
    alim: T,
    mut s1: Complex<T>,
    mut s2: Complex<T>,
    rz: Complex<T>,
    mut ck: Complex<T>,
    inu: i32,
    cssr: &[T; 3],
    csrr: &[T; 3],
    bry: &[T; 3],
) -> Result<(Vec<Complex<T>>, usize), BesselError> {
    let zero = T::zero();
    let _one = T::one();
    let czero = Complex::new(zero, zero);
    let half = T::from(0.5).unwrap();

    let helim = half * elim;
    let elm = (-elim).exp();
    let celmr = elm;
    let ascle = bry[0]; // BRY(1)
    let mut zd = z;
    // Fortran IC=-1 with 1-based loop; 0-based equivalent needs IC=-2
    // so that `ic == i - 1` is false when i=0 (first success).
    let mut ic: i32 = -2;
    let mut j: usize = 1; // toggles between 0 and 1 (Fortran J toggles 1,2)
    let mut cy = [czero; 2];

    let mut found_two = false;
    let mut inub_resume = 0_i32;

    for i in 0..inu {
        let st = s2;
        s2 = Complex::new(
            st.re * ck.re - st.im * ck.im + s1.re,
            st.im * ck.re + st.re * ck.im + s1.im,
        );
        s1 = st;
        ck = ck + rz;

        let az = zabs(s2);
        let alas = az.ln();
        let p2r = -zd.re + alas;

        if p2r >= -elim {
            let s2_log = s2.ln();
            let p2 = Complex::new(-zd.re + s2_log.re, -zd.im + s2_log.im);
            let p2m = p2.re.exp() / tol;
            let p1_test = Complex::new(p2m * p2.im.cos(), p2m * p2.im.sin());
            if !zuchk(p1_test, ascle, tol) {
                j = 1 - j; // toggle: Fortran J = 3 - J
                cy[j] = p1_test;
                if ic == i - 1 {
                    // Two consecutive on-scale values found (label 264)
                    found_two = true;
                    inub_resume = i + 1;
                    break;
                }
                ic = i;
                continue;
            }
        }

        // Value underflowed; rescale if magnitude is getting too large
        if alas >= helim {
            zd = Complex::new(zd.re - elim, zd.im);
            s1 = s1 * celmr;
            s2 = s2 * celmr;
        }
    }

    if found_two {
        // Label 264: resume normal recurrence
        let kflag: usize = 0; // KFLAG=1 in Fortran
        s2 = cy[j];
        j = 1 - j;
        s1 = cy[j];
        let inub = inub_resume;

        if inub < inu {
            // Resume normal forward recurrence at label 225
            let mut p1r = csrr[kflag];
            let mut ascle_inner = bry[kflag];
            let mut kflag_inner = kflag;

            for _i in inub..inu {
                let st = s2;
                s2 = Complex::new(
                    ck.re * st.re - ck.im * st.im + s1.re,
                    ck.im * st.re + ck.re * st.im + s1.im,
                );
                s1 = st;
                ck = ck + rz;

                if kflag_inner < 2 {
                    let p2_test = s2 * p1r;
                    let p2m = p2_test.re.abs().max(p2_test.im.abs());
                    if p2m > ascle_inner {
                        kflag_inner += 1;
                        ascle_inner = bry[kflag_inner];
                        s1 = s1 * p1r;
                        s2 = p2_test;
                        let str_s = cssr[kflag_inner];
                        s1 = s1 * str_s;
                        s2 = s2 * str_s;
                        p1r = csrr[kflag_inner];
                    }
                }
            }

            if n == 1 {
                s1 = s2;
            }

            // Label 240
            let str_scale = csrr[kflag_inner];
            y[0] = s1 * str_scale;
            if n == 1 {
                return Ok((y.to_vec(), *nz));
            }
            y[1] = s2 * str_scale;
            if n == 2 {
                return Ok((y.to_vec(), *nz));
            }

            // Continue with label 250-260 recurrence (with KFLAG tracking)
            let t2_val = fnu + T::from(2).unwrap(); // FNU + (KK-1) where KK=3
            ck = Complex::new(t2_val * rz.re, t2_val * rz.im);
            let mut kk_idx = 2_usize;
            loop {
                if kk_idx >= n {
                    return Ok((y.to_vec(), *nz));
                }
                let mut p1r = csrr[kflag_inner];
                let mut ascle_loop = bry[kflag_inner];

                for i in kk_idx..n {
                    let p2_val = s2;
                    s2 = Complex::new(
                        ck.re * p2_val.re - ck.im * p2_val.im + s1.re,
                        ck.im * p2_val.re + ck.re * p2_val.im + s1.im,
                    );
                    s1 = p2_val;
                    ck = ck + rz;
                    let p2_scaled = s2 * p1r;
                    y[i] = p2_scaled;

                    if kflag_inner < 2 {
                        let p2m = p2_scaled.re.abs().max(p2_scaled.im.abs());
                        if p2m > ascle_loop {
                            kflag_inner += 1;
                            ascle_loop = bry[kflag_inner];
                            s1 = s1 * p1r;
                            s2 = p2_scaled;
                            let str_s = cssr[kflag_inner];
                            s1 = s1 * str_s;
                            s2 = s2 * str_s;
                            p1r = csrr[kflag_inner];
                            kk_idx = i + 1;
                            break;
                        }
                    }

                    if i == n - 1 {
                        return Ok((y.to_vec(), *nz));
                    }
                }
                if kk_idx >= n {
                    return Ok((y.to_vec(), *nz));
                }
            }
        }

        if n == 1 {
            s1 = s2;
        }

        // Label 240
        let str_scale = csrr[kflag];
        y[0] = s1 * str_scale;
        if n == 1 {
            return Ok((y.to_vec(), *nz));
        }
        y[1] = s2 * str_scale;
        if n == 2 {
            return Ok((y.to_vec(), *nz));
        }
        return Ok((y.to_vec(), *nz));
    }

    // Loop completed without finding two consecutive (label 270)
    if n == 1 {
        s1 = s2;
    }

    handle_iflag1_final(
        zd, fnu, n, y, nz, tol, elim, alim, s1, s2, rz, cssr, csrr, bry,
    )
}

/// Handle IFLAG=1 final path: store values and call ZKSCL (labels 270-280).
#[allow(clippy::too_many_arguments)]
fn handle_iflag1_final<T: BesselFloat>(
    zd: Complex<T>,
    fnu: T,
    n: usize,
    y: &mut [Complex<T>],
    nz: &mut usize,
    tol: T,
    elim: T,
    _alim: T,
    s1: Complex<T>,
    s2: Complex<T>,
    rz: Complex<T>,
    cssr: &[T; 3],
    csrr: &[T; 3],
    bry: &[T; 3],
) -> Result<(Vec<Complex<T>>, usize), BesselError> {
    let _zero = T::zero();
    let _one = T::one();

    y[0] = s1;
    if n > 1 {
        y[1] = s2;
    }

    // Label 280: call ZKSCL
    let ascle = bry[0];
    *nz = zkscl(zd, fnu, y, rz, ascle, tol, elim);
    let inu_remaining = n as i32 - *nz as i32;
    if inu_remaining <= 0 {
        return Ok((y.to_vec(), *nz));
    }

    let kk = *nz; // 0-based index of first non-zero value
    let mut s1_local = y[kk];
    y[kk] = s1_local * csrr[0];
    if inu_remaining == 1 {
        return Ok((y.to_vec(), *nz));
    }

    let kk2 = *nz + 1;
    let mut s2_local = y[kk2];
    y[kk2] = s2_local * csrr[0];
    if inu_remaining == 2 {
        return Ok((y.to_vec(), *nz));
    }

    // Continue recurrence for remaining values (label 250-260 with KFLAG tracking)
    let t2_val = fnu + T::from(kk2).unwrap();
    let mut ck = Complex::new(t2_val * rz.re, t2_val * rz.im);
    let mut kflag: usize = 0; // KFLAG=1 in Fortran

    let mut kk_idx = kk2 + 1; // 0-based next index to fill
    loop {
        if kk_idx >= n {
            return Ok((y.to_vec(), *nz));
        }
        let mut p1r = csrr[kflag];
        let mut ascle = bry[kflag];

        for i in kk_idx..n {
            let p2_val = s2_local;
            s2_local = Complex::new(
                ck.re * p2_val.re - ck.im * p2_val.im + s1_local.re,
                ck.im * p2_val.re + ck.re * p2_val.im + s1_local.im,
            );
            s1_local = p2_val;
            ck = ck + rz;
            let p2_scaled = s2_local * p1r;
            y[i] = p2_scaled;

            if kflag < 2 {
                let p2m = p2_scaled.re.abs().max(p2_scaled.im.abs());
                if p2m > ascle {
                    kflag += 1;
                    ascle = bry[kflag];
                    s1_local = s1_local * p1r;
                    s2_local = p2_scaled;
                    let str_s = cssr[kflag];
                    s1_local = s1_local * str_s;
                    s2_local = s2_local * str_s;
                    p1r = csrr[kflag];
                    kk_idx = i + 1;
                    break;
                }
            }

            if i == n - 1 {
                return Ok((y.to_vec(), *nz));
            }
        }
        if kk_idx >= n {
            return Ok((y.to_vec(), *nz));
        }
    }
}
