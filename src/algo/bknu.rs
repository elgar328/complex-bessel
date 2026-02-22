//! Core computation of the K Bessel function in the right half z-plane.
//!
//! Translation of Fortran ZBKNU from TOMS 644 (zbsubs.f lines 2391-2959).
//! This is the largest and most complex routine in the Amos algorithm.

// Exact Fortran constants — preserve verbatim.
#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]

use num_complex::Complex;

use crate::algo::gamln::gamln;
use crate::algo::kscl::zkscl;
use crate::algo::shch::zshch;
use crate::algo::uchk::zuchk;
use crate::machine::BesselFloat;
use crate::types::{Error, Scaling};
use crate::utils::{reciprocal_z, zabs, zdiv};

// ── Constants from Fortran DATA statements ──

const KMAX: i32 = 30;
const R1: f64 = 2.0;
const RTHPI: f64 = 1.25331413731550025; // sqrt(pi/2)
const SPI: f64 = 1.90985931710274403; // 6/π (Fortran ZBKNU DATA SPI)
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
/// Returns `Err(Error::ConvergenceFailure)` if the forward recurrence
/// loop in the Miller algorithm does not converge.
pub(crate) fn zbknu<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    y: &mut [Complex<T>],
    tol: T,
    elim: T,
    alim: T,
) -> Result<usize, Error> {
    // Convenience conversions
    let zero = T::zero();
    let one = T::one();
    let two = T::from_f64(2.0);
    let half = T::from_f64(0.5);
    let czero = Complex::new(zero, zero);

    let n = y.len();
    y.fill(czero);

    let caz = zabs(z);
    let csclr = one / tol;
    let crscr = tol;

    // Scaling/unscaling arrays (1-indexed in Fortran, 0-indexed here)
    let cssr = [csclr, one, crscr];
    let csrr = [crscr, one, csclr];

    // D1MACH(1) = MACH_TINY, D1MACH(2) = MACH_HUGE
    let bry = [
        T::from_f64(1.0e3) * T::MACH_TINY / tol,
        one / (T::from_f64(1.0e3) * T::MACH_TINY / tol),
        T::MACH_HUGE,
    ];
    let mut nz: usize = 0;
    let mut iflag = 0_i32;
    let mut koded = kode;
    let rz = reciprocal_z(z);

    // Safety: fnu is finite and < ~1e15 per upper-interface checks
    let inu = (fnu + half).floor().to_i32().unwrap();
    let dnu = fnu - T::from_f64(inu as f64);
    let dnu2 = if dnu.abs() > tol { dnu * dnu } else { zero };

    // ── Branch: series (|z| ≤ R1) vs Miller/asymptotic (|z| > R1) ──
    if dnu.abs() != half && caz <= T::from_f64(R1) {
        // ══════════════════════════════════════════════════════════════
        // SERIES FOR |Z| ≤ R1  (Fortran lines 82-216)
        // ══════════════════════════════════════════════════════════════
        let mut fc = one;
        let rz_log = rz.ln();
        let fmu = rz_log * dnu;
        let (csh, cch) = zshch(fmu);

        let smu;
        if dnu != zero {
            fc = dnu * T::from_f64(PI);
            fc = fc / fc.sin();
            smu = csh / dnu;
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
        if dnu.abs() > T::from_f64(0.1) {
            // Large |DNU|: direct formula
            g1 = (t1 - t2) / (dnu + dnu);
        } else {
            // Small |DNU|: Chebyshev series for f_0
            let mut ak = one;
            let mut s = T::from_f64(CC[0]);
            for cc_k in &CC[1..] {
                ak = ak * dnu2;
                let tm = T::from_f64(*cc_k) * ak;
                s = s + tm;
                if tm.abs() < tol {
                    break;
                }
            }
            g1 = -s;
        }
        let g2 = (t1 + t2) * half;

        // F = FC * (CCH*G1 + SMU*G2)
        let efmu = fmu.exp();
        let mut f = (cch * g1 + smu * g2) * fc;
        let mut p = efmu * (half / t2);
        let mut q = zdiv(Complex::from(half), efmu) / t1;

        let mut s1 = f;
        let mut s2 = p;
        let mut ak = one;
        let mut a1 = one;
        let mut ck = Complex::from(one);
        let mut bk = one - dnu2;

        let mut kflag: usize; // 0-based index into cssr/csrr arrays

        if inu > 0 || n > 1 {
            // ── Generate K(DNU,Z) and K(DNU+1,Z) for forward recurrence ──
            // (Fortran label 80-100)
            if caz >= tol {
                let cz = z * z * T::from_f64(0.25);
                let t1_sq = T::from_f64(0.25) * caz * caz;
                loop {
                    f = (f * ak + p + q) / bk;
                    p = p / (ak - dnu);
                    q = q / (ak + dnu);
                    let rak = one / ak;
                    ck = ck * cz * rak;
                    s1 = s1 + ck * f;
                    s2 = s2 + ck * (p - f * ak);
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
            s2 = p2 * rz;
            s1 = s1 * scale;

            if koded == Scaling::Exponential {
                let ez = z.exp();
                s1 = s1 * ez;
                s2 = s2 * ez;
            }
        } else {
            // ── Generate K(FNU,Z), 0 ≤ FNU < 0.5 and N=1 ──
            // (Fortran label 60-70)
            if caz >= tol {
                let cz = z * z * T::from_f64(0.25);
                let t1_sq = T::from_f64(0.25) * caz * caz;
                loop {
                    f = (f * ak + p + q) / bk;
                    p = p / (ak - dnu);
                    q = q / (ak + dnu);
                    let rak = one / ak;
                    ck = ck * cz * rak;
                    s1 = s1 + ck * f;
                    a1 = a1 * t1_sq * rak;
                    bk = bk + ak + ak + one;
                    ak = ak + one;
                    if a1 <= tol {
                        break;
                    }
                }
            }

            y[0] = s1;
            if koded == Scaling::Exponential {
                let ez = z.exp();
                y[0] = s1 * ez;
            }
            return Ok(nz);
        }

        // ── Forward recurrence: label 210 ──
        return forward_recurrence(
            z, fnu, &koded, n, y, &mut nz, tol, elim, alim, iflag, s1, s2, rz, inu, kflag, &cssr,
            &csrr, &bry,
        );
    }

    // ══════════════════════════════════════════════════════════════════
    // MILLER ALGORITHM / ASYMPTOTIC PATH (Fortran label 110+)
    // |z| > R1 or |DNU| == 0.5
    // ══════════════════════════════════════════════════════════════════

    // COEF = RTHPI / sqrt(z) = sqrt(π/2) / sqrt(z) = sqrt(π/(2z))
    let sz = z.sqrt();
    let coef_base = zdiv(Complex::from(T::from_f64(RTHPI)), sz);

    let mut kflag: usize = 1; // KFLAG=2 in Fortran

    let coef;
    if koded == Scaling::Unscaled {
        if z.re > alim {
            // Label 290: scale by exp(z), iflag=1
            koded = Scaling::Exponential;
            iflag = 1;
            kflag = 1; // stays KFLAG=2
            coef = coef_base; // will apply exp later
        } else {
            // Multiply coef by exp(-z) * cssr[kflag]
            coef = coef_base * (-z).exp() * cssr[kflag];
        }
    } else {
        coef = coef_base;
    }

    // Check for DNU = ±0.5 special case (label 300)
    if dnu.abs() == half {
        let s1 = coef;
        let s2 = coef;
        return forward_recurrence(
            z, fnu, &koded, n, y, &mut nz, tol, elim, alim, iflag, s1, s2, rz, inu, kflag, &cssr,
            &csrr, &bry,
        );
    }

    // Check if cos(π·DNU) == 0 or FHS == 0
    let ak_cos = T::from_f64(PI) * dnu;
    let ak_val = ak_cos.cos().abs();
    let fhs = (T::from_f64(0.25) - dnu2).abs();
    if ak_val == zero || fhs == zero {
        // DNU = ±0.5 effective (label 300)
        let s1 = coef;
        let s2 = coef;
        return forward_recurrence(
            z, fnu, &koded, n, y, &mut nz, tol, elim, alim, iflag, s1, s2, rz, inu, kflag, &cssr,
            &csrr, &bry,
        );
    }

    // ── Compute R2 = F(E) for determining backward index K ──
    // T1 = (I1MACH(14) - 1) * D1MACH(5) * 3.321928094
    // D1MACH(5) = log10(2) for binary
    let r1m5 = T::from_f64(R1M5);
    let t1_raw = T::from_f64((T::MACH_DIGITS - 1) as f64) * r1m5 * T::from_f64(3.321928094);
    let t1_clamped = t1_raw.max(T::from_f64(12.0)).min(T::from_f64(60.0));
    let t2_miller = T::from_f64(TTH) * t1_clamped - T::from_f64(6.0);

    let t1_angle = if z.re != zero {
        (z.im / z.re).abs().atan()
    } else {
        T::from_f64(HPI)
    };

    let fk: T;
    let mut fhs_miller = fhs;

    if t2_miller <= caz {
        // ── Forward recurrence to find backward index K (label 150) ──
        let etest = ak_val / (T::from_f64(PI) * caz * tol);
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
                return Err(Error::ConvergenceFailure);
            }
            // Label 160
            fk_val = fk_val + T::from_f64(SPI) * t1_angle * (t2_miller / caz).sqrt();
            fhs_miller = (T::from_f64(0.25) - dnu2).abs();
        }
        fk = fk_val;
    } else {
        // ── Compute backward index K for |z| < R2 (label 170) ──
        let a2 = caz.sqrt();
        let ak_inner = T::from_f64(FPI) * ak_val / (tol * a2.sqrt());
        let three = T::from_f64(3.0);
        let aa = three * t1_angle / (one + caz);
        let bb = T::from_f64(14.7) * t1_angle / (T::from_f64(28.0) + caz);
        let ak_log = (ak_inner.ln() + caz * aa.cos() / (one + T::from_f64(0.008) * caz)) / bb.cos();
        fk = T::from_f64(0.12125) * ak_log * ak_log / caz + T::from_f64(1.5);
    }

    // ── Backward recurrence loop (Miller algorithm, label 180-190) ──
    // Safety: fk is a small positive integer derived from log-based formula
    let k = fk.floor().to_i32().unwrap();
    let fk_int = T::from_f64(k as f64);
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
        p2 = (cb * pt - p1) * ak_inner;
        p1 = pt;
        cs = cs + p2;
        fks = a1 - fk_cur + one;
        fk_cur = fk_cur - one;
    }

    // ── Compute (P2/CS) with better scaling (label after 190) ──
    let tm = zabs(cs);
    let ptr = one / tm;
    let s1_raw = p2 * ptr;
    let cs_conj = cs.conj() * ptr;
    // S1 = COEF * (P2/|CS|) * (conj(CS)/|CS|)
    let s1 = coef * s1_raw * cs_conj;

    if inu == 0 && n == 1 {
        // Direct output path
        if iflag == 1 {
            return handle_iflag1_final(
                z, fnu, n, y, &mut nz, tol, elim, alim, s1, s1, rz, &cssr, &csrr, &bry,
            );
        }
        // Label 240: store and return
        let str_scale = csrr[kflag];
        y[0] = s1 * str_scale;
        return Ok(nz);
    }

    // ── Compute P1/P2 ratio for S2 (label 200) ──
    let tm2 = zabs(p2);
    let ptr2 = one / tm2;
    let p1_scaled = p1 * ptr2;
    let p2_conj = p2.conj() * ptr2;
    let pt_ratio = p1_scaled * p2_conj;

    // S2 = S1 * (DNU + 0.5 - P1/P2·conj(P2)/|P2|) / z + S1... wait
    // Fortran: STR = DNU + 0.5 - PTR; STI = -PTI
    //          CALL ZDIV(STR, STI, ZR, ZI, STR, STI)
    //          STR = STR + 1.0
    //          CALL ZMLT(STR, STI, S1R, S1I, S2R, S2I)
    let ratio = Complex::new(dnu + half - pt_ratio.re, -pt_ratio.im);
    let ratio_div_z = zdiv(ratio, z) + Complex::from(one);
    let s2 = ratio_div_z * s1;

    // ── Forward recurrence (label 210) ──
    forward_recurrence(
        z, fnu, &koded, n, y, &mut nz, tol, elim, alim, iflag, s1, s2, rz, inu, kflag, &cssr,
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
#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
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
) -> Result<usize, Error> {
    let one = T::one();

    // CK = (DNU + 1) * RZ
    let str_init = fnu - T::from_f64(inu as f64) + one;
    let mut ck = rz * str_init;

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
            s2 = ck * st + s1;
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
        return Ok(*nz);
    }
    y[1] = s2 * str_scale;
    if n == 2 {
        return Ok(*nz);
    }

    // ── Label 250-260: additional values via recurrence ──
    let mut kk = 2_usize; // next index to fill (0-based)
    loop {
        if kk >= n {
            return Ok(*nz);
        }

        let p1r = csrr[kflag];
        let ascle = bry[kflag];

        for i in kk..n {
            let p2_val = s2;
            s2 = ck * p2_val + s1;
            s1 = p2_val;
            ck = ck + rz;
            let p2_scaled = s2 * p1r;
            y[i] = p2_scaled;

            if kflag < 2 {
                let p2m = p2_scaled.re.abs().max(p2_scaled.im.abs());
                if p2m > ascle {
                    kflag += 1;
                    s1 = s1 * p1r;
                    s2 = p2_scaled;
                    let str_s = cssr[kflag];
                    s1 = s1 * str_s;
                    s2 = s2 * str_s;
                    kk = i + 1;
                    // Fortran: GO TO 250 (restart outer loop)
                    break;
                }
            }

            if i == n - 1 {
                return Ok(*nz);
            }
        }
        // If we didn't break early, we've filled all values
        if kk >= n {
            return Ok(*nz);
        }
    }
}

/// Handle IFLAG=1 forward recurrence on scaled values (Fortran labels 261-280).
#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
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
) -> Result<usize, Error> {
    let zero = T::zero();
    let czero = Complex::new(zero, zero);
    let half = T::from_f64(0.5);

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
        s2 = ck * st + s1;
        s1 = st;
        ck = ck + rz;

        let az = zabs(s2);
        let alas = az.ln();
        let p2r = -zd.re + alas;

        if p2r >= -elim {
            let s2_log = s2.ln();
            let p2 = s2_log - zd;
            let p1_test = p2.exp() / tol;
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
            zd = zd - elim;
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
                s2 = ck * st + s1;
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
                return Ok(*nz);
            }
            y[1] = s2 * str_scale;
            if n == 2 {
                return Ok(*nz);
            }

            // Continue with label 250-260 recurrence (with KFLAG tracking)
            let t2_val = fnu + T::from_f64(2.0); // FNU + (KK-1) where KK=3
            ck = rz * t2_val;
            let mut kk_idx = 2_usize;
            loop {
                if kk_idx >= n {
                    return Ok(*nz);
                }
                let p1r = csrr[kflag_inner];
                let ascle_loop = bry[kflag_inner];

                for i in kk_idx..n {
                    let p2_val = s2;
                    s2 = ck * p2_val + s1;
                    s1 = p2_val;
                    ck = ck + rz;
                    let p2_scaled = s2 * p1r;
                    y[i] = p2_scaled;

                    if kflag_inner < 2 {
                        let p2m = p2_scaled.re.abs().max(p2_scaled.im.abs());
                        if p2m > ascle_loop {
                            kflag_inner += 1;
                            s1 = s1 * p1r;
                            s2 = p2_scaled;
                            let str_s = cssr[kflag_inner];
                            s1 = s1 * str_s;
                            s2 = s2 * str_s;
                            kk_idx = i + 1;
                            break;
                        }
                    }

                    if i == n - 1 {
                        return Ok(*nz);
                    }
                }
                if kk_idx >= n {
                    return Ok(*nz);
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
            return Ok(*nz);
        }
        y[1] = s2 * str_scale;
        if n == 2 {
            return Ok(*nz);
        }
        return Ok(*nz);
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
#[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
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
) -> Result<usize, Error> {
    y[0] = s1;
    if n > 1 {
        y[1] = s2;
    }

    // Label 280: call ZKSCL
    let ascle = bry[0];
    *nz = zkscl(zd, fnu, y, rz, ascle, tol, elim);
    let inu_remaining = n as i32 - *nz as i32;
    if inu_remaining <= 0 {
        return Ok(*nz);
    }

    let kk = *nz; // 0-based index of first non-zero value
    let mut s1_local = y[kk];
    y[kk] = s1_local * csrr[0];
    if inu_remaining == 1 {
        return Ok(*nz);
    }

    let kk2 = *nz + 1;
    let mut s2_local = y[kk2];
    y[kk2] = s2_local * csrr[0];
    if inu_remaining == 2 {
        return Ok(*nz);
    }

    // Continue recurrence for remaining values (label 250-260 with KFLAG tracking)
    let t2_val = fnu + T::from_f64(kk2 as f64);
    let mut ck = rz * t2_val;
    let mut kflag: usize = 0; // KFLAG=1 in Fortran

    let mut kk_idx = kk2 + 1; // 0-based next index to fill
    loop {
        if kk_idx >= n {
            return Ok(*nz);
        }
        let p1r = csrr[kflag];
        let ascle = bry[kflag];

        for i in kk_idx..n {
            let p2_val = s2_local;
            s2_local = ck * p2_val + s1_local;
            s1_local = p2_val;
            ck = ck + rz;
            let p2_scaled = s2_local * p1r;
            y[i] = p2_scaled;

            if kflag < 2 {
                let p2m = p2_scaled.re.abs().max(p2_scaled.im.abs());
                if p2m > ascle {
                    kflag += 1;
                    s1_local = s1_local * p1r;
                    s2_local = p2_scaled;
                    let str_s = cssr[kflag];
                    s1_local = s1_local * str_s;
                    s2_local = s2_local * str_s;
                    kk_idx = i + 1;
                    break;
                }
            }

            if i == n - 1 {
                return Ok(*nz);
            }
        }
        if kk_idx >= n {
            return Ok(*nz);
        }
    }
}
