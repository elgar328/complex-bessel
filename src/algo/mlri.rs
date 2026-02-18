//! Miller algorithm for I Bessel function, normalized by Neumann series.
//!
//! Translation of Fortran ZMLRI from TOMS 644 (zbsubs.f lines 3322-3526).

#![allow(clippy::excessive_precision)]

use num_complex::Complex;

use crate::algo::gamln::gamln;
use crate::machine::BesselFloat;
use crate::types::Scaling;
use crate::utils::zabs;

/// Miller algorithm for I Bessel function.
///
/// Writes results into `y` and returns `nz` where nz = -2 indicates convergence failure.
pub(crate) fn zmlri<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    y: &mut [Complex<T>],
    tol: T,
) -> i32 {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);

    let n = y.len();
    for v in y.iter_mut() {
        *v = czero;
    }

    let scle = T::MACH_TINY / tol;
    let nz: i32 = 0;
    let az = zabs(z);
    let iaz = az.to_i32().unwrap();
    let ifnu = fnu.to_i32().unwrap();
    let inu = ifnu + n as i32 - 1;
    let at = T::from(iaz as f64 + 1.0).unwrap();
    let raz = one / az;
    let str = z.re * raz;
    let sti = -z.im * raz;
    let mut ckr = str * at * raz;
    let mut cki = sti * at * raz;
    let rzr = (str + str) * raz;
    let rzi = (sti + sti) * raz;
    let mut p1r = zero;
    let mut p1i = zero;
    let mut p2r = one;
    let mut p2i = zero;
    let ack = (at + one) * raz;
    let rho = ack + (ack * ack - one).sqrt();
    let rho2 = rho * rho;
    let mut tst = (rho2 + rho2) / ((rho2 - one) * (rho - one));
    tst = tst / tol;

    // Compute relative truncation error index for series (Fortran lines 3367-3381)
    let mut ak = at;
    let mut i_val: i32 = 0;
    for i in 1..=80 {
        let ptr = p2r;
        let pti = p2i;
        p2r = p1r - (ckr * ptr - cki * pti);
        p2i = p1i - (cki * ptr + ckr * pti);
        p1r = ptr;
        p1i = pti;
        ckr = ckr + rzr;
        cki = cki + rzi;
        let ap = zabs(Complex::new(p2r, p2i));
        if ap > tst * ak * ak {
            i_val = i;
            break;
        }
        ak = ak + one;
    }
    if i_val == 0 {
        // Convergence failure (Fortran label 110)
        return -2;
    }

    i_val += 1; // Fortran: I = I + 1 (line 3383)
    let mut k: i32 = 0;

    if inu >= iaz {
        // Compute relative truncation error for ratios (Fortran lines 3388-3420)
        p1r = zero;
        p1i = zero;
        p2r = one;
        p2i = zero;
        let at2 = T::from(inu as f64 + 1.0).unwrap();
        let str2 = z.re * raz;
        let sti2 = -z.im * raz;
        ckr = str2 * at2 * raz;
        cki = sti2 * at2 * raz;
        let ack2 = at2 * raz;
        tst = (ack2 / tol).sqrt();
        let mut itime = 1;

        let mut found = false;
        for kk in 1..=80 {
            let ptr = p2r;
            let pti = p2i;
            p2r = p1r - (ckr * ptr - cki * pti);
            p2i = p1i - (ckr * pti + cki * ptr);
            p1r = ptr;
            p1i = pti;
            ckr = ckr + rzr;
            cki = cki + rzi;
            let ap = zabs(Complex::new(p2r, p2i));
            if ap >= tst {
                if itime == 2 {
                    k = kk;
                    found = true;
                    break;
                }
                let ack3 = zabs(Complex::new(ckr, cki));
                let flam = ack3 + (ack3 * ack3 - one).sqrt();
                let fkap = ap / zabs(Complex::new(p1r, p1i));
                let rho_new = flam.min(fkap);
                tst = tst * (rho_new / (rho_new * rho_new - one)).sqrt();
                itime = 2;
            }
        }
        if !found && k == 0 {
            // Convergence failure (Fortran label 110)
            return -2;
        }
    }

    // Label 40: backward recurrence and sum normalizing relation
    // (Fortran lines 3421-3526)
    k += 1;
    let kk = (i_val + iaz).max(k + inu);
    let mut fkk = T::from(kk as f64).unwrap();
    p1r = zero;
    p1i = zero;
    p2r = scle;
    p2i = zero;
    let fnf = fnu - T::from(ifnu as f64).unwrap();
    let tfnf = fnf + fnf;
    let bk_ln =
        gamln(fkk + tfnf + one).unwrap() - gamln(fkk + one).unwrap() - gamln(tfnf + one).unwrap();
    let mut bk = bk_ln.exp();
    let mut sumr = zero;
    let mut sumi = zero;
    let km = kk - inu;

    // Backward recurrence: first phase (kk down to inu+1) — Fortran DO 50
    for _i in 1..=km {
        let ptr = p2r;
        let pti = p2i;
        p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
        p2i = p1i + (fkk + fnf) * (rzi * ptr + rzr * pti);
        p1r = ptr;
        p1i = pti;
        ak = one - tfnf / (fkk + tfnf);
        let ack_val = bk * ak;
        sumr = sumr + (ack_val + bk) * p1r;
        sumi = sumi + (ack_val + bk) * p1i;
        bk = ack_val;
        fkk = fkk - one;
    }

    // Store Y(N) (Fortran line 3457-3458)
    y[n - 1] = Complex::new(p2r, p2i);

    // Continue backward recurrence storing Y values (Fortran DO 60)
    if n > 1 {
        for i in 2..=n {
            let ptr = p2r;
            let pti = p2i;
            p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
            p2i = p1i + (fkk + fnf) * (rzi * ptr + rzr * pti);
            p1r = ptr;
            p1i = pti;
            ak = one - tfnf / (fkk + tfnf);
            let ack_val = bk * ak;
            sumr = sumr + (ack_val + bk) * p1r;
            sumi = sumi + (ack_val + bk) * p1i;
            bk = ack_val;
            fkk = fkk - one;
            // M = N - I + 1 (Fortran 1-based) → 0-based: n - i
            let m = n - i;
            y[m] = Complex::new(p2r, p2i);
        }
    }

    // Continue backward recurrence for normalization (Fortran DO 80)
    if ifnu > 0 {
        for _i in 1..=ifnu {
            let ptr = p2r;
            let pti = p2i;
            p2r = p1r + (fkk + fnf) * (rzr * ptr - rzi * pti);
            p2i = p1i + (fkk + fnf) * (rzr * pti + rzi * ptr);
            p1r = ptr;
            p1i = pti;
            ak = one - tfnf / (fkk + tfnf);
            let ack_val = bk * ak;
            sumr = sumr + (ack_val + bk) * p1r;
            sumi = sumi + (ack_val + bk) * p1i;
            bk = ack_val;
            fkk = fkk - one;
        }
    }

    // Label 90: compute normalization constant (Fortran lines 3494-3521)
    let mut ptr = z.re;
    let pti = z.im;
    if kode == Scaling::Exponential {
        ptr = zero;
    }
    // zlog(rz)
    let rz_abs = zabs(Complex::new(rzr, rzi));
    let str_log = rz_abs.ln();
    let sti_log = rzi.atan2(rzr);
    let p1r_norm = -fnf * str_log + ptr;
    let p1i_norm = -fnf * sti_log + pti;
    let ap = gamln(one + fnf).unwrap();
    let ptr_norm = p1r_norm - ap;
    let pti_norm = p1i_norm;

    // Division cexp(pt)/(sum+p2) avoiding overflow (Fortran lines 3507-3521)
    p2r = p2r + sumr;
    p2i = p2i + sumi;
    let ap2 = zabs(Complex::new(p2r, p2i));
    let p1_inv = one / ap2;
    // zexp(ptr_norm, pti_norm)
    let exp_r = ptr_norm.exp();
    let str_e = exp_r * pti_norm.cos();
    let sti_e = exp_r * pti_norm.sin();
    let ck_nr = str_e * p1_inv;
    let ck_ni = sti_e * p1_inv;
    let pt_nr = p2r * p1_inv;
    let pt_ni = -p2i * p1_inv;
    // zmlt(ck, pt) → cnorm
    let cnormr = ck_nr * pt_nr - ck_ni * pt_ni;
    let cnormi = ck_nr * pt_ni + ck_ni * pt_nr;

    // Normalize all Y values (Fortran DO 100)
    for y_item in y.iter_mut().take(n) {
        let str_val = y_item.re * cnormr - y_item.im * cnormi;
        let sti_val = y_item.re * cnormi + y_item.im * cnormr;
        *y_item = Complex::new(str_val, sti_val);
    }

    nz
}
