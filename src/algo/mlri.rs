//! Miller algorithm for I Bessel function, normalized by Neumann series.
//!
//! Translation of Fortran ZMLRI from TOMS 644 (zbsubs.f lines 3322-3526).

use num_complex::Complex;

use crate::algo::gamln::gamln;
use crate::machine::BesselFloat;
use crate::types::Scaling;
use crate::utils::{mul_add, mul_add_scalar, reciprocal_z, zabs};

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
    // Note: caller (zbinu) already zeroes the output buffer.

    let scle = T::MACH_TINY / tol;
    let nz: i32 = 0;
    let az = zabs(z);
    // Safety: az (modulus) is finite and bounded
    let iaz = az.to_i32().unwrap();
    // Safety: fnu is finite and < ~1e15 per upper-interface checks
    let ifnu = fnu.to_i32().unwrap();
    let inu = ifnu + n as i32 - 1;
    let at = T::from_f64(iaz as f64 + 1.0);
    let raz = one / az;
    let rz = reciprocal_z(z);
    let mut ck = rz * (at * T::from_f64(0.5));
    let mut p1 = czero;
    let mut p2 = Complex::from(one);
    let ack = (at + one) * raz;
    let rho = ack + (ack * ack - one).sqrt();
    let rho2 = rho * rho;
    let mut tst = (rho2 + rho2) / ((rho2 - one) * (rho - one));
    tst = tst / tol;

    // Compute relative truncation error index for series (Fortran lines 3367-3381)
    let mut ak = at;
    let mut i_val: i32 = 0;
    for i in 1..=80 {
        let pt = p2;
        p2 = mul_add(-ck, pt, p1);
        p1 = pt;
        ck = ck + rz;
        let ap = zabs(p2);
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
        p1 = czero;
        p2 = Complex::from(one);
        let at2 = T::from_f64(inu as f64 + 1.0);
        ck = z.conj() * (at2 * raz * raz);
        let ack2 = at2 * raz;
        tst = (ack2 / tol).sqrt();
        let mut refined = false;

        let mut found = false;
        for kk in 1..=80 {
            let pt = p2;
            p2 = mul_add(-ck, pt, p1);
            p1 = pt;
            ck = ck + rz;
            let ap = zabs(p2);
            if ap >= tst {
                if refined {
                    k = kk;
                    found = true;
                    break;
                }
                let ack3 = zabs(ck);
                let flam = ack3 + (ack3 * ack3 - one).sqrt();
                let fkap = ap / zabs(p1);
                let rho_new = flam.min(fkap);
                tst = tst * (rho_new / (rho_new * rho_new - one)).sqrt();
                refined = true;
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
    let mut fkk = T::from_f64(kk as f64);
    p1 = czero;
    p2 = Complex::new(scle, zero);
    let fnf = fnu - T::from_f64(ifnu as f64);
    let tfnf = fnf + fnf;
    // Safety: all gamln arguments are > 0 guaranteed by algorithm invariant
    let bk_ln =
        gamln(fkk + tfnf + one).unwrap() - gamln(fkk + one).unwrap() - gamln(tfnf + one).unwrap();
    let mut bk = bk_ln.exp();
    let mut sum = czero;
    let km = kk - inu;

    // Backward recurrence: first phase (kk down to inu+1) — Fortran DO 50
    for _i in 1..=km {
        let pt = p2;
        p2 = mul_add_scalar(rz * pt, fkk + fnf, p1);
        p1 = pt;
        ak = one - tfnf / (fkk + tfnf);
        let ack_val = bk * ak;
        sum = mul_add_scalar(p1, ack_val + bk, sum);
        bk = ack_val;
        fkk = fkk - one;
    }

    // Store Y(N) (Fortran line 3457-3458)
    y[n - 1] = p2;

    // Continue backward recurrence storing Y values (Fortran DO 60)
    if n > 1 {
        for i in 2..=n {
            let pt = p2;
            p2 = p1 + rz * pt * (fkk + fnf);
            p1 = pt;
            ak = one - tfnf / (fkk + tfnf);
            let ack_val = bk * ak;
            sum = sum + p1 * (ack_val + bk);
            bk = ack_val;
            fkk = fkk - one;
            // M = N - I + 1 (Fortran 1-based) → 0-based: n - i
            let m = n - i;
            y[m] = p2;
        }
    }

    // Continue backward recurrence for normalization (Fortran DO 80)
    if ifnu > 0 {
        for _i in 1..=ifnu {
            let pt = p2;
            p2 = p1 + rz * pt * (fkk + fnf);
            p1 = pt;
            ak = one - tfnf / (fkk + tfnf);
            let ack_val = bk * ak;
            sum = sum + p1 * (ack_val + bk);
            bk = ack_val;
            fkk = fkk - one;
        }
    }

    // Label 90: compute normalization constant (Fortran lines 3494-3521)
    let pt = if kode == Scaling::Exponential {
        Complex::new(zero, z.im)
    } else {
        z
    };
    let rz_ln = rz.ln();
    let p1_norm = pt - rz_ln * fnf;
    // Safety: 1 + fnf > 0 guaranteed (fnf is fractional part of fnu)
    let ap = gamln(one + fnf).unwrap();
    let pt_norm = Complex::new(p1_norm.re - ap, p1_norm.im);

    // Division cexp(pt)/(sum+p2) avoiding overflow (Fortran lines 3507-3521)
    p2 = p2 + sum;
    let ap2 = zabs(p2);
    let p1_inv = one / ap2;
    let ck_n = pt_norm.exp() * p1_inv;
    let pt_n = p2.conj() * p1_inv;
    // zmlt(ck, pt) → cnorm
    let cnorm = ck_n * pt_n;

    // Normalize all Y values (Fortran DO 100)
    for y_item in y.iter_mut().take(n) {
        *y_item = *y_item * cnorm;
    }

    nz
}
