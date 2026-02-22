//! Airy functions Ai(z), Bi(z) and their derivatives.
//!
//! Translation of Fortran ZAIRY (zbsubs.f lines 1462-1856)
//! and ZBIRY (zbsubs.f lines 1857-2222) from TOMS 644.

#![allow(clippy::excessive_precision)]

use num_complex::Complex;

use crate::algo::acai::zacai;
use crate::algo::binu::zbinu;
use crate::algo::bknu::zbknu;
use crate::machine::BesselFloat;
use crate::types::{Accuracy, AiryDerivative, Error, Scaling};
use crate::utils::{zabs, zdiv};

// ZAIRY constants (Fortran ZAIRY DATA, lines 1600-1603)
use crate::algo::constants::{PI, TTH};
const AI_C1: f64 = 3.55028053887817240e-01;
const AI_C2: f64 = 2.58819403792806799e-01;
const AI_COEF: f64 = 1.83776298473930683e-01; // 1/(pi*sqrt(3))

// ZBIRY constants (Fortran ZBIRY DATA, lines 1989-1992)
const BI_C1: f64 = 6.14926627446000736e-01; // sqrt(3) * AI_C1
const BI_C2: f64 = 4.48288357353826359e-01; // sqrt(3) * AI_C2
const BI_COEF: f64 = 5.77350269189625765e-01; // 1/sqrt(3)

// ─── Shared power series helper ─────────────────────────────────────────────

/// Compute 25-term power series partial sums (s1, s2) for |z| <= 1.
///
/// Shared by ZAIRY and ZBIRY. The caller combines s1, s2 with different
/// constants/signs for Ai vs Bi.
///
/// Caller must guarantee `az >= tol` (the `az < tol` case is handled separately).
fn airy_power_series<T: BesselFloat>(
    z: Complex<T>,
    az: T,
    fid: T,
    tol: T,
) -> (Complex<T>, Complex<T>) {
    let one = T::one();

    let mut s1 = Complex::from(one);
    let mut s2 = Complex::from(one);

    let aa = az * az;
    if aa < tol / az {
        // Series terms too small; s1=s2=1 is sufficient (Fortran label 40)
        return (s1, s2);
    }

    let mut trm1 = Complex::from(one);
    let mut trm2 = Complex::from(one);
    let mut atrm = one;

    // z³ (Fortran lines 1629-1633)
    let z2 = z * z;
    let z3 = z2 * z;
    let az3 = az * aa;

    // Divisor coefficients (Fortran lines 1634-1642)
    let two = T::from_f64(2.0);
    let three = T::from_f64(3.0);
    let four = T::from_f64(4.0);
    let nine = T::from_f64(9.0);
    let eighteen = T::from_f64(18.0);
    let twenty_four = T::from_f64(24.0);
    let thirty = T::from_f64(30.0);

    let mut ak = two + fid;
    let mut bk = three - fid - fid;
    let ck = four - fid;
    let dk = three + fid + fid;
    let mut d1 = ak * dk;
    let mut d2 = bk * ck;
    let mut ad = d1.min(d2);
    ak = twenty_four + nine * fid;
    bk = thirty - nine * fid;

    for _ in 0..25 {
        trm1 = trm1 * z3 / d1;
        s1 = s1 + trm1;

        trm2 = trm2 * z3 / d2;
        s2 = s2 + trm2;

        atrm = atrm * az3 / ad;
        d1 = d1 + ak;
        d2 = d2 + bk;
        ad = d1.min(d2);
        if atrm < tol * ad {
            break;
        }
        ak = ak + eighteen;
        bk = bk + eighteen;
    }

    (s1, s2)
}

// ─── ZAIRY ──────────────────────────────────────────────────────────────────

/// Compute Airy function Ai(z) or its derivative Ai'(z).
///
/// Translation of Fortran ZAIRY (zbsubs.f lines 1462-1856).
///
/// Returns `(result, nz, status)` where nz=0 normal, nz=1 underflow (result=0).
/// `status` is `Accuracy::Reduced` when |z| > sqrt(AA^(2/3))
/// (Fortran IERR=3: more than half of significant digits may be lost).
#[inline]
pub(crate) fn zairy<T: BesselFloat>(
    z: Complex<T>,
    id: AiryDerivative,
    kode: Scaling,
) -> Result<(Complex<T>, i32, Accuracy), Error> {
    let zero = T::zero();
    let one = T::one();
    let tth = T::from_f64(TTH);
    let c1 = T::from_f64(AI_C1);
    let c2 = T::from_f64(AI_C2);
    let coef = T::from_f64(AI_COEF);
    let czero = Complex::new(zero, zero);

    let az = zabs(z);
    let tol = T::tol();
    let fid = if id == AiryDerivative::Derivative {
        one
    } else {
        zero
    };

    if az > one {
        // |z| > 1: K function path (Fortran label 70)
        return zairy_large_z(z, az, id, kode, fid, tol, tth, coef);
    }

    // ── |z| <= 1: power series (Fortran lines 1614-1693) ──

    if az < tol {
        // Tiny z (Fortran label 170, lines 1815-1837)
        let aa = T::from_f64(1.0e3) * T::MACH_TINY;

        if id == AiryDerivative::Value {
            // Ai(z) ≈ C1 - C2*z (Fortran lines 1819-1826)
            let mut s1 = czero;
            if az > aa {
                s1 = z * c2;
            }
            let ai = Complex::new(c1, zero) - s1;
            return if kode == Scaling::Exponential {
                Ok((zairy_scale_exp_zta(z, ai, tth), 0, Accuracy::Normal))
            } else {
                Ok((ai, 0, Accuracy::Normal))
            };
        }
        // Ai'(z) ≈ -C2 + C1*z²/2 (Fortran lines 1827-1837)
        let mut ai = Complex::new(-c2, zero);
        let aa_sqrt = aa.sqrt();
        if az > aa_sqrt {
            let s1 = z * z * T::from_f64(0.5);
            ai = ai + s1 * c1;
        }
        return if kode == Scaling::Exponential {
            Ok((zairy_scale_exp_zta(z, ai, tth), 0, Accuracy::Normal))
        } else {
            Ok((ai, 0, Accuracy::Normal))
        };
    }

    // Non-tiny series (Fortran lines 1622-1693)
    let (s1, s2) = airy_power_series(z, az, fid, tol);

    if id == AiryDerivative::Value {
        // Ai = C1*S1 - C2*(z*S2) (Fortran lines 1664-1665)
        let ai = s1 * c1 - z * s2 * c2;
        if kode == Scaling::Unscaled {
            return Ok((ai, 0, Accuracy::Normal));
        }
        Ok((zairy_scale_exp_zta(z, ai, tth), 0, Accuracy::Normal))
    } else {
        // Ai' = -C2*S2 (Fortran lines 1676-1677)
        let mut ai = -(s2 * c2);
        if az > tol {
            // Ai' += C1/(1+fid) * z² * S1 (Fortran lines 1679-1683)
            let cc = c1 / (one + fid);
            ai = ai + z * s1 * z * cc;
        }
        if kode == Scaling::Unscaled {
            return Ok((ai, 0, Accuracy::Normal));
        }
        Ok((zairy_scale_exp_zta(z, ai, tth), 0, Accuracy::Normal))
    }
}

/// Multiply `val` by exp(zta) where zta = (2/3)*z*sqrt(z).
/// KODE=2 scaling for ZAIRY power series branch (Fortran lines 1667-1673).
#[inline]
fn zairy_scale_exp_zta<T: BesselFloat>(z: Complex<T>, val: Complex<T>, tth: T) -> Complex<T> {
    let zta = z * z.sqrt() * tth;
    val * zta.exp()
}

/// Compute zta = (2/3)*z*sqrt(z) with branch cut correction.
/// Shared by zairy_large_z and zbiry_large_z.
/// (Fortran ZAIRY lines 1731-1749, ZBIRY lines 2123-2139)
#[inline]
fn airy_zta<T: BesselFloat>(z: Complex<T>, tth: T) -> (Complex<T>, Complex<T>) {
    let zero = T::zero();
    let csq = z.sqrt();
    let mut zta = z * csq * tth;
    let ak = zta.im;
    if z.re < zero {
        zta = Complex::new(-zta.re.abs(), ak);
    }
    if z.im == zero && z.re <= zero {
        zta = Complex::new(zero, ak);
    }
    (csq, zta)
}

/// ZAIRY |z| > 1 branch: compute Ai via K Bessel function.
/// Fortran lines 1697-1856.
#[allow(clippy::too_many_arguments)]
fn zairy_large_z<T: BesselFloat>(
    z: Complex<T>,
    az: T,
    id: AiryDerivative,
    kode: Scaling,
    fid: T,
    tol: T,
    tth: T,
    coef: T,
) -> Result<(Complex<T>, i32, Accuracy), Error> {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);

    let fnu = (one + fid) / T::from_f64(3.0);

    // Machine constants (Fortran lines 1709-1720)
    let elim = T::elim();
    let alim = T::alim();
    let rl = T::rl();
    let alaz = az.ln();

    let mut nz: i32 = 0;

    // Range test (Fortran lines 1724-1728)
    let half = T::from_f64(0.5);
    // Fortran: 0.5D0/TOL is capped at R1M5*(K1-1) where K1=I1MACH(15)
    // For IEEE 754 f64: R1M5*(1023) = 307.95..., AA = 10^307.95 ≈ 8.9e307
    // The literal 1073741823.5 = 2^30 - 0.5, used as a machine-safe upper bound
    // (Fortran ZAIRY line 1724 / ZBIRY line 2116)
    let mut aa = (half / tol).min(T::from_f64(1073741823.5));
    aa = aa.powf(tth);
    if az > aa {
        return Err(Error::TotalPrecisionLoss);
    }
    // IERR=3 precision warning (Fortran ZAIRY lines 1728-1729):
    // |z| > sqrt(AA^(2/3)) means more than half of significant digits lost
    let aa_sqrt = aa.sqrt();
    let status = if az > aa_sqrt {
        Accuracy::Reduced
    } else {
        Accuracy::Normal
    };

    // ZTA = (2/3)*z*sqrt(z) with branch cut correction (Fortran lines 1731-1749)
    let (csq, zta) = airy_zta(z, tth);
    let mut iflag: i32 = 0;
    let mut sfac = one;

    // Dispatch based on Re(zta) and Re(z) (Fortran line 1752)
    let aa_zta = zta.re;
    if aa_zta >= zero && z.re > zero {
        // Label 110: direct ZBKNU path (right half plane)
        // Underflow test (Fortran lines 1774-1782)
        if kode != Scaling::Exponential && aa_zta >= alim {
            let test_val = -aa_zta - T::from_f64(0.25) * alaz;
            iflag = 2;
            sfac = one / tol;
            if test_val < -elim {
                // Underflow: return zero (Fortran label 210)
                return Ok((czero, 1, status));
            }
        }

        // Label 120: ZBKNU (Fortran lines 1784-1785)
        let mut cy_buf = [czero];
        let nz_k = zbknu(zta, fnu, kode, &mut cy_buf, tol, elim, alim)?;
        nz += nz_k as i32;

        let s1 = cy_buf[0] * coef;
        return zairy_form_result(z, csq, s1, id, iflag, sfac, nz, status);
    }

    // Re(zta) < 0 or Re(z) <= 0: analytic continuation via ZACAI
    if kode != Scaling::Exponential {
        // Overflow test (Fortran lines 1757-1761)
        if -aa_zta > alim {
            let test_val = -aa_zta + T::from_f64(0.25) * alaz;
            iflag = 1;
            sfac = tol;
            if test_val > elim {
                return Err(Error::Overflow);
            }
        }
    }

    // Label 100: ZACAI call (Fortran lines 1766-1772)
    let mr: i32 = if z.im < zero { -1 } else { 1 };
    let mut cy_buf = [czero];
    let nn = zacai(zta, fnu, kode, mr, &mut cy_buf, rl, tol, elim, alim)?;
    if nn < 0 {
        // Fortran label 280
        return if nn == -1 {
            Err(Error::Overflow)
        } else {
            Err(Error::ConvergenceFailure)
        };
    }
    nz += nn;

    let s1 = cy_buf[0] * coef;
    zairy_form_result(z, csq, s1, id, iflag, sfac, nz, status)
}

/// Form final Airy result from s1 = cy[0]*coef.
/// Handles IFLAG scaling and ID dispatch (Fortran lines 1787-1813).
#[allow(clippy::too_many_arguments)]
fn zairy_form_result<T: BesselFloat>(
    z: Complex<T>,
    csq: Complex<T>,
    s1: Complex<T>,
    id: AiryDerivative,
    iflag: i32,
    sfac: T,
    nz: i32,
    status: Accuracy,
) -> Result<(Complex<T>, i32, Accuracy), Error> {
    if iflag == 0 {
        // Normal case
        if id == AiryDerivative::Value {
            // Ai = sqrt(z) * s1 (Fortran lines 1791-1792)
            Ok((csq * s1, nz, status))
        } else {
            // Ai' = -(z * s1) (Fortran lines 1795-1796)
            Ok((-(z * s1), nz, status))
        }
    } else {
        // IFLAG != 0: scaled computation (Fortran lines 1798-1813)
        let s1s = s1 * sfac;
        if id == AiryDerivative::Value {
            Ok((s1s * csq / sfac, nz, status))
        } else {
            Ok((-(s1s * z) / sfac, nz, status))
        }
    }
}

// ─── ZBIRY ──────────────────────────────────────────────────────────────────

/// Compute Airy function Bi(z) or its derivative Bi'(z).
///
/// Translation of Fortran ZBIRY (zbsubs.f lines 1857-2222).
///
/// KODE=2 returns exp(-|Re(zta)|)*Bi(z) where zta=(2/3)*z*sqrt(z).
///
/// Returns `(result, status)` where `status` is `Accuracy::Reduced`
/// when |z| > sqrt(AA^(2/3)) (Fortran IERR=3).
#[inline]
pub(crate) fn zbiry<T: BesselFloat>(
    z: Complex<T>,
    id: AiryDerivative,
    kode: Scaling,
) -> Result<(Complex<T>, Accuracy), Error> {
    let zero = T::zero();
    let one = T::one();
    let tth = T::from_f64(TTH);
    let c1 = T::from_f64(BI_C1);
    let c2 = T::from_f64(BI_C2);
    let coef = T::from_f64(BI_COEF);
    let pi_t = T::from_f64(PI);

    let az = zabs(z);
    let tol = T::tol();
    let fid = if id == AiryDerivative::Derivative {
        one
    } else {
        zero
    };

    if az > one {
        // |z| > 1: I function path (Fortran label 70)
        return zbiry_large_z(z, az, id, kode, fid, tol, tth, coef, pi_t);
    }

    // ── |z| <= 1: power series (Fortran lines 2002-2084) ──

    if az < tol {
        // Tiny z (Fortran label 130, lines 2204-2208)
        // Bi = C1*(1-fid) + fid*C2 (real-valued)
        let bi = Complex::new(c1 * (one - fid) + fid * c2, zero);
        return Ok((bi, Accuracy::Normal));
    }

    let (s1, s2) = airy_power_series(z, az, fid, tol);

    if id == AiryDerivative::Value {
        // Bi = C1*S1 + C2*(z*S2) (Fortran lines 2053-2054, note + sign)
        let bi = s1 * c1 + z * s2 * c2;
        if kode == Scaling::Unscaled {
            return Ok((bi, Accuracy::Normal));
        }
        // KODE=2: multiply by exp(-|Re(zta)|) (Fortran lines 2056-2063)
        Ok((zbiry_scale_exp(z, bi, tth), Accuracy::Normal))
    } else {
        // Bi' = C2*S2 (Fortran lines 2066-2067, no minus sign)
        let mut bi = s2 * c2;
        if az > tol {
            // Bi' += C1/(1+fid) * z² * S1 (Fortran lines 2069-2073)
            let cc = c1 / (one + fid);
            bi = bi + z * s1 * z * cc;
        }
        if kode == Scaling::Unscaled {
            return Ok((bi, Accuracy::Normal));
        }
        Ok((zbiry_scale_exp(z, bi, tth), Accuracy::Normal))
    }
}

/// Multiply `val` by exp(-|Re(zta)|) for ZBIRY KODE=2 in the power series branch.
/// (Fortran lines 2056-2063 / 2076-2083)
#[inline]
fn zbiry_scale_exp<T: BesselFloat>(z: Complex<T>, val: Complex<T>, tth: T) -> Complex<T> {
    let zta_re = (z * z.sqrt() * tth).re;
    val * (-zta_re.abs()).exp()
}

/// ZBIRY |z| > 1 branch: compute Bi via I Bessel function.
/// Fortran lines 2088-2221.
#[allow(clippy::too_many_arguments)]
fn zbiry_large_z<T: BesselFloat>(
    z: Complex<T>,
    az: T,
    id: AiryDerivative,
    kode: Scaling,
    fid: T,
    tol: T,
    tth: T,
    coef: T,
    pi_t: T,
) -> Result<(Complex<T>, Accuracy), Error> {
    let zero = T::zero();
    let one = T::one();
    let two = T::from_f64(2.0);
    let three = T::from_f64(3.0);
    let half = T::from_f64(0.5);

    let fnu = (one + fid) / three;

    // Machine constants (Fortran lines 2101-2112)
    let elim = T::elim();
    let alim = T::alim();
    let rl = T::rl();
    let fnul = T::fnul();

    // Range test (Fortran lines 2116-2122)
    // Fortran: 0.5D0/TOL is capped at R1M5*(K1-1) where K1=I1MACH(15)
    // For IEEE 754 f64: R1M5*(1023) = 307.95..., AA = 10^307.95 ≈ 8.9e307
    // The literal 1073741823.5 = 2^30 - 0.5, used as a machine-safe upper bound
    // (Fortran ZAIRY line 1724 / ZBIRY line 2116)
    let mut aa = (half / tol).min(T::from_f64(1073741823.5));
    aa = aa.powf(tth);
    if az > aa {
        return Err(Error::TotalPrecisionLoss);
    }
    // IERR=3 precision warning (Fortran ZBIRY lines 2118-2119):
    // |z| > sqrt(AA^(2/3)) means more than half of significant digits lost
    let aa_sqrt = aa.sqrt();
    let status = if az > aa_sqrt {
        Accuracy::Reduced
    } else {
        Accuracy::Normal
    };

    // ZTA = (2/3)*z*sqrt(z) with branch cut correction (Fortran lines 2123-2139)
    let (csq, mut zta) = airy_zta(z, tth);
    let mut sfac = one;

    // Overflow test (Fortran lines 2141-2150)
    let aa_zta = zta.re;
    if kode != Scaling::Exponential {
        let bb_abs = aa_zta.abs();
        if bb_abs >= alim {
            let bb_test = bb_abs + T::from_f64(0.25) * az.ln();
            sfac = tol;
            if bb_test > elim {
                return Err(Error::Overflow);
            }
        }
    }

    // FMR determination (Fortran lines 2152-2157)
    let mut fmr = zero;
    if !(aa_zta >= zero && z.re > zero) {
        fmr = pi_t;
        if z.im < zero {
            fmr = -pi_t;
        }
        zta = -zta;
    }

    // First ZBINU: I_{fnu}(zta) (Fortran lines 2163-2164)
    let czero = Complex::new(zero, zero);
    let mut cy1_buf = [czero];
    zbinu(zta, fnu, kode, &mut cy1_buf, rl, fnul, tol, elim, alim)?;

    // S1 = exp(i*FMR*FNU) * CY(1) * SFAC (Fortran lines 2166-2171)
    let aa_fmr = fmr * fnu;
    let phase = Complex::new(aa_fmr.cos(), aa_fmr.sin());
    let mut s1 = phase * cy1_buf[0] * sfac;

    // Second ZBINU: I_{fnu2}(zta), I_{fnu2+1}(zta) (Fortran lines 2172-2178)
    let fnu2 = (two - fid) / three;
    let mut cy2_buf = [czero; 2];
    zbinu(zta, fnu2, kode, &mut cy2_buf, rl, fnul, tol, elim, alim)?;
    let cy2_0 = cy2_buf[0] * sfac;
    let cy2_1 = cy2_buf[1] * sfac;

    // Backward recurrence: I_{fnu2-1} = 2*fnu2*(CY(1)/ZTA) + CY(2)
    // (Fortran lines 2182-2184)
    let s2 = zdiv(cy2_0, zta) * (fnu2 + fnu2) + cy2_1;

    // Final combination (Fortran lines 2185-2189)
    let aa_rot = fmr * (fnu2 - one);
    let phase2 = Complex::new(aa_rot.cos(), aa_rot.sin());
    s1 = (s1 + phase2 * s2) * coef;

    // Form result (Fortran lines 2190-2202)
    if id == AiryDerivative::Value {
        // Bi = sqrt(z) * S1 / SFAC (Fortran lines 2191-2195)
        Ok((csq * s1 / sfac, status))
    } else {
        // Bi' = z * S1 / SFAC (Fortran lines 2198-2202)
        Ok((z * s1 / sfac, status))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    // ── ZAIRY tests ──

    #[test]
    fn zairy_zero() {
        // Ai(0) = C1 ≈ 0.35502805388781724
        let (ai, nz, _status) = zairy(
            Complex64::new(0.0, 0.0),
            AiryDerivative::Value,
            Scaling::Unscaled,
        )
        .unwrap();
        assert_eq!(nz, 0);
        assert!((ai.re - AI_C1).abs() < 1e-15);
        assert!(ai.im.abs() < 1e-15);
    }

    #[test]
    fn zairy_zero_deriv() {
        // Ai'(0) = -C2 ≈ -0.25881940379280680
        let (ai, nz, _status) = zairy(
            Complex64::new(0.0, 0.0),
            AiryDerivative::Derivative,
            Scaling::Unscaled,
        )
        .unwrap();
        assert_eq!(nz, 0);
        assert!((ai.re - (-AI_C2)).abs() < 1e-15);
        assert!(ai.im.abs() < 1e-15);
    }

    #[test]
    fn zairy_small_real() {
        // Ai(0.5) via power series — compare to known value
        // Ai(0.5) ≈ 0.23169360648083349
        let (ai, nz, _status) = zairy(
            Complex64::new(0.5, 0.0),
            AiryDerivative::Value,
            Scaling::Unscaled,
        )
        .unwrap();
        assert_eq!(nz, 0);
        assert!((ai.re - 0.23169360648083349).abs() < 1e-13);
        assert!(ai.im.abs() < 1e-15);
    }

    #[test]
    fn zairy_negative_real() {
        // Ai(-1.0): oscillatory region, via ZACAI
        // Ai(-1.0) ≈ 0.53556088329235176
        let (ai, nz, _status) = zairy(
            Complex64::new(-1.0, 0.0),
            AiryDerivative::Value,
            Scaling::Unscaled,
        )
        .unwrap();
        assert_eq!(nz, 0);
        assert!((ai.re - 0.53556088329235176).abs() < 1e-13);
        assert!(ai.im.abs() < 1e-15);
    }

    #[test]
    fn zairy_complex() {
        // Ai(1+i): complex argument (mpmath 30-digit reference)
        // Ai(1+i) ≈ 0.060458308371838149 - 0.15188956587718140i
        let (ai, nz, _status) = zairy(
            Complex64::new(1.0, 1.0),
            AiryDerivative::Value,
            Scaling::Unscaled,
        )
        .unwrap();
        assert_eq!(nz, 0);
        assert!((ai.re - 0.060458308371838149).abs() < 1e-13);
        assert!((ai.im - (-0.15188956587718140)).abs() < 1e-13);
    }

    // ── ZBIRY tests ──

    #[test]
    fn zbiry_zero() {
        // Bi(0) = BI_C1 ≈ 0.61492662744600074
        let (bi, _status) = zbiry(
            Complex64::new(0.0, 0.0),
            AiryDerivative::Value,
            Scaling::Unscaled,
        )
        .unwrap();
        assert!((bi.re - BI_C1).abs() < 1e-15);
        assert!(bi.im.abs() < 1e-15);
    }

    #[test]
    fn zbiry_zero_deriv() {
        // Bi'(0) = BI_C2 ≈ 0.44828835735382636
        let (bi, _status) = zbiry(
            Complex64::new(0.0, 0.0),
            AiryDerivative::Derivative,
            Scaling::Unscaled,
        )
        .unwrap();
        assert!((bi.re - BI_C2).abs() < 1e-15);
        assert!(bi.im.abs() < 1e-15);
    }

    #[test]
    fn zbiry_small_real() {
        // Bi(0.5) ≈ 0.85427704310315549 (mpmath 30-digit reference)
        let (bi, _status) = zbiry(
            Complex64::new(0.5, 0.0),
            AiryDerivative::Value,
            Scaling::Unscaled,
        )
        .unwrap();
        assert!((bi.re - 0.85427704310315549).abs() < 1e-13);
        assert!(bi.im.abs() < 1e-15);
    }

    #[test]
    fn airy_wronskian() {
        // Wronskian: Ai(z)*Bi'(z) - Ai'(z)*Bi(z) = 1/pi
        let z = Complex64::new(0.5, 0.3);
        let (ai, _, _) = zairy(z, AiryDerivative::Value, Scaling::Unscaled).unwrap();
        let (ai_p, _, _) = zairy(z, AiryDerivative::Derivative, Scaling::Unscaled).unwrap();
        let (bi, _) = zbiry(z, AiryDerivative::Value, Scaling::Unscaled).unwrap();
        let (bi_p, _) = zbiry(z, AiryDerivative::Derivative, Scaling::Unscaled).unwrap();

        // W = Ai*Bi' - Ai'*Bi
        let w = ai * bi_p - ai_p * bi;
        let inv_pi = 1.0 / core::f64::consts::PI;
        assert!(
            (w.re - inv_pi).abs() < 1e-13,
            "Wronskian real part: got {}, expected {}",
            w.re,
            inv_pi
        );
        assert!(w.im.abs() < 1e-13, "Wronskian imaginary part: got {}", w.im);
    }
}
