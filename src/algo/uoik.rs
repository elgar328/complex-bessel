//! Overflow/underflow test on I and K sequences using uniform asymptotic expansions.
//!
//! Translation of Fortran ZUOIK from TOMS 644 / SLATEC (zbsubs.f lines 3979-4173).
//!
//! IFORM=1 (ZUNIK path) for |Im(z)| <= |Re(z)|*sqrt(3).
//! IFORM=2 (ZUNHJ path) for |Im(z)| > |Re(z)|*sqrt(3).

#![allow(clippy::too_many_arguments)]
#![allow(clippy::excessive_precision)]

use num_complex::Complex;

use crate::algo::constants::AIC;
use crate::algo::uchk::zuchk;
use crate::algo::unhj::zunhj;
use crate::algo::unik::zunik;
use crate::machine::BesselFloat;
use crate::types::{IkFlag, Scaling, SumOption};
use crate::utils::zabs;

/// Overflow/underflow pre-check using uniform asymptotic leading terms.
///
/// Returns `nuf` where:
/// - `nuf = 0`: last member is on scale (no overflow/underflow detected)
/// - `nuf = -1`: overflow would occur
/// - `nuf > 0` (ikflg=1 only): last `nuf` Y values set to zero (underflow)
/// - `nuf = n` (ikflg=2): all Y values set to zero (underflow)
///
/// # Parameters
/// - `ikflg`: `IkFlag::I` for I function test, `IkFlag::K` for K function test
/// - `y`: output slice (length determines number of sequence members)
pub(crate) fn zuoik<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    ikflg: IkFlag,
    y: &mut [Complex<T>],
    tol: T,
    elim: T,
    alim: T,
) -> i32 {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);
    let n = y.len();
    y.fill(czero);
    let mut nuf: i32 = 0;
    let mut nn = n;

    // AIC = ln(sqrt(pi/2)) = 0.5*ln(pi/2) (Fortran line 4017)
    let aic = T::from_f64(AIC);

    // Reflect to right half plane (Fortran lines 4020-4026)
    let zr = if z.re >= zero { z } else { -z };

    // Determine form (Fortran lines 4028-4031)
    let ax = z.re.abs() * T::from_f64(1.7321);
    let ay = z.im.abs();
    let iform: i32 = if ay > ax { 2 } else { 1 };

    // Determine order for test (Fortran lines 4032-4036)
    let mut gnu = fnu.max(one);
    if ikflg == IkFlag::K {
        let fnn = T::from_f64(nn as f64);
        let gnn = fnu + fnn - one;
        gnu = gnn.max(fnn);
    }

    // ── Compute ZN for iform=2 (Fortran lines 4051-4055) ──
    let zn = if iform == 2 {
        // ZNR = ZRI, ZNI = -ZRR; if ZI <= 0 then ZNR = -ZNR
        let znr = if z.im <= zero { -zr.im } else { zr.im };
        Complex::new(znr, -zr.re)
    } else {
        czero // unused for iform=1
    };

    // ── Compute leading terms (Fortran lines 4043-4060) ──
    let (mut cz, phi, aarg) = if iform == 1 {
        // ZUNIK path (Fortran lines 4044-4048)
        let result = zunik(zr, gnu, ikflg, SumOption::SkipSum, tol, None);
        (result.zeta2 - result.zeta1, result.phi, zero)
    } else {
        // ZUNHJ path (Fortran lines 4056-4060)
        let result = zunhj(zn, gnu, SumOption::SkipSum, tol);
        (result.zeta2 - result.zeta1, result.phi, zabs(result.arg))
    };

    // ── Apply KODE and IKFLG adjustments (Fortran lines 4062-4071) ──
    if kode == Scaling::Exponential {
        cz = cz - zr;
    }
    if ikflg == IkFlag::K {
        cz = -cz;
    }
    let aphi = zabs(phi);
    let mut rcz = cz.re;

    // ── Overflow test (Fortran lines 4075-4080) ──
    if rcz > elim {
        return -1; // label 210
    }
    if rcz >= alim {
        // Between ALIM and ELIM: refine (Fortran lines 4077-4079)
        rcz = rcz + aphi.ln();
        if iform == 2 {
            rcz = rcz - T::from_f64(0.25) * aarg.ln() - aic;
        }
        if rcz > elim {
            return -1; // label 210
        }
        // Falls through to label 130 (on scale)
    } else {
        // RCZ < ALIM: check underflow (label 80, Fortran lines 4085-4096)
        if rcz < -elim {
            // Complete underflow (label 90)
            y[..nn].fill(czero);
            return nn as i32;
        }
        if rcz <= -alim {
            // Between -ELIM and -ALIM: refine (Fortran lines 4087-4089)
            rcz = rcz + aphi.ln();
            if iform == 2 {
                rcz = rcz - T::from_f64(0.25) * aarg.ln() - aic;
            }
            if rcz <= -elim {
                // Complete underflow (label 90)
                y[..nn].fill(czero);
                return nn as i32;
            }
            // Refined check with ZUCHK (label 110, Fortran lines 4098-4112)
            let ascle = T::from_f64(1.0e3) * T::MACH_TINY / tol;
            let log_phi = phi.ln();
            cz.im = cz.im + log_phi.im;
            if iform == 2 {
                // Fortran lines 4103-4105: subtract 0.25*ln(arg) + AIC
                let log_arg = Complex::new(aarg, zero).ln();
                cz.im = cz.im - T::from_f64(0.25) * log_arg.im;
            }
            // label 120 (Fortran lines 4107-4111)
            let cz_check = Complex::new(rcz, cz.im).exp() / tol;
            if zuchk(cz_check, ascle, tol) {
                // Underflow (label 90)
                y[..nn].fill(czero);
                return nn as i32;
            }
        }
        // else: -ALIM < RCZ < ALIM → on scale (label 130)
    }

    // ── Label 130: on scale (Fortran lines 4114-4115) ──
    if ikflg == IkFlag::K {
        return nuf;
    }
    if n == 1 {
        return nuf;
    }

    // ── I function per-member underflow chain (label 140, Fortran lines 4119-4168) ──
    loop {
        gnu = fnu + T::from_f64((nn - 1) as f64);

        let (phi2, aarg2);
        if iform == 1 {
            // ZUNIK call (Fortran lines 4122-4126)
            let result2 = zunik(zr, gnu, ikflg, SumOption::SkipSum, tol, None);
            cz = result2.zeta2 - result2.zeta1;
            phi2 = result2.phi;
            aarg2 = zero;
        } else {
            // ZUNHJ call (Fortran lines 4129-4133)
            let result2 = zunhj(zn, gnu, SumOption::SkipSum, tol);
            cz = result2.zeta2 - result2.zeta1;
            phi2 = result2.phi;
            aarg2 = zabs(result2.arg);
        }

        // Apply KODE adjustment (Fortran lines 4135-4137)
        if kode == Scaling::Exponential {
            cz = cz - zr;
        }

        let aphi2 = zabs(phi2);
        rcz = cz.re;

        // Underflow test (Fortran lines 4141-4145)
        if rcz < -elim {
            // label 180: underflow this member
            y[nn - 1] = czero;
            nn -= 1;
            nuf += 1;
            if nn == 0 {
                return nuf;
            }
            continue;
        }
        if rcz > -alim {
            // On scale, done
            return nuf;
        }

        // Between -ELIM and -ALIM: refine (Fortran lines 4143-4145)
        rcz = rcz + aphi2.ln();
        if iform == 2 {
            rcz = rcz - T::from_f64(0.25) * aarg2.ln() - aic;
        }
        if rcz <= -elim {
            // label 180: underflow this member
            y[nn - 1] = czero;
            nn -= 1;
            nuf += 1;
            if nn == 0 {
                return nuf;
            }
            continue;
        }

        // Refined check with ZUCHK (label 190, Fortran lines 4154-4168)
        let ascle = T::from_f64(1.0e3) * T::MACH_TINY / tol;
        let log_phi2 = phi2.ln();
        cz.im = cz.im + log_phi2.im;
        if iform == 2 {
            // Fortran lines 4159-4161
            let log_arg2 = Complex::new(aarg2, zero).ln();
            cz.im = cz.im - T::from_f64(0.25) * log_arg2.im;
        }
        // label 200 (Fortran lines 4163-4167)
        let cz_check = Complex::new(rcz, cz.im).exp() / tol;
        if zuchk(cz_check, ascle, tol) {
            // label 180: underflow this member
            y[nn - 1] = czero;
            nn -= 1;
            nuf += 1;
            if nn == 0 {
                return nuf;
            }
            continue;
        }

        // On scale, done
        return nuf;
    }
}
