//! Overflow/underflow test on I and K sequences using uniform asymptotic expansions.
//!
//! Translation of Fortran ZUOIK from TOMS 644 / SLATEC (zbsubs.f lines 3979-4173).
//!
//! IFORM=1 (ZUNIK path) for |Im(z)| <= |Re(z)|*sqrt(3).
//! IFORM=2 (ZUNHJ path) for |Im(z)| > |Re(z)|*sqrt(3).

#![allow(clippy::too_many_arguments)]
#![allow(clippy::excessive_precision)]
#![allow(unused_assignments)]

use num_complex::Complex;

use crate::algo::constants::AIC;
use crate::algo::uchk::zuchk;
use crate::algo::unhj::zunhj;
use crate::algo::unik::zunik;
use crate::machine::BesselFloat;
use crate::types::Scaling;
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
/// - `ikflg`: 1 for I function test, 2 for K function test
/// - `y`: output slice (length determines number of sequence members)
pub(crate) fn zuoik<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    ikflg: i32,
    y: &mut [Complex<T>],
    tol: T,
    elim: T,
    alim: T,
) -> i32 {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);
    let n = y.len();
    for v in y.iter_mut() {
        *v = czero;
    }
    let mut nuf: i32 = 0;
    let mut nn = n;

    // AIC = ln(sqrt(pi/2)) = 0.5*ln(pi/2) (Fortran line 4017)
    let aic = T::from(AIC).unwrap();

    // Reflect to right half plane (Fortran lines 4020-4026)
    let (zrr, zri) = if z.re >= zero {
        (z.re, z.im)
    } else {
        (-z.re, -z.im)
    };
    let zbr = zrr;
    let zbi = zri;

    // Determine form (Fortran lines 4028-4031)
    let ax = z.re.abs() * T::from(1.7321).unwrap();
    let ay = z.im.abs();
    let iform: i32 = if ay > ax { 2 } else { 1 };

    // Determine order for test (Fortran lines 4032-4036)
    let mut gnu = fnu.max(one);
    if ikflg == 2 {
        let fnn = T::from(nn as f64).unwrap();
        let gnn = fnu + fnn - one;
        gnu = gnn.max(fnn);
    }

    // ── Compute ZN for iform=2 (Fortran lines 4051-4055) ──
    let (znr, zni) = if iform == 2 {
        // ZNR = ZRI, ZNI = -ZRR; if ZI <= 0 then ZNR = -ZNR
        let mut znr_val = zri;
        let zni_val = -zrr;
        if z.im <= zero {
            znr_val = -znr_val;
        }
        (znr_val, zni_val)
    } else {
        (zero, zero) // unused for iform=1
    };

    // ── Compute leading terms (Fortran lines 4043-4060) ──
    let zr_arg = Complex::new(zrr, zri);
    let (mut czr, mut czi, phi, aarg) = if iform == 1 {
        // ZUNIK path (Fortran lines 4044-4048)
        let result = zunik(zr_arg, gnu, ikflg, 1, tol, None);
        let czr_val = -result.zeta1.re + result.zeta2.re;
        let czi_val = -result.zeta1.im + result.zeta2.im;
        (czr_val, czi_val, result.phi, zero)
    } else {
        // ZUNHJ path (Fortran lines 4056-4060)
        let zn = Complex::new(znr, zni);
        let result = zunhj(zn, gnu, 1, tol);
        let czr_val = -result.zeta1.re + result.zeta2.re;
        let czi_val = -result.zeta1.im + result.zeta2.im;
        let aarg_val = zabs(result.arg);
        (czr_val, czi_val, result.phi, aarg_val)
    };

    // ── Apply KODE and IKFLG adjustments (Fortran lines 4062-4071) ──
    if kode == Scaling::Exponential {
        czr = czr - zbr;
        czi = czi - zbi;
    }
    if ikflg == 2 {
        czr = -czr;
        czi = -czi;
    }
    let aphi = zabs(phi);
    let mut rcz = czr;

    // ── Overflow test (Fortran lines 4075-4080) ──
    if rcz > elim {
        return -1; // label 210
    }
    if rcz >= alim {
        // Between ALIM and ELIM: refine (Fortran lines 4077-4079)
        rcz = rcz + aphi.ln();
        if iform == 2 {
            rcz = rcz - T::from(0.25).unwrap() * aarg.ln() - aic;
        }
        if rcz > elim {
            return -1; // label 210
        }
        // Falls through to label 130 (on scale)
    } else {
        // RCZ < ALIM: check underflow (label 80, Fortran lines 4085-4096)
        if rcz < -elim {
            // Complete underflow (label 90)
            for item in y.iter_mut().take(nn) {
                *item = czero;
            }
            return nn as i32;
        }
        if rcz <= -alim {
            // Between -ELIM and -ALIM: refine (Fortran lines 4087-4089)
            rcz = rcz + aphi.ln();
            if iform == 2 {
                rcz = rcz - T::from(0.25).unwrap() * aarg.ln() - aic;
            }
            if rcz <= -elim {
                // Complete underflow (label 90)
                for item in y.iter_mut().take(nn) {
                    *item = czero;
                }
                return nn as i32;
            }
            // Refined check with ZUCHK (label 110, Fortran lines 4098-4112)
            let ascle = T::from(1.0e3).unwrap() * T::MACH_TINY / tol;
            let log_phi = phi.ln();
            czr = czr + log_phi.re;
            czi = czi + log_phi.im;
            if iform == 2 {
                // Fortran lines 4103-4105: subtract 0.25*ln(arg) + AIC
                let log_arg = Complex::new(aarg, zero).ln();
                czr = czr - T::from(0.25).unwrap() * log_arg.re - aic;
                czi = czi - T::from(0.25).unwrap() * log_arg.im;
            }
            // label 120 (Fortran lines 4107-4111)
            let ax_val = rcz.exp() / tol;
            let cz_check = Complex::new(ax_val * czi.cos(), ax_val * czi.sin());
            if zuchk(cz_check, ascle, tol) {
                // Underflow (label 90)
                for item in y.iter_mut().take(nn) {
                    *item = czero;
                }
                return nn as i32;
            }
        }
        // else: -ALIM < RCZ < ALIM → on scale (label 130)
    }

    // ── Label 130: on scale (Fortran lines 4114-4115) ──
    if ikflg == 2 {
        return nuf;
    }
    if n == 1 {
        return nuf;
    }

    // ── I function per-member underflow chain (label 140, Fortran lines 4119-4168) ──
    loop {
        gnu = fnu + T::from((nn - 1) as f64).unwrap();

        let (phi2, aarg2);
        if iform == 1 {
            // ZUNIK call (Fortran lines 4122-4126)
            let result2 = zunik(zr_arg, gnu, ikflg, 1, tol, None);
            czr = -result2.zeta1.re + result2.zeta2.re;
            czi = -result2.zeta1.im + result2.zeta2.im;
            phi2 = result2.phi;
            aarg2 = zero;
        } else {
            // ZUNHJ call (Fortran lines 4129-4133)
            let zn = Complex::new(znr, zni);
            let result2 = zunhj(zn, gnu, 1, tol);
            czr = -result2.zeta1.re + result2.zeta2.re;
            czi = -result2.zeta1.im + result2.zeta2.im;
            phi2 = result2.phi;
            aarg2 = zabs(result2.arg);
        }

        // Apply KODE adjustment (Fortran lines 4135-4137)
        if kode == Scaling::Exponential {
            czr = czr - zbr;
            czi = czi - zbi;
        }

        let aphi2 = zabs(phi2);
        rcz = czr;

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
            rcz = rcz - T::from(0.25).unwrap() * aarg2.ln() - aic;
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
        let ascle = T::from(1.0e3).unwrap() * T::MACH_TINY / tol;
        let log_phi2 = phi2.ln();
        czr = czr + log_phi2.re;
        czi = czi + log_phi2.im;
        if iform == 2 {
            // Fortran lines 4159-4161
            let log_arg2 = Complex::new(aarg2, zero).ln();
            czr = czr - T::from(0.25).unwrap() * log_arg2.re - aic;
            czi = czi - T::from(0.25).unwrap() * log_arg2.im;
        }
        // label 200 (Fortran lines 4163-4167)
        let ax_val = rcz.exp() / tol;
        let cz_check = Complex::new(ax_val * czi.cos(), ax_val * czi.sin());
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
