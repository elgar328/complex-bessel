//! Overflow/underflow test on I and K sequences using uniform asymptotic expansions.
//!
//! Translation of Fortran ZUOIK from TOMS 644 / SLATEC (zbsubs.f lines 3979-4173).
//!
//! IFORM=1 (ZUNIK path) is fully implemented.
//! IFORM=2 (ZUNHJ path) returns nuf=0 (stub, requires Phase 5b).

#![allow(clippy::too_many_arguments)]
#![allow(clippy::excessive_precision)]

use num_complex::Complex;
use num_traits::Float;

use crate::algo::uchk::zuchk;
use crate::algo::unik::zunik;
use crate::machine::BesselFloat;
use crate::types::Scaling;
use crate::utils::zabs;

/// Overflow/underflow pre-check using uniform asymptotic leading terms.
///
/// Returns `(y, nuf)` where:
/// - `nuf = 0`: last member is on scale (no overflow/underflow detected)
/// - `nuf = -1`: overflow would occur
/// - `nuf > 0` (ikflg=1 only): last `nuf` Y values set to zero (underflow)
/// - `nuf = n` (ikflg=2): all Y values set to zero (underflow)
///
/// # Parameters
/// - `ikflg`: 1 for I function test, 2 for K function test
/// - `n`: number of sequence members
pub(crate) fn zuoik<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    ikflg: i32,
    n: usize,
    tol: T,
    elim: T,
    alim: T,
) -> (Vec<Complex<T>>, i32) {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);
    let mut y = vec![czero; n];
    let mut nuf: i32 = 0;
    let mut nn = n;

    // AIC = ln(sqrt(pi/2)) = 0.5*ln(pi/2) (Fortran line 4017)
    let aic = T::from(1.265512123484645396).unwrap();

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

    // ── Compute leading terms (Fortran lines 4043-4060) ──
    if iform == 2 {
        // ZUNHJ path — not yet implemented (Phase 5b)
        return (y, 0);
    }

    // IFORM=1: ZUNIK path (Fortran lines 4044-4048)
    let zr_arg = Complex::new(zrr, zri);
    let result = zunik(zr_arg, gnu, ikflg, 1, tol, None);
    let mut czr = -result.zeta1.re + result.zeta2.re;
    let mut czi = -result.zeta1.im + result.zeta2.im;
    let phi = result.phi;

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
        return (y, -1); // label 210
    }
    if rcz >= alim {
        // Between ALIM and ELIM: refine (Fortran lines 4077-4079)
        rcz = rcz + aphi.ln();
        // iform==2 adjustment skipped (we checked iform==1 above)
        if rcz > elim {
            return (y, -1); // label 210
        }
        // Falls through to label 130 (on scale)
    } else {
        // RCZ < ALIM: check underflow (label 80, Fortran lines 4085-4096)
        if rcz < -elim {
            // Complete underflow (label 90)
            for item in y.iter_mut().take(nn) {
                *item = czero;
            }
            return (y, nn as i32);
        }
        if rcz <= -alim {
            // Between -ELIM and -ALIM: refine (Fortran lines 4087-4089)
            rcz = rcz + aphi.ln();
            if rcz <= -elim {
                // Complete underflow (label 90)
                for item in y.iter_mut().take(nn) {
                    *item = czero;
                }
                return (y, nn as i32);
            }
            // Refined check with ZUCHK (label 110, Fortran lines 4098-4112)
            let ascle = T::from(1.0e3).unwrap() * T::MACH_TINY / tol;
            let log_phi = phi.ln();
            // CZR += Re(log(phi)), CZI += Im(log(phi))
            // (RCZ was already updated with ln(APHI) = Re(log(phi)))
            let czi_full = czi + log_phi.im;
            // iform==2 adjustment skipped
            let ax_val = rcz.exp() / tol;
            let cz_check = Complex::new(ax_val * czi_full.cos(), ax_val * czi_full.sin());
            if zuchk(cz_check, ascle, tol) {
                // Underflow (label 90)
                for item in y.iter_mut().take(nn) {
                    *item = czero;
                }
                return (y, nn as i32);
            }
        }
        // else: -ALIM < RCZ < ALIM → on scale (label 130)
    }

    // ── Label 130: on scale (Fortran lines 4114-4115) ──
    if ikflg == 2 {
        return (y, nuf);
    }
    if n == 1 {
        return (y, nuf);
    }

    // ── I function per-member underflow chain (label 140, Fortran lines 4119-4168) ──
    loop {
        gnu = fnu + T::from((nn - 1) as f64).unwrap();

        // IFORM=1: ZUNIK call (Fortran lines 4122-4126)
        let result2 = zunik(zr_arg, gnu, ikflg, 1, tol, None);
        czr = -result2.zeta1.re + result2.zeta2.re;
        czi = -result2.zeta1.im + result2.zeta2.im;

        // Apply KODE adjustment (Fortran lines 4135-4137)
        if kode == Scaling::Exponential {
            czr = czr - zbr;
            czi = czi - zbi;
        }

        let aphi2 = zabs(result2.phi);
        rcz = czr;

        // Underflow test (Fortran lines 4141-4145)
        if rcz < -elim {
            // label 180: underflow this member
            y[nn - 1] = czero;
            nn -= 1;
            nuf += 1;
            if nn == 0 {
                return (y, nuf);
            }
            continue;
        }
        if rcz > -alim {
            // On scale, done
            return (y, nuf);
        }

        // Between -ELIM and -ALIM: refine (Fortran lines 4143-4145)
        rcz = rcz + aphi2.ln();
        if rcz <= -elim {
            // label 180: underflow this member
            y[nn - 1] = czero;
            nn -= 1;
            nuf += 1;
            if nn == 0 {
                return (y, nuf);
            }
            continue;
        }

        // Refined check with ZUCHK (label 190, Fortran lines 4154-4168)
        let ascle = T::from(1.0e3).unwrap() * T::MACH_TINY / tol;
        let log_phi2 = result2.phi.ln();
        let czi_full = czi + log_phi2.im;
        let ax_val = rcz.exp() / tol;
        let cz_check = Complex::new(ax_val * czi_full.cos(), ax_val * czi_full.sin());
        if zuchk(cz_check, ascle, tol) {
            // label 180: underflow this member
            y[nn - 1] = czero;
            nn -= 1;
            nuf += 1;
            if nn == 0 {
                return (y, nuf);
            }
            continue;
        }

        // On scale, done
        return (y, nuf);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    const TOL: f64 = 2.220446049250313e-16;
    const ELIM: f64 = 700.9217936944459;
    const ALIM: f64 = 664.8716455337102;

    #[test]
    fn zuoik_no_overflow_underflow() {
        // Normal case: moderate z, moderate fnu → on scale
        let z = Complex64::new(5.0, 1.0);
        let (_, nuf) = zuoik(z, 10.0, Scaling::Unscaled, 1, 3, TOL, ELIM, ALIM);
        assert_eq!(nuf, 0);
    }

    #[test]
    fn zuoik_k_function_no_overflow() {
        let z = Complex64::new(5.0, 1.0);
        let (_, nuf) = zuoik(z, 10.0, Scaling::Unscaled, 2, 3, TOL, ELIM, ALIM);
        assert_eq!(nuf, 0);
    }

    #[test]
    fn zuoik_iform2_stub() {
        // |Im(z)| > |Re(z)|*1.7321 → iform=2 → stub returns 0
        let z = Complex64::new(1.0, 10.0);
        let (_, nuf) = zuoik(z, 10.0, Scaling::Unscaled, 1, 3, TOL, ELIM, ALIM);
        assert_eq!(nuf, 0);
    }

    #[test]
    fn zuoik_regression_existing_tests() {
        // Ensure the function still works for cases used by Phase 2-4
        let z = Complex64::new(1.0, 0.5);
        let (_, nuf) = zuoik(z, 2.5, Scaling::Unscaled, 1, 2, TOL, ELIM, ALIM);
        // Should be on scale (no underflow/overflow)
        assert!(nuf >= 0);
    }
}
