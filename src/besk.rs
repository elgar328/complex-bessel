//! K Bessel function upper interface.
//!
//! Translation of Fortran ZBESK from TOMS 644 / SLATEC.
//! Dispatches to `zbknu` (right half-plane), `zacon` (left half-plane),
//! and `zbunk` (large order uniform asymptotics).

use num_complex::Complex;
use num_traits::Float;

use crate::algo::acon::zacon;
use crate::algo::bknu::zbknu;
use crate::algo::bunk::zbunk;
use crate::algo::uoik::zuoik;
use crate::machine::BesselFloat;
use crate::types::{BesselError, BesselResult, Scaling};
use crate::utils::zabs;

/// Compute K_{ν+j}(z) for j = 0, 1, ..., n-1.
///
/// This is the main entry point for K Bessel function computation,
/// equivalent to Fortran ZBESK in TOMS 644.
///
/// # Parameters
/// - `fnu`: starting order ν ≥ 0
/// - `z`: complex argument (z ≠ 0)
/// - `n`: number of sequence members ≥ 1
/// - `scaling`: `Unscaled` = K_ν(z), `Exponential` = e^z · K_ν(z)
///
/// # Returns
/// `BesselResult` with computed values and underflow count.
///
/// # Errors
/// - `InvalidInput`: z = 0, ν < 0, or n < 1
/// - `Overflow`: |z| too small or result magnitude overflow
/// - `PrecisionLoss`: |z| or ν caused loss of more than half significant digits
/// - `TotalPrecisionLoss`: complete loss of significance
/// - `ConvergenceFailure`: algorithm did not converge
pub(crate) fn zbesk<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    scaling: Scaling,
    n: usize,
) -> Result<BesselResult<T>, BesselError> {
    let zero = T::zero();
    let czero = Complex::new(zero, zero);

    // ── Input validation (Fortran IERR=1) ──
    if n < 1 {
        return Err(BesselError::InvalidInput);
    }
    if fnu < zero {
        return Err(BesselError::InvalidInput);
    }
    if z == czero {
        return Err(BesselError::InvalidInput);
    }

    // ── Machine constants ──
    let tol = T::tol();
    let elim = T::elim();
    let alim = T::alim();
    let fnul = T::fnul();

    // AA = min(0.5/TOL, I1MACH(9)*0.5) — upper bound for range checks
    // I1MACH(9) = largest machine integer = 2147483647 for 32-bit
    let aa_tol = T::from(0.5).unwrap() / tol;
    let bb = T::from(2147483647.0 * 0.5).unwrap();
    let aa = aa_tol.min(bb);

    let az = zabs(z);
    let fn_val = fnu + T::from(n as f64 - 1.0).unwrap(); // FN = FNU + N - 1

    // ── Range check: total precision loss (IERR=4) ──
    if az > aa || fn_val > aa {
        return Err(BesselError::TotalPrecisionLoss);
    }

    // ── Range check: partial precision loss (IERR=3) ──
    let aa_sqrt = aa.sqrt();
    let mut precision_warning = false;
    if az > aa_sqrt || fn_val > aa_sqrt {
        precision_warning = true;
    }

    // ── Underflow limit: |z| too small ──
    let ufl = T::MACH_TINY * T::from(1.0e3).unwrap();
    if az < ufl {
        return Err(BesselError::Overflow);
    }

    // ── Dispatch based on FNU and z ──
    let mut nn = n;
    let mut nz: usize = 0;

    if fn_val > fnul {
        // ── Large order path: uniform asymptotics (ZBUNK) ──
        // Fortran ZBESK label 80-90 (zbsubs.f lines 1149-1159)
        let mr = if z.re >= zero {
            0i32
        } else if z.im < zero {
            -1i32
        } else {
            1i32
        };
        let (y, nw) = zbunk(z, fnu, scaling, mr, nn, tol, elim, alim);
        if nw < 0 {
            return if nw == -1 {
                Err(BesselError::Overflow)
            } else {
                Err(BesselError::ConvergenceFailure)
            };
        }
        nz += nw as usize;
        return if precision_warning {
            Err(BesselError::PrecisionLoss)
        } else {
            Ok(BesselResult {
                values: y,
                underflow_count: nz,
            })
        };
    }

    // ── Small-to-moderate order path ──
    // Fortran: check for overflow when FN > 1 and |z| is tiny
    let rl = T::rl();
    let fnul = T::fnul();

    if fn_val > T::one() {
        if fn_val > T::from(2.0).unwrap() {
            // FN > 2: ZUOIK overflow/underflow pre-check (Fortran lines 89-95)
            let (_uoik_y, nuf) = zuoik(z, fnu, scaling, 2, nn, tol, elim, alim);
            if nuf < 0 {
                return Err(BesselError::Overflow);
            }
            nz += nuf as usize;
            nn -= nuf as usize;
            if nn == 0 {
                return Ok(BesselResult {
                    values: vec![Complex::new(T::zero(), T::zero()); n],
                    underflow_count: nz,
                });
            }
        }

        // Check for overflow: -FN * ln(0.5 * |z|) > ELIM
        if az <= tol {
            let half = T::from(0.5).unwrap();
            let aln = -(fn_val * (half * az).ln());
            if aln > elim {
                return Err(BesselError::Overflow);
            }
        }
    }

    // ── Main computation dispatch ──
    let y = if z.re >= zero {
        // Right half-plane: direct ZBKNU
        let (y, nw) = zbknu(z, fnu, scaling, nn, tol, elim, alim)?;
        nz += nw;
        y
    } else {
        // Left half-plane: analytic continuation (ZACON)
        let mr = if z.im < T::zero() { -1i32 } else { 1i32 };
        let (y_acon, nw) = zacon(z, fnu, scaling, mr, nn, rl, fnul, tol, elim, alim)?;
        nz += nw;
        y_acon
    };

    if precision_warning {
        // Fortran IERR=3: computation completed but with precision loss.
        // We return Ok with the result; the caller should check the values.
        // TODO: consider adding a warning field to BesselResult
    }

    Ok(BesselResult {
        values: y,
        underflow_count: nz,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    // Reference values: Fortran TOMS 644 (zbsubs.f, revision 930101).
    // Full reference set: tests/reference_values/besk_f64.json
    // Comprehensive tests: tests/phase2_milestone.rs

    // ── Input validation tests ──

    #[test]
    fn besk_z_zero_returns_error() {
        let z = Complex64::new(0.0, 0.0);
        assert!(matches!(
            zbesk(z, 0.0, Scaling::Unscaled, 1),
            Err(BesselError::InvalidInput)
        ));
    }

    #[test]
    fn besk_negative_order_returns_error() {
        let z = Complex64::new(1.0, 0.0);
        assert!(matches!(
            zbesk(z, -1.0, Scaling::Unscaled, 1),
            Err(BesselError::InvalidInput)
        ));
    }

    #[test]
    fn besk_n_zero_returns_error() {
        let z = Complex64::new(1.0, 0.0);
        assert!(matches!(
            zbesk(z, 0.0, Scaling::Unscaled, 0),
            Err(BesselError::InvalidInput)
        ));
    }

    #[test]
    fn besk_left_half_plane() {
        // K_0(-1+i) via ZACON analytic continuation
        let z = Complex64::new(-1.0, 1.0);
        let result = zbesk(z, 0.0, Scaling::Unscaled, 1);
        assert!(result.is_ok(), "K_0(-1+i) should succeed with ZACON");
    }

    // ── Smoke tests (basic correctness) ──

    #[test]
    fn besk_k0_real() {
        // Fortran TOMS 644: K_0(1.0) = 0.42102443824070834
        let z = Complex64::new(1.0, 0.0);
        let result = zbesk(z, 0.0, Scaling::Unscaled, 1).unwrap();
        assert_eq!(result.underflow_count, 0);
        assert!((result.values[0].re - 0.42102443824070834).abs() < 1e-14);
        assert!(result.values[0].im.abs() < 1e-14);
    }

    #[test]
    fn besk_complex_arg() {
        // Fortran TOMS 644: K_0(1+i) = (0.08019772694651774, -0.35727745928533017)
        let z = Complex64::new(1.0, 1.0);
        let result = zbesk(z, 0.0, Scaling::Unscaled, 1).unwrap();
        let expected = Complex64::new(0.08019772694651774, -0.35727745928533017);
        let err = (result.values[0] - expected).norm() / expected.norm();
        assert!(err < 2e-14, "K_0(1+i) rel err = {err:.2e}");
    }

    #[test]
    fn besk_sequence() {
        // Fortran TOMS 644: K_{0,1,2}(2.0)
        let z = Complex64::new(2.0, 0.0);
        let result = zbesk(z, 0.0, Scaling::Unscaled, 3).unwrap();
        assert_eq!(result.values.len(), 3);
        assert_eq!(result.underflow_count, 0);
        assert!((result.values[0].re - 0.11389387274953341).abs() < 1e-14);
        assert!((result.values[1].re - 0.13986588181652246).abs() < 1e-14);
        assert!((result.values[2].re - 0.2537597545660559).abs() < 1e-13);
    }
}
