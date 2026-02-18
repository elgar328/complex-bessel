//! K Bessel function upper interface.
//!
//! Translation of Fortran ZBESK from TOMS 644 / SLATEC.
//! Dispatches to `zbknu` (right half-plane), `zacon` (left half-plane),
//! and `zbunk` (large order uniform asymptotics).

use num_complex::Complex;

use crate::algo::acon::zacon;
use crate::algo::bknu::zbknu;
use crate::algo::bunk::zbunk;
use crate::algo::uoik::zuoik;
use crate::machine::BesselFloat;
use crate::types::{BesselError, BesselStatus, Scaling};
use crate::utils::zabs;

/// Compute K_{ν+j}(z) for j = 0, 1, ..., n-1 into the provided slice.
///
/// This is the main entry point for K Bessel function computation,
/// equivalent to Fortran ZBESK in TOMS 644.
///
/// # Parameters
/// - `fnu`: starting order ν ≥ 0
/// - `z`: complex argument (z ≠ 0)
/// - `y`: output slice of length n ≥ 1; results are written here
/// - `scaling`: `Unscaled` = K_ν(z), `Exponential` = e^z · K_ν(z)
///
/// # Returns
/// `(nz, status)` where `nz` is the underflow count and `status`
/// indicates precision quality.
///
/// # Errors
/// - `InvalidInput`: z = 0, ν < 0, or n < 1
/// - `Overflow`: |z| too small or result magnitude overflow
/// - `TotalPrecisionLoss`: complete loss of significance
/// - `ConvergenceFailure`: algorithm did not converge
pub(crate) fn zbesk<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    scaling: Scaling,
    y: &mut [Complex<T>],
) -> Result<(usize, BesselStatus), BesselError> {
    let n = y.len();
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
        let nw = zbunk(z, fnu, scaling, mr, &mut y[..nn], tol, elim, alim);
        if nw < 0 {
            return if nw == -1 {
                Err(BesselError::Overflow)
            } else {
                Err(BesselError::ConvergenceFailure)
            };
        }
        nz += nw as usize;
        let status = if precision_warning {
            BesselStatus::ReducedPrecision
        } else {
            BesselStatus::Normal
        };
        return Ok((nz, status));
    }

    // ── Small-to-moderate order path ──
    // Fortran: check for overflow when FN > 1 and |z| is tiny
    let rl = T::rl();

    if fn_val > T::one() {
        if fn_val > T::from(2.0).unwrap() {
            // FN > 2: ZUOIK overflow/underflow pre-check (Fortran lines 89-95)
            // zuoik zeros the buffer; y[..nn] will be overwritten by actual computation later.
            let nuf = zuoik(z, fnu, scaling, 2, &mut y[..nn], tol, elim, alim);
            if nuf < 0 {
                return Err(BesselError::Overflow);
            }
            nz += nuf as usize;
            nn -= nuf as usize;
            if nn == 0 {
                // All members underflowed; y is already zeroed by zuoik.
                let status = if precision_warning {
                    BesselStatus::ReducedPrecision
                } else {
                    BesselStatus::Normal
                };
                return Ok((nz, status));
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
    if z.re >= zero {
        // Right half-plane: direct ZBKNU
        let nw = zbknu(z, fnu, scaling, &mut y[..nn], tol, elim, alim)?;
        nz += nw;
    } else {
        // Left half-plane: analytic continuation (ZACON)
        let mr = if z.im < T::zero() { -1i32 } else { 1i32 };
        let nw = zacon(z, fnu, scaling, mr, &mut y[..nn], rl, fnul, tol, elim, alim)?;
        nz += nw;
    }

    let status = if precision_warning {
        BesselStatus::ReducedPrecision
    } else {
        BesselStatus::Normal
    };

    Ok((nz, status))
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn besk_z_zero_returns_error() {
        let z = Complex64::new(0.0, 0.0);
        let mut buf = [Complex64::new(0.0, 0.0)];
        assert!(matches!(
            zbesk(z, 0.0, Scaling::Unscaled, &mut buf),
            Err(BesselError::InvalidInput)
        ));
    }

    #[test]
    fn besk_negative_order_returns_error() {
        let z = Complex64::new(1.0, 0.0);
        let mut buf = [Complex64::new(0.0, 0.0)];
        assert!(matches!(
            zbesk(z, -1.0, Scaling::Unscaled, &mut buf),
            Err(BesselError::InvalidInput)
        ));
    }

    #[test]
    fn besk_n_zero_returns_error() {
        let z = Complex64::new(1.0, 0.0);
        let mut buf: [Complex64; 0] = [];
        assert!(matches!(
            zbesk(z, 0.0, Scaling::Unscaled, &mut buf),
            Err(BesselError::InvalidInput)
        ));
    }
}
