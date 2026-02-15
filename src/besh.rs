//! Hankel function upper interface.
//!
//! Translation of Fortran ZBESH from TOMS 644 / SLATEC.
//! Dispatches to `zbknu` (right half-plane of rotated argument),
//! `zacon` (left half-plane), and `zbunk` (large order uniform asymptotics).
//!
//! The Hankel function is computed via:
//!   H(m, fnu, z) = (1/mp)*exp(-mp*fnu)*K(fnu, z*exp(-mp))
//! where mp = (3-2*m)*i*pi/2, m = 1 or 2.

// Exact Fortran constants — preserve verbatim.
#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]

use num_complex::Complex;
use num_traits::Float;

use crate::algo::bknu::zbknu;
use crate::machine::BesselFloat;
use crate::types::{BesselError, BesselResult, HankelKind, Scaling};
use crate::utils::zabs;

/// pi/2, Fortran DATA constant (zbsubs.f line 173)
const HPI: f64 = 1.57079632679489662;

/// Compute H_{fnu+j}^(m)(z) for j = 0, 1, ..., n-1.
///
/// This is the main entry point for Hankel function computation,
/// equivalent to Fortran ZBESH in TOMS 644.
///
/// # Parameters
/// - `z`: complex argument (z != 0)
/// - `fnu`: starting order (>= 0)
/// - `kind`: `First` (m=1) or `Second` (m=2)
/// - `scaling`: `Unscaled` = H(m,v,z), `Exponential` = H(m,v,z)*exp(-(3-2m)*i*z)
/// - `n`: number of sequence members (>= 1)
///
/// # Current Limitations
/// - ZACON (Phase 4) not implemented: H^(1) with Im(z) < 0 returns error
/// - ZBUNK (Phase 5b) not implemented: large orders with Re(zn) < 0 returns error
/// - ZUOIK (Phase 5a) not implemented: underflow pre-check skipped for fn > 2
pub(crate) fn zbesh<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kind: HankelKind,
    scaling: Scaling,
    n: usize,
) -> Result<BesselResult<T>, BesselError> {
    let zero = T::zero();
    let one = T::one();
    let two = T::from(2.0).unwrap();

    // ── Input validation (Fortran IERR=1, lines 176-183) ──
    if n < 1 {
        return Err(BesselError::InvalidInput);
    }
    if fnu < zero {
        return Err(BesselError::InvalidInput);
    }
    if z == Complex::new(zero, zero) {
        return Err(BesselError::InvalidInput);
    }

    let nn = n;

    // ── Machine constants (Fortran lines 196-208) ──
    let tol = T::tol();
    let elim = T::elim();
    let alim = T::alim();
    let fnul = T::fnul();

    let fn_val = fnu + T::from((nn - 1) as f64).unwrap();

    // M = 1 or 2, MM = 3 - 2*M, FMM = float(MM)
    let m: i32 = match kind {
        HankelKind::First => 1,
        HankelKind::Second => 2,
    };
    let mm = 3 - 2 * m; // +1 for m=1, -1 for m=2
    let fmm = T::from(mm as f64).unwrap();

    // Rotated argument: ZN = z * exp(-mp) where mp = mm*HPI*i
    // ZNR = FMM * ZI, ZNI = -FMM * ZR (Fortran lines 212-213)
    let znr = fmm * z.im;
    let zni = -fmm * z.re;
    let zn = Complex::new(znr, zni);

    // ── Range check (Fortran lines 217-225) ──
    let az = zabs(z);
    let aa_tol = T::from(0.5).unwrap() / tol;
    let bb = T::from(2147483647.0 * 0.5).unwrap();
    let aa = aa_tol.min(bb);

    if az > aa || fn_val > aa {
        return Err(BesselError::TotalPrecisionLoss);
    }

    let aa_sqrt = aa.sqrt();
    let mut precision_warning = false;
    if az > aa_sqrt || fn_val > aa_sqrt {
        precision_warning = true;
    }

    // ── Underflow limit (Fortran line 229-230) ──
    let ufl = T::MACH_TINY * T::from(1.0e3).unwrap();
    if az < ufl {
        return Err(BesselError::Overflow);
    }

    // ── Overflow pre-checks and dispatch (Fortran lines 231-268) ──
    let mut nz: usize = 0;
    let mut nn_eff = nn;

    let cy = if fnu > fnul {
        // Large order: ZBUNK (not implemented)
        // Fallback to zbknu if Re(zn) >= 0 (Fortran label 90 path)
        if znr >= zero {
            let (y, nw) = zbknu(zn, fnu, scaling, nn_eff, tol, elim, alim)?;
            nz += nw;
            y
        } else {
            return Err(BesselError::ConvergenceFailure);
        }
    } else {
        // ── Small-to-moderate order overflow checks (Fortran lines 232-249) ──
        if fn_val > one {
            if fn_val > two {
                // fn > 2: would call ZUOIK for overflow/underflow pre-check
                // TODO: Phase 5a — ZUOIK not yet implemented. Skip pre-check.
            }

            // Check for overflow when az is very small (Fortran lines 234-237)
            if fn_val <= two && az <= tol {
                let half = T::from(0.5).unwrap();
                let arg = half * az;
                let aln = -fn_val * arg.ln();
                if aln > elim {
                    return Err(BesselError::Overflow);
                }
            }
        }

        // ── Main dispatch (Fortran label 70, lines 251-268) ──
        if znr < zero || (znr == zero && zni < zero && m == 2) {
            // Left half-plane of rotated argument: ZACON needed
            // TODO: Phase 4 — ZACON not yet implemented
            return Err(BesselError::ConvergenceFailure);
        }

        // Right half-plane: direct zbknu on rotated argument
        let (y, nw) = zbknu(zn, fnu, scaling, nn_eff, tol, elim, alim)?;
        nz += nw;
        y
    };

    // ── Phase post-processing (Fortran label 110, lines 285-336) ──
    // H(m, fnu, z) = -FMM*(i/HPI)*exp(-FMM*fnu*HPI*i)*K(fnu, zn)
    let hpi_t = T::from(HPI).unwrap();
    let sgn = hpi_t.copysign(-fmm); // sign(HPI, -FMM)

    // Decompose FNU for precision: inu, inuh, ir (Fortran lines 296-299)
    let inu = fnu.to_i32().unwrap();
    let inuh = inu / 2;
    let ir = inu - 2 * inuh; // 0 or 1
    let arg = (fnu - T::from((inu - ir) as f64).unwrap()) * sgn;
    let rhpi = one / sgn;

    // CSGN = RHPI * (-sin(arg) + i*cos(arg)) (Fortran lines 303-304)
    let mut csgnr = -rhpi * arg.sin();
    let mut csgni = rhpi * arg.cos();

    // If inuh is odd, negate CSGN (Fortran lines 305-309)
    if inuh % 2 != 0 {
        csgnr = -csgnr;
        csgni = -csgni;
    }

    // ZTI = -FMM (Fortran line 311)
    let zti = -fmm;
    let rtol = one / tol;
    let ascle = ufl * rtol;

    // Apply phase factor to each component (Fortran lines 314-336)
    let mut result = cy;
    for item in result.iter_mut().take(nn_eff) {
        let mut aa_val = item.re;
        let mut bb_val = item.im;
        let mut atol = one;

        // Underflow protection: scale up tiny values (Fortran lines 323-327)
        if aa_val.abs().max(bb_val.abs()) <= ascle {
            aa_val = aa_val * rtol;
            bb_val = bb_val * rtol;
            atol = tol;
        }

        // CY(I) = (AA + i*BB) * CSGN * ATOL (Fortran lines 329-332)
        let str = aa_val * csgnr - bb_val * csgni;
        let sti = aa_val * csgni + bb_val * csgnr;
        *item = Complex::new(str * atol, sti * atol);

        // Advance CSGN: CSGN *= (0, ZTI) (Fortran lines 333-335)
        let str_new = -csgni * zti;
        csgni = csgnr * zti;
        csgnr = str_new;
    }

    if precision_warning {
        // Fortran IERR=3: computation completed but with precision loss
    }

    Ok(BesselResult {
        values: result,
        underflow_count: nz,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    // ── Input validation ──

    #[test]
    fn besh_z_zero_returns_error() {
        let z = Complex64::new(0.0, 0.0);
        assert!(matches!(
            zbesh(z, 0.0, HankelKind::First, Scaling::Unscaled, 1),
            Err(BesselError::InvalidInput)
        ));
    }

    #[test]
    fn besh_negative_order_returns_error() {
        let z = Complex64::new(1.0, 0.0);
        assert!(matches!(
            zbesh(z, -1.0, HankelKind::First, Scaling::Unscaled, 1),
            Err(BesselError::InvalidInput)
        ));
    }

    #[test]
    fn besh_n_zero_returns_error() {
        let z = Complex64::new(1.0, 0.0);
        assert!(matches!(
            zbesh(z, 0.0, HankelKind::First, Scaling::Unscaled, 0),
            Err(BesselError::InvalidInput)
        ));
    }

    // ── Basic correctness ──

    #[test]
    fn besh_h1_conjugate_h2_real_axis() {
        // For real x > 0: H^(1)(v, x) = conj(H^(2)(v, x))
        let z = Complex64::new(2.0, 0.0);
        let fnu = 0.5;
        let h1 = zbesh(z, fnu, HankelKind::First, Scaling::Unscaled, 1).unwrap();
        let h2 = zbesh(z, fnu, HankelKind::Second, Scaling::Unscaled, 1).unwrap();
        let diff = (h1.values[0] - h2.values[0].conj()).norm();
        let scale = h1.values[0].norm();
        assert!(
            diff / scale < 1e-14,
            "H^(1) != conj(H^(2)) on real axis, err = {:.2e}",
            diff / scale
        );
    }

    #[test]
    fn besh_h1_upper_half_plane() {
        // H^(1)(0, 1+i) should be computable (Im(z) >= 0)
        let z = Complex64::new(1.0, 1.0);
        let result = zbesh(z, 0.0, HankelKind::First, Scaling::Unscaled, 1);
        assert!(result.is_ok(), "H^(1)(0, 1+i) should succeed");
    }

    #[test]
    fn besh_h2_lower_half_plane() {
        // H^(2)(0, 1-i) should be computable (Im(z) <= 0)
        let z = Complex64::new(1.0, -1.0);
        let result = zbesh(z, 0.0, HankelKind::Second, Scaling::Unscaled, 1);
        assert!(result.is_ok(), "H^(2)(0, 1-i) should succeed");
    }

    #[test]
    fn besh_h1_lower_half_returns_error() {
        // H^(1)(0, 1-i) needs ZACON (Im(z) < 0, m=1)
        let z = Complex64::new(1.0, -1.0);
        let result = zbesh(z, 0.0, HankelKind::First, Scaling::Unscaled, 1);
        assert!(result.is_err(), "H^(1) with Im(z)<0 should fail (no ZACON)");
    }
}
