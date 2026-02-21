//! Hankel function upper interface.
//!
//! Translation of Fortran ZBESH from TOMS 644 (zbsubs.f lines 5-354).
//! Dispatches to `zbknu` (right half-plane of rotated argument),
//! `zacon` (left half-plane), and `zbunk` (large order uniform asymptotics).
//!
//! The Hankel function is computed via:
//!   H(m, fnu, z) = (1/mp)*exp(-mp*fnu)*K(fnu, z*exp(-mp))
//! where mp = (3-2*m)*i*pi/2, m = 1 or 2.

use num_complex::Complex;

use crate::algo::acon::zacon;
use crate::algo::bknu::zbknu;
use crate::algo::bunk::zbunk;
use crate::algo::constants::HPI;
use crate::algo::uoik::zuoik;
use crate::machine::BesselFloat;
use crate::types::{Accuracy, Error, HankelKind, IkFlag, Scaling};
use crate::utils::{mul_i, mul_neg_i, zabs};

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
/// - `y`: output slice; length determines n (must be >= 1)
///
/// # Returns
/// `(nz, status)` where `nz` is the underflow count.
#[inline]
pub(crate) fn zbesh<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kind: HankelKind,
    scaling: Scaling,
    y: &mut [Complex<T>],
) -> Result<(usize, Accuracy), Error> {
    let zero = T::zero();
    let one = T::one();
    let two = T::from_f64(2.0);
    let czero = Complex::new(zero, zero);

    let n = y.len();

    // ── Input validation (Fortran IERR=1, lines 176-183) ──
    if n < 1 {
        return Err(Error::InvalidInput);
    }
    if fnu < zero {
        return Err(Error::InvalidInput);
    }
    if z == czero {
        return Err(Error::InvalidInput);
    }

    let nn = n;

    // ── Machine constants (Fortran lines 196-208) ──
    let tol = T::tol();
    let elim = T::elim();
    let alim = T::alim();
    let fnul = T::fnul();

    let fn_val = fnu + T::from_f64((nn - 1) as f64);

    // M = 1 or 2, MM = 3 - 2*M, FMM = float(MM)
    let m: i32 = match kind {
        HankelKind::First => 1,
        HankelKind::Second => 2,
    };
    let mm = 3 - 2 * m; // +1 for m=1, -1 for m=2
    let fmm = T::from_f64(mm as f64);

    // Rotated argument: ZN = z * exp(-mp) where mp = mm*HPI*i
    // ZNR = FMM * ZI, ZNI = -FMM * ZR (Fortran lines 212-213)
    let zn = Complex::new(fmm * z.im, -fmm * z.re);

    // ── Range check (Fortran lines 217-225) ──
    let az = zabs(z);
    let aa_tol = T::from_f64(0.5) / tol;
    let bb = T::from_f64(2147483647.0 * 0.5);
    let aa = aa_tol.min(bb);

    if az > aa || fn_val > aa {
        return Err(Error::TotalPrecisionLoss);
    }

    let aa_sqrt = aa.sqrt();
    let precision_warning = az > aa_sqrt || fn_val > aa_sqrt;

    // ── Underflow limit (Fortran line 229-230) ──
    let ufl = T::MACH_TINY * T::from_f64(1.0e3);
    if az < ufl {
        return Err(Error::Overflow);
    }

    // ── Overflow pre-checks and dispatch (Fortran lines 231-268) ──
    let mut nz: usize = 0;
    let mut nn_eff = nn;

    if fnu > fnul {
        // Large order: uniform asymptotic expansions (ZBUNK)
        // Fortran ZBESH label 90-100 (zbsubs.f lines 269-284)
        let mut mr = 0i32;
        let mut zn_call = zn;
        // Check if ZN is in the left half plane (needs analytic continuation)
        if !((zn.re >= zero) && (zn.re != zero || zn.im >= zero || m != 2)) {
            mr = -mm;
            if zn.re == zero && zn.im < zero {
                // Negate ZN to put it in RHP (Fortran lines 278-279)
                zn_call = -zn;
            }
        }
        let nw = zbunk(zn_call, fnu, scaling, mr, &mut y[..nn_eff], tol, elim, alim);
        if nw < 0 {
            return if nw == -1 {
                Err(Error::Overflow)
            } else {
                Err(Error::ConvergenceFailure)
            };
        }
        nz += nw as usize;
    } else {
        // ── Small-to-moderate order overflow checks (Fortran lines 232-249) ──
        if fn_val > one {
            if fn_val > two {
                // fn > 2: ZUOIK overflow/underflow pre-check (Fortran lines 238-248)
                let nuf = zuoik(
                    zn,
                    fnu,
                    scaling,
                    IkFlag::K,
                    &mut y[..nn_eff],
                    tol,
                    elim,
                    alim,
                );
                if nuf < 0 {
                    return Err(Error::Overflow);
                }
                nz += nuf as usize;
                nn_eff -= nuf as usize;
                if nn_eff == 0 {
                    // Fortran zbsubs.f lines 249, 338-344:
                    // IF(NN.EQ.0) GO TO 140 → IF(ZNR.LT.0) GO TO 230 (IERR=2)
                    if zn.re < zero {
                        return Err(Error::Overflow);
                    }
                    // All underflowed — y is already zeroed, return early
                    let status = if precision_warning {
                        Accuracy::Reduced
                    } else {
                        Accuracy::Normal
                    };
                    return Ok((nz, status));
                }
            }

            // Check for overflow when az is very small (Fortran lines 234-237)
            if fn_val <= two && az <= tol {
                let half = T::from_f64(0.5);
                let arg = half * az;
                let aln = -fn_val * arg.ln();
                if aln > elim {
                    return Err(Error::Overflow);
                }
            }
        }

        // ── Main dispatch (Fortran label 70, lines 251-268) ──
        if zn.re < zero || (zn.re == zero && zn.im < zero && m == 2) {
            // Left half-plane of rotated argument: analytic continuation (ZACON)
            let mr = -mm; // Fortran: MR = -MM (zbsubs.f line 263)
            let rl = T::rl();
            let nw = zacon(
                zn,
                fnu,
                scaling,
                mr,
                &mut y[..nn_eff],
                rl,
                fnul,
                tol,
                elim,
                alim,
            )?;
            nz = nw; // Fortran: NZ=NW (zbsubs.f line 267, assignment not accumulation)
        } else {
            // Right half-plane: direct zbknu on rotated argument
            let nw = zbknu(zn, fnu, scaling, &mut y[..nn_eff], tol, elim, alim)?;
            nz += nw;
        }
    };

    // ── Phase post-processing (Fortran label 110, lines 285-336) ──
    // H(m, fnu, z) = -FMM*(i/HPI)*exp(-FMM*fnu*HPI*i)*K(fnu, zn)
    let hpi_t = T::from_f64(HPI);
    let sgn = hpi_t.copysign(-fmm); // sign(HPI, -FMM)

    // Decompose FNU for precision: inu, inuh, ir (Fortran lines 296-299)
    // Safety: fnu is finite and < ~1e15 per upper-interface checks
    let inu = fnu.to_i32().unwrap();
    let inuh = inu / 2;
    let ir = inu - 2 * inuh; // 0 or 1
    let arg = (fnu - T::from_f64((inu - ir) as f64)) * sgn;
    let rhpi = one / sgn;

    // CSGN = RHPI * (-sin(arg) + i*cos(arg)) (Fortran lines 303-304)
    let mut csgn = Complex::new(-rhpi * arg.sin(), rhpi * arg.cos());

    // If inuh is odd, negate CSGN (Fortran lines 305-309)
    if inuh % 2 != 0 {
        csgn = -csgn;
    }

    // ZTI = -FMM (Fortran line 311)
    let zti = -fmm;
    let rtol = one / tol;
    let ascle = ufl * rtol;

    // Apply phase factor to each component (Fortran lines 314-336)
    for item in y.iter_mut().take(nn_eff) {
        // Underflow protection: scale up tiny values (Fortran lines 323-327)
        let (scaled, atol) = if item.re.abs().max(item.im.abs()) <= ascle {
            (*item * rtol, tol)
        } else {
            (*item, one)
        };

        // CY(I) = scaled * CSGN * ATOL (Fortran lines 329-332)
        *item = scaled * csgn * atol;

        // Advance CSGN: CSGN *= (0, ZTI) (Fortran lines 333-335)
        csgn = if zti > zero {
            mul_i(csgn)
        } else {
            mul_neg_i(csgn)
        };
    }

    let status = if precision_warning {
        Accuracy::Reduced
    } else {
        Accuracy::Normal
    };

    Ok((nz, status))
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    // ── Input validation ──

    #[test]
    fn besh_z_zero_returns_error() {
        let z = Complex64::new(0.0, 0.0);
        let mut y = [Complex64::new(0.0, 0.0)];
        assert!(matches!(
            zbesh(z, 0.0, HankelKind::First, Scaling::Unscaled, &mut y),
            Err(Error::InvalidInput)
        ));
    }

    #[test]
    fn besh_negative_order_returns_error() {
        let z = Complex64::new(1.0, 0.0);
        let mut y = [Complex64::new(0.0, 0.0)];
        assert!(matches!(
            zbesh(z, -1.0, HankelKind::First, Scaling::Unscaled, &mut y),
            Err(Error::InvalidInput)
        ));
    }

    #[test]
    fn besh_n_zero_returns_error() {
        let z = Complex64::new(1.0, 0.0);
        let mut y: [Complex64; 0] = [];
        assert!(matches!(
            zbesh(z, 0.0, HankelKind::First, Scaling::Unscaled, &mut y),
            Err(Error::InvalidInput)
        ));
    }

    // ── Basic correctness ──

    #[test]
    fn besh_h1_conjugate_h2_real_axis() {
        // For real x > 0: H^(1)(v, x) = conj(H^(2)(v, x))
        let z = Complex64::new(2.0, 0.0);
        let fnu = 0.5;
        let mut y1 = [Complex64::new(0.0, 0.0)];
        let mut y2 = [Complex64::new(0.0, 0.0)];
        zbesh(z, fnu, HankelKind::First, Scaling::Unscaled, &mut y1).unwrap();
        zbesh(z, fnu, HankelKind::Second, Scaling::Unscaled, &mut y2).unwrap();
        let diff = (y1[0] - y2[0].conj()).norm();
        let scale = y1[0].norm();
        assert!(
            diff / scale < 1e-14,
            "H^(1) != conj(H^(2)) on real axis, err = {:.2e}",
            diff / scale
        );
    }
}
