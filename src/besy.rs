//! Y Bessel function upper interface.
//!
//! Translation of Fortran ZBESY from TOMS 644 (zbsubs.f lines 1177-1461).
//! Y(fnu, z) = i*CC*I(fnu, arg) - (2/pi)*conj(CC)*K(fnu, arg)
//! where CC = exp(-i*pi*fnu/2), arg = z*exp(-i*pi/2).

#[cfg(feature = "alloc")]
use alloc::vec;

use num_complex::Complex;

use crate::algo::constants::HPI;
use crate::besi::zbesi;
use crate::besk::zbesk;
use crate::machine::BesselFloat;
use crate::types::{Accuracy, Error, Scaling};
use crate::utils::{mul_i, mul_neg_i, zabs};

/// Compute Y_{fnu+j}(z) for j = 0, 1, ..., n-1.
///
/// Results are written into the `y` slice. `n` is derived from `y.len()`.
///
/// Uses the identity Y(v,z) = i*CC*I(v,arg) - (2/pi)*conj(CC)*K(v,arg)
/// where CC = exp(-i*pi*v/2) and arg = z*exp(-i*pi/2).
///
/// Equivalent to Fortran ZBESY in TOMS 644.
#[inline]
pub(crate) fn zbesy<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    scaling: Scaling,
    y: &mut [Complex<T>],
) -> Result<(usize, Accuracy), Error> {
    let n = y.len();
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);
    let hpi_t = T::from_f64(HPI);

    // Zero the output buffer
    y.fill(czero);

    // Input validation (Fortran lines 1342-1348)
    if n < 1 {
        return Err(Error::InvalidInput);
    }
    if fnu < zero {
        return Err(Error::InvalidInput);
    }
    if z == czero {
        return Err(Error::InvalidInput);
    }
    // Reject non-finite inputs: NaN passes every ordered comparison below,
    // then panics at to_i32() deep in the algorithm chain (e.g. zmlri).
    if !z.re.is_finite() || !z.im.is_finite() || !fnu.is_finite() {
        return Err(Error::InvalidInput);
    }

    // Range check, as in the other upper interfaces. The delegates zbesi/zbesk
    // repeat it, but fnu.to_i32() below must already be protected here.
    let az = zabs(z);
    let fn_val = fnu + T::from_f64((n - 1) as f64);
    let aa_tol = T::from_f64(0.5) / T::tol();
    let bb = T::from_f64(2147483647.0 * 0.5);
    let aa = aa_tol.min(bb);

    if az > aa || fn_val > aa {
        return Err(Error::TotalPrecisionLoss);
    }

    // Partial precision loss (Fortran IERR=3). The delegates zbesi/zbesk hit
    // the same condition (|zn| == |z|), but their statuses are not returned,
    // so compute it here like the other upper interfaces do.
    let aa_sqrt = aa.sqrt();
    let status = if az > aa_sqrt || fn_val > aa_sqrt {
        Accuracy::Reduced
    } else {
        Accuracy::Normal
    };

    // Rotate argument: zn = (zz.im, -zz.re) where zz = z with Im >= 0
    // (Fortran lines 1349-1353)
    let zn = Complex::new(z.im.abs(), -z.re);

    // Compute coefficients CC and CSPN (Fortran lines 1359-1373)
    // Safety: fnu is finite (checked above) and <= aa < i32::MAX per range check
    let ifnu = fnu.to_i32().unwrap();
    let ffnu = fnu - T::from_f64(ifnu as f64);
    let arg = hpi_t * ffnu;
    let mut csgn = Complex::new(arg.cos(), arg.sin());

    // Multiply by i^(ifnu mod 4) (Fortran lines 1364-1367)
    csgn = match (ifnu % 4) as usize {
        0 => csgn,
        1 => mul_i(csgn),
        2 => -csgn,
        _ => mul_neg_i(csgn),
    };

    let rhpi = one / hpi_t;
    let mut cspn = csgn.conj() * rhpi;
    // CSGN *= i (Fortran lines 1371-1373)
    csgn = mul_i(csgn);

    // For n == 1: use stack arrays (no alloc needed)
    if n == 1 {
        let mut i_buf = [czero];
        let mut k_buf = [czero];

        // Compute I(fnu, zn) via ZBESI (Fortran line 1354)
        let (nz1, _) = zbesi(zn, fnu, scaling, &mut i_buf)?;
        // Compute K(fnu, zn) via ZBESK (Fortran line 1356)
        let (nz2, _) = zbesk(zn, fnu, scaling, &mut k_buf)?;
        let nz = nz1.min(nz2);

        if scaling == Scaling::Unscaled {
            // KODE=1: simple combination (Fortran DO 50, lines 1375-1394)
            let cy_val = csgn * i_buf[0] - cspn * k_buf[0];

            // Conjugate if original Im(z) < 0 (Fortran lines 1390-1394)
            y[0] = if z.im < zero { cy_val.conj() } else { cy_val };

            return Ok((nz, status));
        }

        // KODE=2: scaled version with underflow protection (Fortran lines 1396-1456)
        let tol = T::tol();
        let elim = T::elim();

        let exr = z.re.cos();
        let exi = z.re.sin();
        let mut ey = zero;
        let tay = (z.im + z.im).abs();
        if tay < elim {
            ey = (-tay).exp();
        }

        // Apply exponential scaling to CSPN (Fortran lines 1411-1413)
        cspn = Complex::new(exr, exi) * cspn * ey;

        let rtol = one / tol;
        let ascle = T::MACH_TINY * rtol * T::from_f64(1.0e3);

        // Scale K values if near underflow (Fortran lines 1423-1433)
        let (zv_s, zv_atol) = if k_buf[0].re.abs().max(k_buf[0].im.abs()) <= ascle {
            (k_buf[0] * rtol, tol)
        } else {
            (k_buf[0], one)
        };
        let zv = zv_s * cspn * zv_atol;

        // Scale I values if near underflow (Fortran lines 1434-1444)
        let (zu_s, zu_atol) = if i_buf[0].re.abs().max(i_buf[0].im.abs()) <= ascle {
            (i_buf[0] * rtol, tol)
        } else {
            (i_buf[0], one)
        };
        let zu = zu_s * csgn * zu_atol;

        // CY(I) = ZU - ZV (Fortran lines 1445-1449)
        let cy_val = zu - zv;

        y[0] = if z.im < zero { cy_val.conj() } else { cy_val };

        let nz_out = if cy_val == czero && ey == zero { 1 } else { 0 };

        return Ok((nz_out, status));
    }

    // For n > 1: requires alloc feature
    #[cfg(feature = "alloc")]
    {
        let mut i_buf = vec![czero; n];
        let mut k_buf = vec![czero; n];

        // Compute I(fnu, zn) via ZBESI (Fortran line 1354)
        let (nz1, _) = zbesi(zn, fnu, scaling, &mut i_buf)?;
        // Compute K(fnu, zn) via ZBESK (Fortran line 1356)
        let (nz2, _) = zbesk(zn, fnu, scaling, &mut k_buf)?;
        let nz = nz1.min(nz2);

        if scaling == Scaling::Unscaled {
            // KODE=1: simple combination (Fortran DO 50, lines 1375-1394)
            for i in 0..n {
                // CY(I) = CSGN*I(I) - CSPN*K(I)
                let cy_val = csgn * i_buf[i] - cspn * k_buf[i];
                y[i] = cy_val;

                // Advance CSGN *= i (rotate by pi/2) (Fortran lines 1383-1385)
                csgn = mul_i(csgn);
                // Advance CSPN *= -i (rotate by -pi/2) (Fortran lines 1386-1388)
                cspn = mul_neg_i(cspn);
            }

            // Conjugate if original Im(z) < 0 (Fortran lines 1390-1394)
            if z.im < zero {
                for yi in y.iter_mut().take(n) {
                    *yi = yi.conj();
                }
            }

            return Ok((nz, status));
        }

        // KODE=2: scaled version with underflow protection (Fortran lines 1396-1456)
        let tol = T::tol();
        let elim = T::elim();

        let exr = z.re.cos();
        let exi = z.re.sin();
        let mut ey = zero;
        let tay = (z.im + z.im).abs();
        if tay < elim {
            ey = (-tay).exp();
        }

        // Apply exponential scaling to CSPN (Fortran lines 1411-1413)
        cspn = Complex::new(exr, exi) * cspn * ey;

        let mut nz_out: usize = 0;
        let rtol = one / tol;
        let ascle = T::MACH_TINY * rtol * T::from_f64(1.0e3);

        for i in 0..n {
            // Scale K values if near underflow (Fortran lines 1423-1433)
            let (zv_s, zv_atol) = if k_buf[i].re.abs().max(k_buf[i].im.abs()) <= ascle {
                (k_buf[i] * rtol, tol)
            } else {
                (k_buf[i], one)
            };
            let zv = zv_s * cspn * zv_atol;

            // Scale I values if near underflow (Fortran lines 1434-1444)
            let (zu_s, zu_atol) = if i_buf[i].re.abs().max(i_buf[i].im.abs()) <= ascle {
                (i_buf[i] * rtol, tol)
            } else {
                (i_buf[i], one)
            };
            let zu = zu_s * csgn * zu_atol;

            // CY(I) = ZU - ZV (Fortran lines 1445-1449)
            let cy_val = zu - zv;

            y[i] = if z.im < zero { cy_val.conj() } else { cy_val };

            if cy_val == czero && ey == zero {
                nz_out += 1;
            }

            // Advance CSGN *= i, CSPN *= -i (Fortran lines 1450-1455)
            csgn = mul_i(csgn);
            cspn = mul_neg_i(cspn);
        }

        Ok((nz_out, status))
    }

    #[cfg(not(feature = "alloc"))]
    {
        // This path should never be reached in practice:
        // - Public single-value functions always use n==1
        // - _seq functions are #[cfg(feature = "alloc")]
        Err(Error::InvalidInput)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn besy_z_zero_returns_error() {
        let z = Complex64::new(0.0, 0.0);
        let mut y = [Complex64::new(0.0, 0.0)];
        assert!(zbesy(z, 0.0, Scaling::Unscaled, &mut y).is_err());
    }

    #[test]
    fn besy_non_finite_inputs_return_invalid_input() {
        let nan = f64::NAN;
        let inf = f64::INFINITY;
        let mut y = [Complex64::new(0.0, 0.0)];
        for z in [
            Complex64::new(nan, 0.0),
            Complex64::new(0.0, nan),
            Complex64::new(nan, nan),
            Complex64::new(inf, 0.0),
        ] {
            assert!(matches!(
                zbesy(z, 1.0, Scaling::Unscaled, &mut y),
                Err(Error::InvalidInput)
            ));
        }
        let z = Complex64::new(1.0, 1.0);
        for fnu in [nan, inf] {
            assert!(matches!(
                zbesy(z, fnu, Scaling::Unscaled, &mut y),
                Err(Error::InvalidInput)
            ));
        }
    }

    #[test]
    fn besy_reports_reduced_accuracy_for_large_argument() {
        // |z| > sqrt(aa) ≈ 32767 (f64) means more than half of the significant
        // digits may be lost (Fortran IERR=3). zbesy used to discard the
        // delegate statuses and always report Normal.
        let mut y = [Complex64::new(0.0, 0.0)];
        let (_, status) =
            zbesy(Complex64::new(40000.0, 0.0), 0.0, Scaling::Unscaled, &mut y).unwrap();
        assert_eq!(status, Accuracy::Reduced);

        let (_, status) = zbesy(Complex64::new(1.0, 1.0), 0.0, Scaling::Unscaled, &mut y).unwrap();
        assert_eq!(status, Accuracy::Normal);
    }

    #[test]
    fn besy_large_finite_order_returns_total_precision_loss() {
        // fnu >= 2^31 used to panic at fnu.to_i32().unwrap() because zbesy
        // had no range check of its own before delegating to zbesi/zbesk.
        let z = Complex64::new(1.0, 1.0);
        let mut y = [Complex64::new(0.0, 0.0)];
        assert!(matches!(
            zbesy(z, 3.0e9, Scaling::Unscaled, &mut y),
            Err(Error::TotalPrecisionLoss)
        ));
    }
}
