//! Y Bessel function upper interface.
//!
//! Translation of Fortran ZBESY from TOMS 644 (zbsubs.f lines 1177-1461).
//! Y(fnu, z) = i*CC*I(fnu, arg) - (2/pi)*conj(CC)*K(fnu, arg)
//! where CC = exp(-i*pi*fnu/2), arg = z*exp(-i*pi/2).

#[cfg(all(feature = "alloc", not(feature = "std")))]
use alloc::vec;

use num_complex::Complex;

use crate::algo::constants::HPI;
use crate::besi::zbesi;
use crate::besk::zbesk;
use crate::machine::BesselFloat;
use crate::types::{BesselError, BesselStatus, Scaling};

/// Compute Y_{fnu+j}(z) for j = 0, 1, ..., n-1.
///
/// Results are written into the `y` slice. `n` is derived from `y.len()`.
///
/// Uses the identity Y(v,z) = i*CC*I(v,arg) - (2/pi)*conj(CC)*K(v,arg)
/// where CC = exp(-i*pi*v/2) and arg = z*exp(-i*pi/2).
///
/// Equivalent to Fortran ZBESY in TOMS 644.
pub(crate) fn zbesy<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    scaling: Scaling,
    y: &mut [Complex<T>],
) -> Result<(usize, BesselStatus), BesselError> {
    let n = y.len();
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);
    let hpi_t = T::from_f64(HPI);

    // Zero the output buffer
    for yi in y.iter_mut() {
        *yi = czero;
    }

    // Input validation (Fortran lines 1342-1348)
    if n < 1 {
        return Err(BesselError::InvalidInput);
    }
    if fnu < zero {
        return Err(BesselError::InvalidInput);
    }
    if z.re == zero && z.im == zero {
        return Err(BesselError::InvalidInput);
    }

    // Rotate argument: zn = (zz.im, -zz.re) where zz = z with Im >= 0
    // (Fortran lines 1349-1353)
    let zzi = if z.im < zero { -z.im } else { z.im };
    let zzr = z.re;
    let znr = zzi;
    let zni = -zzr;
    let zn = Complex::new(znr, zni);

    // Compute coefficients CC and CSPN (Fortran lines 1359-1373)
    // CIPR = [1, 0, -1, 0], CIPI = [0, 1, 0, -1] (powers of i)
    let cipr: [T; 4] = [one, zero, -one, zero];
    let cipi: [T; 4] = [zero, one, zero, -one];

    let ifnu = fnu.to_i32().unwrap();
    let ffnu = fnu - T::from_f64(ifnu as f64);
    let arg = hpi_t * ffnu;
    let mut csgnr = arg.cos();
    let mut csgni = arg.sin();

    // Multiply by i^(ifnu mod 4) (Fortran lines 1364-1367)
    let i4 = (ifnu % 4) as usize;
    let str = csgnr * cipr[i4] - csgni * cipi[i4];
    csgni = csgnr * cipi[i4] + csgni * cipr[i4];
    csgnr = str;

    let rhpi = one / hpi_t;
    let mut cspnr = csgnr * rhpi;
    let mut cspni = -csgni * rhpi;
    // CSGN *= i (Fortran lines 1371-1373)
    let str2 = -csgni;
    csgni = csgnr;
    csgnr = str2;

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
            let str_val = csgnr * i_buf[0].re
                - csgni * i_buf[0].im
                - (cspnr * k_buf[0].re - cspni * k_buf[0].im);
            let sti_val = csgnr * i_buf[0].im + csgni * i_buf[0].re
                - (cspnr * k_buf[0].im + cspni * k_buf[0].re);

            // Conjugate if original Im(z) < 0 (Fortran lines 1390-1394)
            y[0] = if z.im < zero {
                Complex::new(str_val, -sti_val)
            } else {
                Complex::new(str_val, sti_val)
            };

            return Ok((nz, BesselStatus::Normal));
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
        let str_ex = (exr * cspnr - exi * cspni) * ey;
        cspni = (exr * cspni + exi * cspnr) * ey;
        cspnr = str_ex;

        let rtol = one / tol;
        let ascle = T::MACH_TINY * rtol * T::from_f64(1.0e3);

        // Scale K values if near underflow (Fortran lines 1423-1433)
        let mut zvr = k_buf[0].re;
        let mut zvi = k_buf[0].im;
        let mut atol = one;
        if zvr.abs().max(zvi.abs()) <= ascle {
            zvr = zvr * rtol;
            zvi = zvi * rtol;
            atol = tol;
        }
        let str_zv = (zvr * cspnr - zvi * cspni) * atol;
        let zvi_new = (zvr * cspni + zvi * cspnr) * atol;
        zvr = str_zv;
        zvi = zvi_new;

        // Scale I values if near underflow (Fortran lines 1434-1444)
        let mut zur = i_buf[0].re;
        let mut zui = i_buf[0].im;
        atol = one;
        if zur.abs().max(zui.abs()) <= ascle {
            zur = zur * rtol;
            zui = zui * rtol;
            atol = tol;
        }
        let str_zu = (zur * csgnr - zui * csgni) * atol;
        let zui_new = (zur * csgni + zui * csgnr) * atol;
        zur = str_zu;
        zui = zui_new;

        // CY(I) = ZU - ZV (Fortran lines 1445-1449)
        let cyr = zur - zvr;
        let cyi_val = zui - zvi;

        y[0] = if z.im < zero {
            Complex::new(cyr, -cyi_val)
        } else {
            Complex::new(cyr, cyi_val)
        };

        let nz_out = if cyr == zero && cyi_val == zero && ey == zero {
            1
        } else {
            0
        };

        return Ok((nz_out, BesselStatus::Normal));
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
                let str_val = csgnr * i_buf[i].re
                    - csgni * i_buf[i].im
                    - (cspnr * k_buf[i].re - cspni * k_buf[i].im);
                let sti_val = csgnr * i_buf[i].im + csgni * i_buf[i].re
                    - (cspnr * k_buf[i].im + cspni * k_buf[i].re);
                y[i] = Complex::new(str_val, sti_val);

                // Advance CSGN *= i (rotate by pi/2) (Fortran lines 1383-1385)
                let str_csgn = -csgni;
                csgni = csgnr;
                csgnr = str_csgn;
                // Advance CSPN *= -i (rotate by -pi/2) (Fortran lines 1386-1388)
                let str_cspn = cspni;
                cspni = -cspnr;
                cspnr = str_cspn;
            }

            // Conjugate if original Im(z) < 0 (Fortran lines 1390-1394)
            if z.im < zero {
                for yi in y.iter_mut().take(n) {
                    *yi = Complex::new(yi.re, -yi.im);
                }
            }

            return Ok((nz, BesselStatus::Normal));
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
        let str_ex = (exr * cspnr - exi * cspni) * ey;
        cspni = (exr * cspni + exi * cspnr) * ey;
        cspnr = str_ex;

        let mut nz_out: usize = 0;
        let rtol = one / tol;
        let ascle = T::MACH_TINY * rtol * T::from_f64(1.0e3);

        for i in 0..n {
            // Scale K values if near underflow (Fortran lines 1423-1433)
            let mut zvr = k_buf[i].re;
            let mut zvi = k_buf[i].im;
            let mut atol = one;
            if zvr.abs().max(zvi.abs()) <= ascle {
                zvr = zvr * rtol;
                zvi = zvi * rtol;
                atol = tol;
            }
            let str_zv = (zvr * cspnr - zvi * cspni) * atol;
            let zvi_new = (zvr * cspni + zvi * cspnr) * atol;
            zvr = str_zv;
            zvi = zvi_new;

            // Scale I values if near underflow (Fortran lines 1434-1444)
            let mut zur = i_buf[i].re;
            let mut zui = i_buf[i].im;
            atol = one;
            if zur.abs().max(zui.abs()) <= ascle {
                zur = zur * rtol;
                zui = zui * rtol;
                atol = tol;
            }
            let str_zu = (zur * csgnr - zui * csgni) * atol;
            let zui_new = (zur * csgni + zui * csgnr) * atol;
            zur = str_zu;
            zui = zui_new;

            // CY(I) = ZU - ZV (Fortran lines 1445-1449)
            let cyr = zur - zvr;
            let cyi_val = zui - zvi;

            y[i] = if z.im < zero {
                Complex::new(cyr, -cyi_val)
            } else {
                Complex::new(cyr, cyi_val)
            };

            if cyr == zero && cyi_val == zero && ey == zero {
                nz_out += 1;
            }

            // Advance CSGN *= i, CSPN *= -i (Fortran lines 1450-1455)
            let str_csgn = -csgni;
            csgni = csgnr;
            csgnr = str_csgn;
            let str_cspn = cspni;
            cspni = -cspnr;
            cspnr = str_cspn;
        }

        Ok((nz_out, BesselStatus::Normal))
    }

    #[cfg(not(feature = "alloc"))]
    {
        // This path should never be reached in practice:
        // - Public single-value functions always use n==1
        // - _seq functions are #[cfg(feature = "alloc")]
        Err(BesselError::InvalidInput)
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
}
