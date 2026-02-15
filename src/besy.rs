//! Y Bessel function upper interface.
//!
//! Translation of Fortran ZBESY from TOMS 644 (zbsubs.f lines 1177-1461).
//! Y(fnu, z) = i*CC*I(fnu, arg) - (2/pi)*conj(CC)*K(fnu, arg)
//! where CC = exp(-i*pi*fnu/2), arg = z*exp(-i*pi/2).

#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]

use num_complex::Complex;

use crate::besi::zbesi;
use crate::besk::zbesk;
use crate::machine::BesselFloat;
use crate::types::{BesselError, BesselResult, Scaling};

const HPI: f64 = 1.57079632679489662; // pi/2

/// Compute Y_{fnu+j}(z) for j = 0, 1, ..., n-1.
///
/// Uses the identity Y(v,z) = i*CC*I(v,arg) - (2/pi)*conj(CC)*K(v,arg)
/// where CC = exp(-i*pi*v/2) and arg = z*exp(-i*pi/2).
///
/// Equivalent to Fortran ZBESY in TOMS 644.
pub(crate) fn zbesy<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    scaling: Scaling,
    n: usize,
) -> Result<BesselResult<T>, BesselError> {
    let zero = T::zero();
    let one = T::one();
    let _half = T::from(0.5).unwrap();
    let hpi_t = T::from(HPI).unwrap();

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

    // Compute I(fnu, zn) via ZBESI (Fortran line 1354)
    let i_result = zbesi(zn, fnu, scaling, n)?;
    let nz1 = i_result.underflow_count;

    // Compute K(fnu, zn) via ZBESK (Fortran line 1356)
    let k_result = zbesk(zn, fnu, scaling, n)?;
    let nz2 = k_result.underflow_count;

    let nz = nz1.min(nz2);

    // Compute coefficients CC and CSPN (Fortran lines 1359-1373)
    // CIPR = [1, 0, -1, 0], CIPI = [0, 1, 0, -1] (powers of i)
    let cipr: [T; 4] = [one, zero, -one, zero];
    let cipi: [T; 4] = [zero, one, zero, -one];

    let ifnu = fnu.to_i32().unwrap();
    let ffnu = fnu - T::from(ifnu as f64).unwrap();
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

    let mut cy = vec![Complex::new(zero, zero); n];

    if scaling == Scaling::Unscaled {
        // KODE=1: simple combination (Fortran DO 50, lines 1375-1394)
        for (i, cy_item) in cy.iter_mut().enumerate().take(n) {
            // CY(I) = CSGN*I(I) - CSPN*K(I)
            let str_val = csgnr * i_result.values[i].re
                - csgni * i_result.values[i].im
                - (cspnr * k_result.values[i].re - cspni * k_result.values[i].im);
            let sti_val = csgnr * i_result.values[i].im + csgni * i_result.values[i].re
                - (cspnr * k_result.values[i].im + cspni * k_result.values[i].re);
            *cy_item = Complex::new(str_val, sti_val);

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
            for cy_item in cy.iter_mut().take(n) {
                *cy_item = Complex::new(cy_item.re, -cy_item.im);
            }
        }

        return Ok(BesselResult {
            values: cy,
            underflow_count: nz,
        });
    }

    // KODE=2: scaled version with underflow protection (Fortran lines 1396-1456)
    let tol = T::tol();
    let k1_abs = T::MACH_MIN_EXP.abs();
    let k2_abs = T::MACH_MAX_EXP.abs();
    let k_val = if k1_abs < k2_abs {
        T::MACH_MIN_EXP.abs()
    } else {
        T::MACH_MAX_EXP.abs()
    };
    let d1m5 = T::from(0.30102999566398120).unwrap(); // log10(2)
    let elim =
        T::from(2.303).unwrap() * (T::from(k_val as f64).unwrap() * d1m5 - T::from(3.0).unwrap());

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
    let ascle = T::MACH_TINY * rtol * T::from(1.0e3).unwrap();

    for (i, cy_item) in cy.iter_mut().enumerate().take(n) {
        // Scale K values if near underflow (Fortran lines 1423-1433)
        let mut zvr = k_result.values[i].re;
        let mut zvi = k_result.values[i].im;
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
        let mut zur = i_result.values[i].re;
        let mut zui = i_result.values[i].im;
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

        *cy_item = if z.im < zero {
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

    Ok(BesselResult {
        values: cy,
        underflow_count: nz_out,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn besy_z_zero_returns_error() {
        let z = Complex64::new(0.0, 0.0);
        assert!(zbesy(z, 0.0, Scaling::Unscaled, 1).is_err());
    }
}
