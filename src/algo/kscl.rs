//! K function scaling to handle underflow in the recurrence.
//!
//! Translation of Fortran ZKSCL from TOMS 644 (zbsubs.f lines 2960-3081).

use num_complex::Complex;

use crate::algo::uchk::zuchk;
use crate::machine::BesselFloat;
use crate::utils::{mul_add, zabs};

/// Rescale K function values that may underflow.
///
/// Sets K function values to zero on underflow. Continues the three-term
/// recurrence on scaled (exponentially large) values until two consecutive
/// members come on scale, then returns with the values scaled by `1/tol`.
///
/// Equivalent to Fortran ZKSCL in TOMS 644.
///
/// Parameters:
/// - `zr`: the complex argument z
/// - `fnu`: starting order ν
/// - `y`: K function values for orders ν, ν+1, ...; modified in place
/// - `rz`: 2/z, used in recurrence
/// - `ascle`: underflow threshold
/// - `tol`: tolerance
/// - `elim`: underflow elimination threshold
///
/// Returns `nz`: number of underflowed (zeroed) leading components.
pub(crate) fn zkscl<T: BesselFloat>(
    zr: Complex<T>,
    fnu: T,
    y: &mut [Complex<T>],
    rz: Complex<T>,
    ascle: T,
    tol: T,
    elim: T,
) -> usize {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);
    let n = y.len();
    if n == 0 {
        return 0;
    }

    let mut nz: usize = 0;
    // ic: last 1-based index whose value was successfully unscaled (0 = none)
    let mut ic: usize = 0;

    // Save copies of first min(2,N) values for later recurrence
    let nn = n.min(2);
    let mut cy = [czero; 2];

    // ── Phase 1: Try to unscale the first 1 or 2 values ──
    for i in 0..nn {
        let s1 = y[i];
        cy[i] = s1;
        let az = zabs(s1);
        let acs = -zr.re + az.ln();
        nz += 1;
        y[i] = czero;

        if acs >= -elim {
            // Unscale: exp(ln(s1) - zr) / tol
            let cs = s1.ln() - zr;
            let cs_scaled = cs.exp() / tol;
            if !zuchk(cs_scaled, ascle, tol) {
                y[i] = cs_scaled;
                ic = i + 1; // 1-based
                nz -= 1;
            }
        }
    }

    if n == 1 {
        return nz;
    }

    // If the second value (ic>1) was NOT on scale, zero first and set nz=2
    if ic <= 1 {
        y[0] = czero;
        nz = 2;
    }

    if n == 2 || nz == 0 {
        return nz;
    }

    // ── Phase 2: Three-term recurrence K_{ν+i}(z) ──
    // Recurrence: K_{ν+k+1}(z) = ck·K_{ν+k}(z) + K_{ν+k-1}(z)
    // where ck = (ν+k+1) · rz
    let fn_val = fnu + one;
    let mut ck = rz * fn_val;
    let mut s1 = cy[0];
    let mut s2 = cy[1];
    let helim = T::from_f64(0.5) * elim;
    let celmr = (-elim).exp();
    let mut zd = zr;

    let mut found_two = false;
    let mut kk_found = 0_usize; // 1-based index where two consecutive were found

    #[allow(clippy::needless_range_loop)]
    for i in 2..n {
        let kk = i + 1; // 1-based index (Fortran KK = I, where I starts at 3)

        // Recurrence step: s2_new = ck*s2 + s1
        let cs = s2;
        s2 = mul_add(ck, cs, s1);
        s1 = cs;
        ck = ck + rz;

        let az = zabs(s2);
        let alas = az.ln();
        let acs = -zd.re + alas;
        y[i] = czero;

        if acs >= -elim {
            // Try to unscale
            let cs_log = s2.ln() - zd;
            let cs_scaled = cs_log.exp() / tol;
            if !zuchk(cs_scaled, ascle, tol) {
                y[i] = cs_scaled;
                // Check if previous was also on scale: IC == KK-1 (both 1-based)
                if ic == kk - 1 {
                    found_two = true;
                    kk_found = kk;
                    break;
                }
                ic = kk;
                continue;
            }
        }

        // Value underflowed; if recurrence values are getting too large, rescale
        if alas >= helim {
            zd = zd - elim;
            s1 = s1 * celmr;
            s2 = s2 * celmr;
        }
    }

    if found_two {
        // Fortran: NZ = KK - 2 (KK is 1-based)
        nz = kk_found - 2;
    } else {
        // Loop completed without finding two consecutive on-scale values
        nz = n;
        if ic == n {
            nz = n - 1;
        }
    }

    // Zero out the leading NZ values
    for val in y.iter_mut().take(nz) {
        *val = czero;
    }

    nz
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn kscl_single_element_large() {
        // Single element with large scaled value
        let zr = Complex64::new(1.0, 0.0);
        let rz = Complex64::new(2.0, 0.0);
        let mut y = [Complex64::new(1e100, 0.0)];
        let nz = zkscl(zr, 0.5, &mut y, rz, 1e-300, 1e-16, 700.0);
        // After unscaling: exp(ln(1e100) - 1) / 1e-16 = exp(229.26) / 1e-16
        // This is huge, should be on scale
        assert_eq!(nz, 0);
        assert!(y[0] != Complex64::new(0.0, 0.0));
    }

    #[test]
    fn kscl_all_underflow() {
        // Tiny values with large zr so acs = -zr.re + ln(|s1|) < -elim
        // acs = -700 + ln(1e-100) = -700 - 230 = -930 < -700 → underflow
        let zr = Complex64::new(700.0, 0.0);
        let rz = Complex64::new(2.0 / 700.0, 0.0);
        let mut y = [Complex64::new(1e-100, 0.0), Complex64::new(1e-100, 0.0)];
        let nz = zkscl(zr, 0.5, &mut y, rz, 1e-290, 1e-16, 700.0);
        assert_eq!(nz, 2);
        assert_eq!(y[0], Complex64::new(0.0, 0.0));
        assert_eq!(y[1], Complex64::new(0.0, 0.0));
    }

    #[test]
    fn kscl_empty() {
        let zr = Complex64::new(1.0, 0.0);
        let rz = Complex64::new(2.0, 0.0);
        let mut y: [Complex64; 0] = [];
        let nz = zkscl(zr, 0.5, &mut y, rz, 1e-300, 1e-16, 700.0);
        assert_eq!(nz, 0);
    }
}
