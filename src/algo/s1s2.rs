//! Stokes multiplier normalization for analytic continuation.
//!
//! Translation of Fortran ZS1S2 from TOMS 644 (zbsubs.f lines 3237-3286).

use num_complex::Complex;

use crate::machine::BesselFloat;
use crate::utils::zabs;

/// Result of the Stokes multiplier normalization.
#[derive(Debug, Clone, Copy)]
pub(crate) struct S1S2Output<T> {
    /// Rescaled K-function component.
    pub s1: Complex<T>,
    /// I-function component (possibly zeroed on underflow).
    pub s2: Complex<T>,
    /// Underflow flag: 1 if both s1 and s2 underflowed, 0 otherwise.
    pub nz: i32,
    /// Updated underflow counter.
    pub iuf: i32,
}

/// Test and adjust for underflow in the I+K analytic continuation formula.
///
/// In the analytic continuation, s1 (the K function) is rescaled by
/// `exp(-2·zr)` using logarithmic arithmetic to avoid overflow. If the
/// rescaled result and s2 (the I function) are both below `ascle`,
/// both are set to zero (underflow).
///
/// Equivalent to Fortran ZS1S2 in TOMS 644.
///
/// Parameters:
/// - `zr`: the complex argument z
/// - `s1`: K-function value (unscaled)
/// - `s2`: I-function value (scaled)
/// - `ascle`: underflow threshold (typically 1e3 * TINY / TOL)
/// - `alim`: overflow elimination threshold
/// - `iuf`: running underflow counter from caller
pub(crate) fn zs1s2<T: BesselFloat>(
    zr: Complex<T>,
    s1: Complex<T>,
    s2: Complex<T>,
    ascle: T,
    alim: T,
    iuf: i32,
) -> S1S2Output<T> {
    let zero = T::zero();
    let czero = Complex::new(zero, zero);

    let mut s1_out = s1;
    let mut iuf_out = iuf;
    let mut as1 = zabs(s1);
    let as2 = zabs(s2);

    // If s1 is nonzero, attempt to rescale: s1 = s1 * exp(-2*zr)
    if s1 != czero && as1 != zero {
        // Quick magnitude check via real part:
        // ALN = -2*Re(zr) + ln(|s1|)
        let aln = -zr.re - zr.re + as1.ln();
        let s1d = s1; // save original s1
        s1_out = czero;
        as1 = zero;

        if aln >= -alim {
            // Compute s1 = exp(ln(s1d) - 2*zr) using complex logarithm
            let c1 = s1d.ln() - (zr + zr);
            s1_out = c1.exp();
            as1 = zabs(s1_out);
            iuf_out += 1;
        }
        // else: ALN < -ALIM → s1 underflows, stays zero
    }

    // If both components are below the underflow threshold, zero everything
    let aa = as1.max(as2);
    let s2_out = if aa > ascle {
        s2
    } else {
        s1_out = czero;
        return S1S2Output {
            s1: s1_out,
            s2: czero,
            nz: 1,
            iuf: 0,
        };
    };

    S1S2Output {
        s1: s1_out,
        s2: s2_out,
        nz: 0,
        iuf: iuf_out,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn s1s2_both_zero_input() {
        let zr = Complex64::new(1.0, 0.5);
        let s1 = Complex64::new(0.0, 0.0);
        let s2 = Complex64::new(0.0, 0.0);
        let result = zs1s2(zr, s1, s2, 1e-300, 700.0, 0);
        // Both zero and below ascle → underflow
        assert_eq!(result.nz, 1);
        assert_eq!(result.s1, Complex64::new(0.0, 0.0));
        assert_eq!(result.s2, Complex64::new(0.0, 0.0));
    }

    #[test]
    fn s1s2_large_values_no_underflow() {
        let zr = Complex64::new(0.1, 0.0);
        let s1 = Complex64::new(10.0, 5.0);
        let s2 = Complex64::new(20.0, 10.0);
        let result = zs1s2(zr, s1, s2, 1e-300, 700.0, 0);
        // Large values → no underflow
        assert_eq!(result.nz, 0);
        assert_eq!(result.iuf, 1); // s1 was nonzero and rescaled
    }

    #[test]
    fn s1s2_s1_underflows() {
        // s1 is tiny and -2*Re(zr) + ln(|s1|) < -ALIM → s1 zeroed
        let zr = Complex64::new(500.0, 0.0);
        let s1 = Complex64::new(1e-10, 0.0);
        let s2 = Complex64::new(1.0, 0.0);
        let result = zs1s2(zr, s1, s2, 1e-300, 700.0, 0);
        // -2*500 + ln(1e-10) = -1000 - 23 = -1023 < -700 → s1 zeroed
        assert_eq!(result.s1, Complex64::new(0.0, 0.0));
        // s2 is large enough → no overall underflow
        assert_eq!(result.nz, 0);
        assert_eq!(result.iuf, 0); // iuf not incremented when s1 underflows
    }

    #[test]
    fn s1s2_rescale_s1() {
        // Test that s1 is correctly rescaled by exp(-2*zr)
        let zr = Complex64::new(0.5, 0.0);
        let s1 = Complex64::new(1.0, 0.0);
        let s2 = Complex64::new(1.0, 0.0);
        let result = zs1s2(zr, s1, s2, 1e-300, 700.0, 0);
        // s1_new = exp(ln(1) - 2*0.5) = exp(-1) ≈ 0.3679
        assert!((result.s1.re - (-1.0_f64).exp()).abs() < 1e-14);
        assert!(result.s1.im.abs() < 1e-14);
        assert_eq!(result.nz, 0);
        assert_eq!(result.iuf, 1);
    }
}
