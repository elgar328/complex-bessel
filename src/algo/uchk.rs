//! Underflow check for scaled complex values.
//!
//! Translation of Fortran ZUCHK from TOMS 644 (zbsubs.f lines 4778-4805).

use num_complex::Complex;

use crate::machine::BesselFloat;

/// Check if a scaled complex value would underflow when unscaled.
///
/// `y` is a scaled quantity whose magnitude is greater than `ascle`.
/// This function tests whether scaling `y` by `tol` (to its proper value)
/// would cause the smaller component to underflow relative to the larger.
///
/// Returns `true` if the value is deemed underflowed (should be set to zero),
/// `false` if it is acceptable.
///
/// The test: if `min(|re|, |im|) ≤ ascle` and `max(|re|, |im|) < min/tol`,
/// then the phase angle has lost all accuracy and underflow is assumed.
///
/// Equivalent to Fortran ZUCHK (NZ=1 → true, NZ=0 → false).
#[inline]
pub(crate) fn zuchk<T: BesselFloat>(y: Complex<T>, ascle: T, tol: T) -> bool {
    let wr = y.re.abs();
    let wi = y.im.abs();
    let st = wr.min(wi);
    if st > ascle {
        return false;
    }
    let ss = wr.max(wi);
    // st/tol: what the smaller component magnitude would be at proper scale
    ss < st / tol
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn uchk_large_value_not_underflow() {
        // Both components well above ascle → no underflow
        let y = Complex64::new(1.0, 1.0);
        assert!(!zuchk(y, 1e-300, 1e-16));
    }

    #[test]
    fn uchk_tiny_value_underflow() {
        // Both components tiny, max < min/tol → underflow
        let y = Complex64::new(1e-310, 1e-320);
        assert!(zuchk(y, 1e-300, 1e-16));
    }

    #[test]
    fn uchk_mixed_no_underflow() {
        // One component large, one small but min > ascle → no underflow
        let y = Complex64::new(1.0, 1e-10);
        assert!(!zuchk(y, 1e-300, 1e-16));
    }

    #[test]
    fn uchk_at_boundary() {
        // min component exactly at ascle: st <= ascle, check if ss < st/tol
        let ascle = 1e-300_f64;
        let tol = 1e-16_f64;
        // st = ascle, st/tol = 1e-284, ss must be < 1e-284 for underflow
        let y = Complex64::new(1e-285, ascle);
        assert!(zuchk(y, ascle, tol));
    }
}
