//! Machine constants and the `BesselFloat` trait.
//!
//! Constants follow the Fortran I1MACH/D1MACH conventions from TOMS 644.

use num_traits::Float;

/// Floating-point trait for Bessel function computation.
///
/// Implemented for `f64` and `f32`. Provides machine constants and
/// derived thresholds used by the Amos algorithm.
pub trait BesselFloat: Float + core::fmt::Debug + 'static {
    /// Machine epsilon (D1MACH(3)).
    const MACH_EPSILON: Self;
    /// Smallest positive normal number (D1MACH(1)).
    const MACH_TINY: Self;
    /// Largest representable number (D1MACH(2)).
    const MACH_HUGE: Self;
    /// Number of binary digits in the mantissa (I1MACH(14)).
    const MACH_DIGITS: i32;
    /// Minimum binary exponent (I1MACH(12)).
    const MACH_MIN_EXP: i32;
    /// Maximum binary exponent (I1MACH(11)).
    const MACH_MAX_EXP: i32;

    /// Tolerance: max(MACH_EPSILON, 1e-18).
    fn tol() -> Self;
    /// Decimal digits: log10(2) * (DIGITS - 1).
    fn dig() -> Self;
    /// Large order threshold: 10 + 6*(DIG - 3).
    fn fnul() -> Self;
    /// Asymptotic region boundary: 1.2*DIG + 3.
    fn rl() -> Self;
    /// Underflow elimination threshold: 2.303*(K*R1M5 - 3), K = min(|MIN_EXP|, MAX_EXP).
    fn elim() -> Self;
    /// Overflow elimination threshold: ELIM + max(-2.303*R1M5*(DIGITS-1), -41.45).
    fn alim() -> Self;
}

impl BesselFloat for f64 {
    const MACH_EPSILON: f64 = 2.220446049250313e-16;
    const MACH_TINY: f64 = 2.2250738585072014e-308;
    const MACH_HUGE: f64 = 1.7976931348623157e+308;
    const MACH_DIGITS: i32 = 53;
    const MACH_MIN_EXP: i32 = -1021;
    const MACH_MAX_EXP: i32 = 1024;

    #[inline]
    fn tol() -> f64 { 2.220446049250313e-16 }       // max(EPSILON, 1e-18)
    #[inline]
    fn dig() -> f64 { 15.653559774527022 }           // R1M5 * (DIGITS - 1)
    #[inline]
    fn fnul() -> f64 { 85.92135864716212 }           // 10 + 6*(DIG - 3)
    #[inline]
    fn rl() -> f64 { 21.784271729432426 }            // 1.2*DIG + 3
    #[inline]
    fn elim() -> f64 { 700.9217936944459 }           // 2.303*(K*R1M5 - 3)
    #[inline]
    fn alim() -> f64 { 664.8716455337102 }           // ELIM + max(-2.303*R1M5*(DIGITS-1), -41.45)
}

impl BesselFloat for f32 {
    const MACH_EPSILON: f32 = 1.1920929e-7;
    const MACH_TINY: f32 = 1.1754944e-38;
    const MACH_HUGE: f32 = 3.4028235e+38;
    const MACH_DIGITS: i32 = 24;
    const MACH_MIN_EXP: i32 = -125;
    const MACH_MAX_EXP: i32 = 128;

    #[inline]
    fn tol() -> f32 { 1.1920929e-7 }                // max(EPSILON, 1e-18)
    #[inline]
    fn dig() -> f32 { 6.923689900271568 }            // R1M5 * (DIGITS - 1)
    #[inline]
    fn fnul() -> f32 { 33.542139401629406 }          // 10 + 6*(DIG - 3)
    #[inline]
    fn rl() -> f32 { 11.308427880325882 }            // 1.2*DIG + 3
    #[inline]
    fn elim() -> f32 { 79.75001000176859 }           // 2.303*(K*R1M5 - 3)
    #[inline]
    fn alim() -> f32 { 63.80475216144317 }           // ELIM + max(-2.303*R1M5*(DIGITS-1), -41.45)
}
