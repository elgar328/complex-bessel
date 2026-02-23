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

    /// Infallible conversion from f64.
    ///
    /// For f64 this is the identity; for f32 it truncates via `as f32`.
    /// All Amos algorithm constants originate as f64 literals, so this
    /// conversion always succeeds for the supported types.
    fn from_f64(x: f64) -> Self;

    /// Tolerance: max(MACH_EPSILON, 1e-18).
    fn tol() -> Self;
    /// Large order threshold: 10 + 6*(DIG - 3), where DIG = log10(2) * (DIGITS - 1).
    fn fnul() -> Self;
    /// Asymptotic region boundary: 1.2*DIG + 3, where DIG = log10(2) * (DIGITS - 1).
    fn rl() -> Self;
    /// Underflow elimination threshold: 2.303*(K*R1M5 - 3), K = min(|MIN_EXP|, MAX_EXP).
    fn elim() -> Self;
    /// Overflow elimination threshold: ELIM + max(-2.303*R1M5*(DIGITS-1), -41.45).
    fn alim() -> Self;

    /// Fused multiply-add: `self * a + b`.
    ///
    /// With `std` enabled, uses hardware FMA via the C library `fma()`.
    /// Without `std`, falls back to plain `self * a + b` to avoid the
    /// slow software FMA in libm.
    ///
    /// Named `fma` to avoid ambiguity with [`Float::mul_add`].
    fn fma(self, a: Self, b: Self) -> Self;
}

impl BesselFloat for f64 {
    const MACH_EPSILON: f64 = 2.220446049250313e-16;
    const MACH_TINY: f64 = 2.2250738585072014e-308;
    const MACH_HUGE: f64 = 1.7976931348623157e+308;
    const MACH_DIGITS: i32 = 53;
    const MACH_MIN_EXP: i32 = -1021;
    const MACH_MAX_EXP: i32 = 1024;

    #[inline]
    fn from_f64(x: f64) -> f64 {
        x
    }
    #[inline]
    fn tol() -> f64 {
        2.220446049250313e-16
    } // max(EPSILON, 1e-18)
    #[inline]
    fn fnul() -> f64 {
        85.92135864716212
    } // 10 + 6*(DIG - 3)
    #[inline]
    fn rl() -> f64 {
        21.784271729432426
    } // 1.2*DIG + 3
    #[inline]
    fn elim() -> f64 {
        700.9217936944459
    } // 2.303*(K*R1M5 - 3)
    #[inline]
    fn alim() -> f64 {
        664.8716455337102
    } // ELIM + max(-2.303*R1M5*(DIGITS-1), -41.45)

    #[cfg(feature = "std")]
    #[inline]
    fn fma(self, a: f64, b: f64) -> f64 {
        Float::mul_add(self, a, b)
    }

    #[cfg(not(feature = "std"))]
    #[inline]
    fn fma(self, a: f64, b: f64) -> f64 {
        self * a + b
    }
}

// Derived constants are written at full f64 precision to document the exact
// formula results; the compiler rounds to f32 at compile time.
#[allow(clippy::excessive_precision)]
impl BesselFloat for f32 {
    const MACH_EPSILON: f32 = 1.1920929e-7;
    const MACH_TINY: f32 = 1.1754944e-38;
    const MACH_HUGE: f32 = 3.4028235e+38;
    const MACH_DIGITS: i32 = 24;
    const MACH_MIN_EXP: i32 = -125;
    const MACH_MAX_EXP: i32 = 128;

    #[inline]
    fn from_f64(x: f64) -> f32 {
        x as f32
    }
    #[inline]
    fn tol() -> f32 {
        1.1920929e-7
    } // max(EPSILON, 1e-18)
    #[inline]
    fn fnul() -> f32 {
        33.542139401629406
    } // 10 + 6*(DIG - 3)
    #[inline]
    fn rl() -> f32 {
        11.308427880325882
    } // 1.2*DIG + 3
    #[inline]
    fn elim() -> f32 {
        79.75001000176859
    } // 2.303*(K*R1M5 - 3)
    #[inline]
    fn alim() -> f32 {
        63.80475216144317
    } // ELIM + max(-2.303*R1M5*(DIGITS-1), -41.45)

    #[cfg(feature = "std")]
    #[inline]
    fn fma(self, a: f32, b: f32) -> f32 {
        Float::mul_add(self, a, b)
    }

    #[cfg(not(feature = "std"))]
    #[inline]
    fn fma(self, a: f32, b: f32) -> f32 {
        self * a + b
    }
}
