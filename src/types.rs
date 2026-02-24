//! Core types for Bessel function computation.

#[cfg(feature = "alloc")]
use alloc::vec::Vec;
use core::fmt;

use num_complex::Complex;

use crate::machine::BesselFloat;

/// Status of the computation result.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Accuracy {
    /// Computation within normal precision bounds (full machine precision).
    Normal,
    /// Result computed but may have lost more than half of significant digits.
    ///
    /// Triggered when:
    /// - Bessel functions: |z| or ν exceeds ~32767 (f64)
    /// - Airy functions: |z| exceeds ~1024 (f64)
    Reduced,
}

/// Result of an Airy function computation, returned by `_raw` functions
/// (e.g., [`airy_raw`](crate::airy_raw)).
///
/// Single-value convenience functions (`airy`, `biry`, …) do not expose
/// this type; they return only the computed value and discard the status.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AiryResult<T: BesselFloat> {
    /// Computed function value.
    pub value: Complex<T>,
    /// Precision status of the computation.
    ///
    /// [`Accuracy::Reduced`] indicates that |z| is large enough
    /// for more than half of significant digits to be lost.
    pub status: Accuracy,
}

/// Result of a sequence computation, returned by `_seq` functions
/// (e.g., [`besselk_seq`](crate::besselk_seq)).
///
/// Single-value convenience functions (`besselj`, `besselk`, …) do not expose
/// this type; they return only the computed value and discard the status.
#[cfg(feature = "alloc")]
#[derive(Debug, Clone, PartialEq)]
pub struct BesselResult<T: BesselFloat> {
    /// Computed function values for orders ν, ν+1, ..., ν+n-1.
    pub values: Vec<Complex<T>>,
    /// Number of leading components set to zero due to underflow.
    pub underflow_count: usize,
    /// Precision status of the computation.
    ///
    /// Single-value convenience functions do not expose this status;
    /// use a `_seq` function to inspect it when needed.
    pub status: Accuracy,
}

/// Scaling option for Bessel function computation.
///
/// The `_scaled` variant returns `factor · f(z)`, where factor is:
/// - J, Y: `exp(-|Im(z)|)`
/// - I: `exp(-|Re(z)|)`
/// - K: `exp(z)`
/// - H^(1): `exp(-iz)`
/// - H^(2): `exp(iz)`
/// - Ai, Ai': `exp(ζ)` where `ζ = (2/3) · z · √z`
/// - Bi, Bi': `exp(-|Re(ζ)|)` where `ζ = (2/3) · z · √z`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Scaling {
    /// No scaling applied.
    Unscaled,
    /// Exponential scaling to prevent overflow/underflow.
    Exponential,
}

/// Kind of Hankel function.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum HankelKind {
    /// H^(1), Hankel function of the first kind.
    First,
    /// H^(2), Hankel function of the second kind.
    Second,
}

/// Selects I or K function path in uniform asymptotic expansion.
///
/// Controls the prefactor (1/sqrt(2π) vs sqrt(π/2)) and whether the
/// asymptotic sum uses straight or alternating signs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum IkFlag {
    /// I function (Fortran ikflg=1): straight sum, sqrt(1/2π) prefactor.
    I,
    /// K function (Fortran ikflg=2): alternating sum, sqrt(π/2) prefactor.
    K,
}

/// Controls whether the full asymptotic sum is computed or only phi/zeta.
///
/// When only the leading-order overflow/underflow test is needed (e.g. in
/// ZUOIK), the expensive sum computation can be skipped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SumOption {
    /// Compute all parameters including the asymptotic sum (Fortran ipmtr=0).
    Full,
    /// Compute phi/zeta only, skip the sum (Fortran ipmtr=1, for overflow pre-check).
    SkipSum,
}

/// Selects Airy function value or its derivative.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AiryDerivative {
    /// Ai(z) or Bi(z).
    Value,
    /// Ai'(z) or Bi'(z).
    Derivative,
}

/// Error type for Bessel function computation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Error {
    /// Invalid input (e.g., z=0 for K/Y/H, n < 1 in sequence functions).
    InvalidInput,
    /// Overflow: |z| or ν too large, or |z| too small.
    Overflow,
    /// Complete loss of significant digits; |z| or ν too large for meaningful computation.
    TotalPrecisionLoss,
    /// Algorithm did not meet termination criteria.
    ConvergenceFailure,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::InvalidInput => {
                write!(f, "invalid input: check z and nu constraints")
            }
            Error::Overflow => {
                write!(f, "overflow: result magnitude exceeds representable range")
            }
            Error::TotalPrecisionLoss => {
                write!(f, "total precision loss: no significant digits remain")
            }
            Error::ConvergenceFailure => {
                write!(
                    f,
                    "convergence failure: algorithm did not meet termination criteria"
                )
            }
        }
    }
}

impl core::error::Error for Error {}
