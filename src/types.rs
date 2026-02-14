//! Core types for Bessel function computation.

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;
use core::fmt;

use num_complex::Complex;

use crate::machine::BesselFloat;

/// Result of a Bessel function sequence computation.
pub struct BesselResult<T: BesselFloat> {
    /// Computed function values for orders ν, ν+1, ..., ν+n-1.
    pub values: Vec<Complex<T>>,
    /// Number of leading components set to zero due to underflow.
    pub underflow_count: usize,
}

/// Scaling option for Bessel function computation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scaling {
    /// No scaling applied.
    Unscaled,
    /// Exponential scaling to prevent overflow/underflow.
    Exponential,
}

/// Kind of Hankel function.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HankelKind {
    /// H^(1), Hankel function of the first kind.
    First,
    /// H^(2), Hankel function of the second kind.
    Second,
}

/// Selects Airy function value or its derivative.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AiryDerivative {
    /// Ai(z) or Bi(z).
    Value,
    /// Ai'(z) or Bi'(z).
    Derivative,
}

/// Error type for Bessel function computation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BesselError {
    /// Invalid input (e.g., z=0 for Hankel, negative order).
    InvalidInput,
    /// Overflow: |z| or ν too large, or |z| too small.
    Overflow,
    /// Argument reduction caused loss of more than half of significant digits.
    PrecisionLoss,
    /// Argument reduction caused loss of all significant digits.
    TotalPrecisionLoss,
    /// Algorithm did not meet termination criteria.
    ConvergenceFailure,
}

impl fmt::Display for BesselError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BesselError::InvalidInput => {
                write!(f, "invalid input: check z and nu constraints")
            }
            BesselError::Overflow => {
                write!(f, "overflow: result magnitude exceeds representable range")
            }
            BesselError::PrecisionLoss => {
                write!(
                    f,
                    "precision loss: more than half of significant digits lost"
                )
            }
            BesselError::TotalPrecisionLoss => {
                write!(f, "total precision loss: no significant digits remain")
            }
            BesselError::ConvergenceFailure => {
                write!(
                    f,
                    "convergence failure: algorithm did not meet termination criteria"
                )
            }
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for BesselError {}
