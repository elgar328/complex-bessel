//! Shared Fortran DATA constants used by multiple algorithm modules.
//!
//! Module-specific constants remain in their respective files.

#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]

/// π (Fortran: DPI, GPI, PI)
pub(crate) const PI: f64 = 3.14159265358979324e+00;

/// π/2 (Fortran: HPI)
pub(crate) const HPI: f64 = 1.57079632679489662e+00;

/// log₁₀(2) = D1MACH(5) for binary IEEE 754 (Fortran: R1M5)
pub(crate) const R1M5: f64 = 0.30102999566398120;

/// 2/3 (Fortran: TTH)
pub(crate) const TTH: f64 = 6.66666666666666667e-01;

/// ln(√(π/2)) (Fortran: AIC in ZUNI2, ZUNK2)
pub(crate) const AIC: f64 = 1.26551212348464539e+00;
