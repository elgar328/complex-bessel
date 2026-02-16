//! Pure Rust implementation of complex Bessel functions based on Amos Algorithm 644 (ACM TOMS 644).
//!
//! Provides Bessel functions J, Y, I, K, Hankel H⁽¹⁾/H⁽²⁾, and Airy functions Ai/Bi
//! for complex arguments and real orders.
//!
//! # Features
//!
//! - **Dual precision** — all functions accept `Complex<f64>` or `Complex<f32>`
//! - **Complete function set** — J, Y, I, K, H⁽¹⁾, H⁽²⁾, Ai, Bi
//! - **Consecutive orders** — `_seq` variants return ν, ν+1, …, ν+n−1 in one call
//! - **Exponential scaling** — `_scaled` variants prevent overflow/underflow
//! - **Negative orders** — single-value functions accept ν < 0 via reflection formulas
//! - **`no_std`** — works with `alloc` only
//!
//! # Quick start
//!
//! ```
//! use complex_bessel::*;
//! use num_complex::Complex;
//!
//! let z = Complex::new(1.0, 2.0);
//!
//! let j = besselj(0.5, z).unwrap();
//! let k = besselk(1.0, z).unwrap();
//! let h = hankel(HankelKind::First, 0.0, z).unwrap();
//!
//! // Scaled versions prevent overflow/underflow
//! let k_scaled = besselk_scaled(1.0, z).unwrap();
//!
//! // Airy functions
//! let ai = airy(z, AiryDerivative::Value).unwrap();
//! ```
//!
//! # Generic types
//!
//! All functions are generic over `f64` and `f32` via the [`BesselFloat`] trait:
//!
//! ```
//! use complex_bessel::besselj;
//! use num_complex::Complex;
//!
//! // f64 (default)
//! let z64 = Complex::new(1.0_f64, 2.0);
//! let j64 = besselj(0.5_f64, z64).unwrap();
//!
//! // f32
//! let z32 = Complex::new(1.0_f32, 2.0);
//! let j32 = besselj(0.5_f32, z32).unwrap();
//! ```
//!
//! # Consecutive orders
//!
//! The `_seq` variants compute values at consecutive orders ν, ν+1, …, ν+n−1
//! in a single call. Internal recurrence is shared, so this is more efficient
//! than calling the single-value function n times.
//!
//! ```
//! use complex_bessel::*;
//! use num_complex::Complex;
//!
//! let z = Complex::new(1.0, 2.0);
//!
//! // K_0(z), K_1(z), K_2(z)
//! let seq = besselk_seq(0.0, z, 3, Scaling::Unscaled).unwrap();
//! assert_eq!(seq.values.len(), 3);
//! ```
//!
//! Sequence results include a [`BesselStatus`] field:
//! - [`BesselStatus::Normal`] — full precision (~14 digits for f64, ~6–7 for f32)
//! - [`BesselStatus::ReducedPrecision`] — some precision lost (|z| or ν very large)
//!
//! Single-value functions silently return the best available result.
//!
//! Sequence variants require ν ≥ 0. Use single-value functions for negative orders.
//!
//! # Exponential scaling
//!
//! The `_scaled` variants multiply by an exponential factor to prevent
//! overflow/underflow for large arguments. See [`Scaling`] for the exact
//! factor applied to each function.
//!
//! ```
//! use complex_bessel::*;
//! use num_complex::Complex;
//!
//! let z = Complex::new(100.0, 0.0);
//!
//! // K_0(100) ≈ 4.66e-45 — unscaled works but close to underflow
//! let k = besselk(0.0, z).unwrap();
//!
//! // exp(100) * K_0(100) ≈ 0.1257 — scaled version stays in normal range
//! let k_s = besselk_scaled(0.0, z).unwrap();
//! ```
//!
//! | Function | Scaled variant returns |
//! |----------|-----------------------|
//! | J, Y | exp(−\|Im(z)\|) · f(z) |
//! | I | exp(−\|Re(z)\|) · I(z) |
//! | K | exp(z) · K(z) |
//! | H<sup>(1)</sup> | exp(−iz) · H<sup>(1)</sup>(z) |
//! | H<sup>(2)</sup> | exp(iz) · H<sup>(2)</sup>(z) |
//! | Ai | exp(ζ) · Ai(z) |
//! | Bi | exp(−\|Re(ζ)\|) · Bi(z) |
//!
//! where ζ = (2/3) z√z.
//!
//! # Negative orders
//!
//! All single-value functions accept any real order, including negative values.
//! DLMF reflection formulas are applied automatically:
//!
//! - **J**: J<sub>−ν</sub>(z) = cos(νπ) J<sub>ν</sub>(z) − sin(νπ) Y<sub>ν</sub>(z) (DLMF 10.4.1)
//! - **Y**: Y<sub>−ν</sub>(z) = sin(νπ) J<sub>ν</sub>(z) + cos(νπ) Y<sub>ν</sub>(z) (DLMF 10.4.2)
//! - **I**: I<sub>−ν</sub>(z) = I<sub>ν</sub>(z) + (2/π) sin(νπ) K<sub>ν</sub>(z) (DLMF 10.27.2)
//! - **K**: K<sub>−ν</sub>(z) = K<sub>ν</sub>(z) (even in ν, DLMF 10.27.3)
//! - **H<sup>(1)</sup>**: H<sup>(1)</sup><sub>−ν</sub>(z) = exp(νπi) H<sup>(1)</sup><sub>ν</sub>(z) (DLMF 10.4.6)
//! - **H<sup>(2)</sup>**: H<sup>(2)</sup><sub>−ν</sub>(z) = exp(−νπi) H<sup>(2)</sup><sub>ν</sub>(z) (DLMF 10.4.7)
//!
//! For integer orders, simplified identities are used (e.g., J<sub>−n</sub>(z) = (−1)<sup>n</sup> J<sub>n</sub>(z)).
//!
//! # `no_std` support
//!
//! Disable the default `std` feature:
//!
//! ```toml
//! [dependencies]
//! complex-bessel = { version = "0.1", default-features = false }
//! ```
//!
//! Requires `alloc`. All functions remain available.

#![cfg_attr(not(feature = "std"), no_std)]

#[cfg(not(feature = "std"))]
extern crate alloc;

pub(crate) mod airy;
pub(crate) mod algo;
pub(crate) mod besh;
pub(crate) mod besi;
pub(crate) mod besj;
pub(crate) mod besk;
pub(crate) mod besy;
pub mod machine;
pub mod types;
pub(crate) mod utils;

pub use machine::BesselFloat;
pub use types::{AiryDerivative, BesselError, BesselResult, BesselStatus, HankelKind, Scaling};

use num_complex::Complex;

// ── Helper: integer order detection ──

/// Check if `nu` is a non-negative integer. Returns `Some(n)` if so.
fn as_integer<T: BesselFloat>(nu: T) -> Option<i64> {
    if nu == nu.floor() {
        // Safe conversion: orders beyond i64 range are not practical
        nu.to_i64()
    } else {
        None
    }
}

// ── Internal: compute with given scaling for negative order support ──

fn besselj_internal<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
) -> Result<Complex<T>, BesselError> {
    let zero = T::zero();
    if nu >= zero {
        let result = besj::zbesj(z, nu, scaling, 1)?;
        return Ok(result.values[0]);
    }

    // Negative order: J_{-ν}(z) = cos(νπ)*J_ν(z) - sin(νπ)*Y_ν(z) (DLMF 10.4.1)
    let abs_nu = nu.abs();

    // Integer shortcut: J_{-n}(z) = (-1)^n * J_n(z)
    if let Some(n) = as_integer(abs_nu) {
        let result = besj::zbesj(z, abs_nu, scaling, 1)?;
        let sign = if n % 2 == 0 { T::one() } else { -T::one() };
        return Ok(result.values[0] * sign);
    }

    // General case: need both J and Y at positive |ν|
    let pi = T::from(core::f64::consts::PI).unwrap();
    let nu_pi = abs_nu * pi;
    let cos_nu_pi = nu_pi.cos();
    let sin_nu_pi = nu_pi.sin();

    let j_result = besj::zbesj(z, abs_nu, scaling, 1)?;
    let y_result = besy::zbesy(z, abs_nu, scaling, 1)?;

    let j_val = j_result.values[0];
    let y_val = y_result.values[0];

    Ok(j_val * cos_nu_pi - y_val * sin_nu_pi)
}

fn bessely_internal<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
) -> Result<Complex<T>, BesselError> {
    let zero = T::zero();
    if nu >= zero {
        let result = besy::zbesy(z, nu, scaling, 1)?;
        return Ok(result.values[0]);
    }

    // Negative order: Y_{-ν}(z) = sin(νπ)*J_ν(z) + cos(νπ)*Y_ν(z) (DLMF 10.4.2)
    let abs_nu = nu.abs();

    // Integer shortcut: Y_{-n}(z) = (-1)^n * Y_n(z)
    if let Some(n) = as_integer(abs_nu) {
        let result = besy::zbesy(z, abs_nu, scaling, 1)?;
        let sign = if n % 2 == 0 { T::one() } else { -T::one() };
        return Ok(result.values[0] * sign);
    }

    // General case: need both J and Y at positive |ν|
    let pi = T::from(core::f64::consts::PI).unwrap();
    let nu_pi = abs_nu * pi;
    let cos_nu_pi = nu_pi.cos();
    let sin_nu_pi = nu_pi.sin();

    let j_result = besj::zbesj(z, abs_nu, scaling, 1)?;
    let y_result = besy::zbesy(z, abs_nu, scaling, 1)?;

    let j_val = j_result.values[0];
    let y_val = y_result.values[0];

    Ok(j_val * sin_nu_pi + y_val * cos_nu_pi)
}

fn besseli_internal<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
) -> Result<Complex<T>, BesselError> {
    let zero = T::zero();
    if nu >= zero {
        let result = besi::zbesi(z, nu, scaling, 1)?;
        return Ok(result.values[0]);
    }

    // Negative order: I_{-ν}(z) = I_ν(z) + (2/π)*sin(νπ)*K_ν(z) (DLMF 10.27.2)
    let abs_nu = nu.abs();

    // Integer shortcut: I_{-n}(z) = I_n(z)
    if as_integer(abs_nu).is_some() {
        let result = besi::zbesi(z, abs_nu, scaling, 1)?;
        return Ok(result.values[0]);
    }

    // General case: need both I and K at positive |ν|
    let pi = T::from(core::f64::consts::PI).unwrap();
    let two = T::from(2.0).unwrap();
    let nu_pi = abs_nu * pi;
    let sin_nu_pi = nu_pi.sin();

    let i_result = besi::zbesi(z, abs_nu, scaling, 1)?;
    let k_result = besk::zbesk(z, abs_nu, scaling, 1)?;

    let i_val = i_result.values[0];
    let k_val = k_result.values[0];

    Ok(i_val + k_val * (two / pi * sin_nu_pi))
}

fn besselk_internal<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
) -> Result<Complex<T>, BesselError> {
    // K_{-ν}(z) = K_ν(z) (DLMF 10.27.3) — K is even in ν
    let abs_nu = nu.abs();
    let result = besk::zbesk(z, abs_nu, scaling, 1)?;
    Ok(result.values[0])
}

fn hankel_internal<T: BesselFloat>(
    kind: HankelKind,
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
) -> Result<Complex<T>, BesselError> {
    let zero = T::zero();
    if nu >= zero {
        let result = besh::zbesh(z, nu, kind, scaling, 1)?;
        return Ok(result.values[0]);
    }

    // Negative order (DLMF 10.4.6, 10.4.7):
    //   H^(1)_{-ν}(z) = exp(νπi) * H^(1)_ν(z)
    //   H^(2)_{-ν}(z) = exp(-νπi) * H^(2)_ν(z)
    let abs_nu = nu.abs();
    let result = besh::zbesh(z, abs_nu, kind, scaling, 1)?;
    let h_val = result.values[0];

    let pi = T::from(core::f64::consts::PI).unwrap();
    let nu_pi = abs_nu * pi;

    let rotation = match kind {
        // exp(νπi) = cos(νπ) + i*sin(νπ)
        HankelKind::First => Complex::new(nu_pi.cos(), nu_pi.sin()),
        // exp(-νπi) = cos(νπ) - i*sin(νπ)
        HankelKind::Second => Complex::new(nu_pi.cos(), -nu_pi.sin()),
    };

    Ok(h_val * rotation)
}

// ── Single-value convenience functions ──

/// Bessel function of the first kind, J_ν(z).
///
/// Computes a single value of the Bessel function J_ν(z) for complex z
/// and real order ν (any real value, including negative).
///
/// For negative ν, the DLMF 10.4.1 reflection formula is applied:
/// `J_{-ν}(z) = cos(νπ) J_ν(z) - sin(νπ) Y_ν(z)`.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails (overflow, precision loss, etc.).
pub fn besselj<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    besselj_internal(nu, z, Scaling::Unscaled)
}

/// Bessel function of the second kind, Y_ν(z).
///
/// Computes a single value of the Bessel function Y_ν(z) for complex z
/// and real order ν (any real value, including negative).
///
/// For negative ν, the DLMF 10.4.2 reflection formula is applied:
/// `Y_{-ν}(z) = sin(νπ) J_ν(z) + cos(νπ) Y_ν(z)`.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails (overflow, z = 0, etc.).
pub fn bessely<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    bessely_internal(nu, z, Scaling::Unscaled)
}

/// Modified Bessel function of the first kind, I_ν(z).
///
/// Computes a single value of I_ν(z) for complex z and real order ν
/// (any real value, including negative).
///
/// For negative ν, the DLMF 10.27.2 reflection formula is applied:
/// `I_{-ν}(z) = I_ν(z) + (2/π) sin(νπ) K_ν(z)`.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails (overflow, precision loss, etc.).
pub fn besseli<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    besseli_internal(nu, z, Scaling::Unscaled)
}

/// Modified Bessel function of the second kind, K_ν(z).
///
/// Computes a single value of K_ν(z) for complex z and real order ν
/// (any real value, including negative). K is even in ν: K_{-ν}(z) = K_ν(z).
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails (overflow, z = 0, etc.).
pub fn besselk<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    besselk_internal(nu, z, Scaling::Unscaled)
}

/// Hankel function H_ν^(m)(z), m = 1 or 2.
///
/// Computes a single value of the Hankel function for complex z and real order ν
/// (any real value, including negative).
///
/// For negative ν, the DLMF 10.4.6–7 reflection formulas are applied:
/// - `H^(1)_{-ν}(z) = exp(νπi) H^(1)_ν(z)`
/// - `H^(2)_{-ν}(z) = exp(-νπi) H^(2)_ν(z)`
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails (overflow, z = 0, etc.).
pub fn hankel<T: BesselFloat>(
    kind: HankelKind,
    nu: T,
    z: Complex<T>,
) -> Result<Complex<T>, BesselError> {
    hankel_internal(kind, nu, z, Scaling::Unscaled)
}

/// Airy function Ai(z) or its derivative Ai'(z).
///
/// Computes the Airy function for complex z. Airy functions are solutions
/// to the differential equation `w'' - z·w = 0`.
///
/// Use [`AiryDerivative::Value`] for Ai(z) or [`AiryDerivative::Derivative`] for Ai'(z).
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
pub fn airy<T: BesselFloat>(
    z: Complex<T>,
    deriv: AiryDerivative,
) -> Result<Complex<T>, BesselError> {
    let (result, _nz) = airy::zairy(z, deriv, Scaling::Unscaled)?;
    Ok(result)
}

/// Airy function Bi(z) or its derivative Bi'(z).
///
/// Computes the Airy function of the second kind for complex z.
/// Bi(z) is the solution to `w'' - z·w = 0` that grows exponentially
/// for large positive real z.
///
/// Use [`AiryDerivative::Value`] for Bi(z) or [`AiryDerivative::Derivative`] for Bi'(z).
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
pub fn biry<T: BesselFloat>(
    z: Complex<T>,
    deriv: AiryDerivative,
) -> Result<Complex<T>, BesselError> {
    airy::zbiry(z, deriv, Scaling::Unscaled)
}

// ── Scaled single-value functions ──

/// Scaled Bessel function of the first kind: `exp(-|Im(z)|) · J_ν(z)`.
///
/// The exponential factor cancels the asymptotic growth of J for large imaginary
/// arguments, keeping results in a representable floating-point range.
/// This is especially useful when |Im(z)| is large.
///
/// Supports negative ν via the same reflection formula as [`besselj`].
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
pub fn besselj_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    besselj_internal(nu, z, Scaling::Exponential)
}

/// Scaled Bessel function of the second kind: `exp(-|Im(z)|) · Y_ν(z)`.
///
/// The exponential factor cancels the asymptotic growth of Y for large imaginary
/// arguments, keeping results in a representable floating-point range.
///
/// Supports negative ν via the same reflection formula as [`bessely`].
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
pub fn bessely_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    bessely_internal(nu, z, Scaling::Exponential)
}

/// Scaled modified Bessel function of the first kind: `exp(-|Re(z)|) · I_ν(z)`.
///
/// I_ν(z) grows exponentially for large |Re(z)|, so the unscaled value can
/// easily overflow. The scaling factor `exp(-|Re(z)|)` keeps the result finite.
///
/// Supports negative ν via the same reflection formula as [`besseli`].
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
pub fn besseli_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    besseli_internal(nu, z, Scaling::Exponential)
}

/// Scaled modified Bessel function of the second kind: `exp(z) · K_ν(z)`.
///
/// K_ν(z) decays exponentially for large Re(z), so unscaled values can underflow
/// to zero. The scaling factor `exp(z)` keeps the result in a normal range.
///
/// Supports negative ν (K is even in ν: K_{-ν} = K_ν).
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Example
///
/// ```
/// use complex_bessel::besselk_scaled;
/// use num_complex::Complex;
///
/// let z = Complex::new(500.0, 0.0);
///
/// // K_0(500) would underflow to 0, but the scaled version stays finite:
/// let k_s = besselk_scaled(0.0, z).unwrap();
/// assert!(k_s.re > 0.0); // exp(500) * K_0(500) ≈ 0.0564
/// ```
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
pub fn besselk_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    besselk_internal(nu, z, Scaling::Exponential)
}

/// Scaled Hankel function H_ν^(m)(z).
///
/// - H^(1): returns `exp(-iz) · H_ν^(1)(z)`
/// - H^(2): returns `exp(iz) · H_ν^(2)(z)`
///
/// The Hankel functions grow exponentially in the complex plane;
/// the scaling factor removes this growth, preventing overflow.
///
/// Supports negative ν via the same reflection formulas as [`hankel`].
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
pub fn hankel_scaled<T: BesselFloat>(
    kind: HankelKind,
    nu: T,
    z: Complex<T>,
) -> Result<Complex<T>, BesselError> {
    hankel_internal(kind, nu, z, Scaling::Exponential)
}

/// Scaled Airy function: `exp(ζ) · Ai(z)` or `exp(ζ) · Ai'(z)`,
/// where ζ = (2/3) z√z.
///
/// Ai(z) decays super-exponentially for large positive real z.
/// The scaling factor `exp(ζ)` keeps the result representable.
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
pub fn airy_scaled<T: BesselFloat>(
    z: Complex<T>,
    deriv: AiryDerivative,
) -> Result<Complex<T>, BesselError> {
    let (result, _nz) = airy::zairy(z, deriv, Scaling::Exponential)?;
    Ok(result)
}

/// Scaled Airy function: `exp(-|Re(ζ)|) · Bi(z)` or `exp(-|Re(ζ)|) · Bi'(z)`,
/// where ζ = (2/3) z√z.
///
/// Bi(z) grows super-exponentially for large positive real z.
/// The scaling factor `exp(-|Re(ζ)|)` keeps the result representable.
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
pub fn biry_scaled<T: BesselFloat>(
    z: Complex<T>,
    deriv: AiryDerivative,
) -> Result<Complex<T>, BesselError> {
    airy::zbiry(z, deriv, Scaling::Exponential)
}

// ── Sequence functions with scaling option ──

/// Compute J_{ν+j}(z) for j = 0, 1, …, n−1 in a single call.
///
/// Returns a [`BesselResult`] containing `n` values and a [`BesselStatus`]:
/// - [`BesselStatus::Normal`] — full precision (~14 digits for f64)
/// - [`BesselStatus::ReducedPrecision`] — some precision lost (|z| or ν very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Requires ν ≥ 0. Use [`besselj`] for negative orders.
///
/// See [crate-level docs](crate#consecutive-orders) for more on sequence functions.
///
/// # Errors
///
/// Returns [`BesselError::InvalidInput`] if ν < 0 or n < 1.
pub fn besselj_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, BesselError> {
    besj::zbesj(z, nu, scaling, n)
}

/// Compute Y_{ν+j}(z) for j = 0, 1, …, n−1 in a single call.
///
/// Returns a [`BesselResult`] containing `n` values and a [`BesselStatus`]:
/// - [`BesselStatus::Normal`] — full precision (~14 digits for f64)
/// - [`BesselStatus::ReducedPrecision`] — some precision lost (|z| or ν very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Requires ν ≥ 0. Use [`bessely`] for negative orders.
///
/// See [crate-level docs](crate#consecutive-orders) for more on sequence functions.
///
/// # Errors
///
/// Returns [`BesselError::InvalidInput`] if ν < 0 or n < 1.
pub fn bessely_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, BesselError> {
    besy::zbesy(z, nu, scaling, n)
}

/// Compute I_{ν+j}(z) for j = 0, 1, …, n−1 in a single call.
///
/// Returns a [`BesselResult`] containing `n` values and a [`BesselStatus`]:
/// - [`BesselStatus::Normal`] — full precision (~14 digits for f64)
/// - [`BesselStatus::ReducedPrecision`] — some precision lost (|z| or ν very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Requires ν ≥ 0. Use [`besseli`] for negative orders.
///
/// See [crate-level docs](crate#consecutive-orders) for more on sequence functions.
///
/// # Errors
///
/// Returns [`BesselError::InvalidInput`] if ν < 0 or n < 1.
pub fn besseli_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, BesselError> {
    besi::zbesi(z, nu, scaling, n)
}

/// Compute K_{ν+j}(z) for j = 0, 1, …, n−1 in a single call.
///
/// Returns a [`BesselResult`] containing `n` values and a [`BesselStatus`]:
/// - [`BesselStatus::Normal`] — full precision (~14 digits for f64)
/// - [`BesselStatus::ReducedPrecision`] — some precision lost (|z| or ν very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Requires ν ≥ 0. Use [`besselk`] for negative orders.
///
/// See [crate-level docs](crate#consecutive-orders) for more on sequence functions.
///
/// # Example
///
/// ```
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0, 2.0);
///
/// // K_0(z), K_1(z), K_2(z) in one call
/// let result = besselk_seq(0.0, z, 3, Scaling::Unscaled).unwrap();
/// assert_eq!(result.values.len(), 3);
/// assert!(matches!(result.status, BesselStatus::Normal));
/// ```
///
/// # Errors
///
/// Returns [`BesselError::InvalidInput`] if ν < 0 or n < 1.
pub fn besselk_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, BesselError> {
    besk::zbesk(z, nu, scaling, n)
}

/// Compute H_{ν+j}^(m)(z) for j = 0, 1, …, n−1 in a single call.
///
/// Returns a [`BesselResult`] containing `n` values and a [`BesselStatus`]:
/// - [`BesselStatus::Normal`] — full precision (~14 digits for f64)
/// - [`BesselStatus::ReducedPrecision`] — some precision lost (|z| or ν very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Requires ν ≥ 0. Use [`hankel`] for negative orders.
///
/// See [crate-level docs](crate#consecutive-orders) for more on sequence functions.
///
/// # Errors
///
/// Returns [`BesselError::InvalidInput`] if ν < 0 or n < 1.
pub fn hankel_seq<T: BesselFloat>(
    kind: HankelKind,
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, BesselError> {
    besh::zbesh(z, nu, kind, scaling, n)
}
