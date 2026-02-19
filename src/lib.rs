//! Pure Rust implementation of complex Bessel functions based on Amos Algorithm 644 (ACM TOMS 644).
//!
//! Provides Bessel functions J, Y, I, K, Hankel H⁽¹⁾/H⁽²⁾, and Airy functions Ai/Bi
//! for complex arguments and real orders.
//!
//! # Features
//!
//! - **Dual precision** — all functions accept `Complex<f64>` or `Complex<f32>`
//! - **Full TOMS 644 coverage** — J, Y, I, K, H⁽¹⁾, H⁽²⁾, Ai, Bi
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
//! let h = hankel1(0.0, z).unwrap();
//!
//! // Scaled versions prevent overflow/underflow
//! let k_scaled = besselk_scaled(1.0, z).unwrap();
//!
//! // Airy functions
//! let ai = airy(z).unwrap();
//! let ai_prime = airyprime(z).unwrap();
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
//! # #[cfg(feature = "alloc")] {
//! use complex_bessel::*;
//! use num_complex::Complex;
//!
//! let z = Complex::new(1.0, 2.0);
//!
//! // K_0(z), K_1(z), K_2(z)
//! let seq = besselk_seq(0.0, z, 3, Scaling::Unscaled).unwrap();
//! assert_eq!(seq.values.len(), 3);
//! # }
//! ```
//!
//! Sequence results include a [`BesselStatus`] field:
//! - [`BesselStatus::Normal`] — full precision (~14 digits for f64)
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
//! | J, Y | exp(−\|Im(z)\|) · J(z), Y(z) |
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
//! Three tiers of feature support:
//!
//! | Feature | API | Allocator |
//! |---------|-----|-----------|
//! | (none) | 30 single-value + 8 Airy | **Not required** |
//! | `alloc` | + `_seq` variants + `BesselResult` | Required |
//! | `std` (default) | Full + `impl Error` | Required |
//!
//! ```toml
//! # Pure no_std, no allocator needed:
//! complex-bessel = { version = "0.1", default-features = false }
//!
//! # no_std with alloc (adds _seq functions):
//! complex-bessel = { version = "0.1", default-features = false, features = ["alloc"] }
//! ```

#![cfg_attr(not(feature = "std"), no_std)]

#[cfg(all(feature = "alloc", not(feature = "std")))]
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
#[cfg(feature = "alloc")]
pub use types::BesselResult;
pub use types::{BesselError, BesselStatus, Scaling};

use num_complex::Complex;
use types::{AiryDerivative, HankelKind};

// ── Helper: integer order detection ──

/// Check if `nu` is a non-negative integer. Returns `Some(n)` if so.
#[inline]
fn as_integer<T: BesselFloat>(nu: T) -> Option<i64> {
    if nu == nu.floor() {
        // Safe conversion: orders beyond i64 range are not practical
        nu.to_i64()
    } else {
        None
    }
}

// ── Internal: compute with given scaling for negative order support ──

#[inline]
fn besselj_internal<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
) -> Result<Complex<T>, BesselError> {
    let zero = T::zero();
    let czero = Complex::new(zero, zero);
    if nu >= zero {
        let mut y = [czero];
        besj::zbesj(z, nu, scaling, &mut y)?;
        return Ok(y[0]);
    }

    // Negative order: J_{-ν}(z) = cos(νπ)*J_ν(z) - sin(νπ)*Y_ν(z) (DLMF 10.4.1)
    let abs_nu = nu.abs();

    // Integer shortcut: J_{-n}(z) = (-1)^n * J_n(z)
    if let Some(n) = as_integer(abs_nu) {
        let mut y = [czero];
        besj::zbesj(z, abs_nu, scaling, &mut y)?;
        let sign = if n % 2 == 0 { T::one() } else { -T::one() };
        return Ok(y[0] * sign);
    }

    // General case: need both J and Y at positive |ν|
    let cos_nu_pi = utils::cospi(abs_nu);
    let sin_nu_pi = utils::sinpi(abs_nu);

    let mut j_buf = [czero];
    let mut y_buf = [czero];
    besj::zbesj(z, abs_nu, scaling, &mut j_buf)?;
    besy::zbesy(z, abs_nu, scaling, &mut y_buf)?;

    Ok(j_buf[0] * cos_nu_pi - y_buf[0] * sin_nu_pi)
}

#[inline]
fn bessely_internal<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
) -> Result<Complex<T>, BesselError> {
    let zero = T::zero();
    let czero = Complex::new(zero, zero);
    if nu >= zero {
        let mut y = [czero];
        besy::zbesy(z, nu, scaling, &mut y)?;
        return Ok(y[0]);
    }

    // Negative order: Y_{-ν}(z) = sin(νπ)*J_ν(z) + cos(νπ)*Y_ν(z) (DLMF 10.4.2)
    let abs_nu = nu.abs();

    // Integer shortcut: Y_{-n}(z) = (-1)^n * Y_n(z)
    if let Some(n) = as_integer(abs_nu) {
        let mut y = [czero];
        besy::zbesy(z, abs_nu, scaling, &mut y)?;
        let sign = if n % 2 == 0 { T::one() } else { -T::one() };
        return Ok(y[0] * sign);
    }

    // General case: need both J and Y at positive |ν|
    let cos_nu_pi = utils::cospi(abs_nu);
    let sin_nu_pi = utils::sinpi(abs_nu);

    let mut j_buf = [czero];
    let mut y_buf = [czero];
    besj::zbesj(z, abs_nu, scaling, &mut j_buf)?;
    besy::zbesy(z, abs_nu, scaling, &mut y_buf)?;

    Ok(j_buf[0] * sin_nu_pi + y_buf[0] * cos_nu_pi)
}

#[inline]
fn besseli_internal<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
) -> Result<Complex<T>, BesselError> {
    let zero = T::zero();
    let czero = Complex::new(zero, zero);
    if nu >= zero {
        let mut y = [czero];
        besi::zbesi(z, nu, scaling, &mut y)?;
        return Ok(y[0]);
    }

    // Negative order: I_{-ν}(z) = I_ν(z) + (2/π)*sin(νπ)*K_ν(z) (DLMF 10.27.2)
    let abs_nu = nu.abs();

    // Integer shortcut: I_{-n}(z) = I_n(z)
    if as_integer(abs_nu).is_some() {
        let mut y = [czero];
        besi::zbesi(z, abs_nu, scaling, &mut y)?;
        return Ok(y[0]);
    }

    // General case: need both I and K at positive |ν|
    let pi = T::from_f64(core::f64::consts::PI);
    let two = T::from_f64(2.0);
    let sin_nu_pi = utils::sinpi(abs_nu);

    let mut i_buf = [czero];
    let mut k_buf = [czero];
    besi::zbesi(z, abs_nu, scaling, &mut i_buf)?;
    besk::zbesk(z, abs_nu, scaling, &mut k_buf)?;

    let i_val = i_buf[0];
    let mut k_val = k_buf[0];

    // Scaling correction: I_scaled = exp(-|Re(z)|)*I, K_scaled = exp(z)*K.
    // To combine them we must convert K to the same scaling as I:
    //   factor = exp(-|Re(z)|) / exp(z) = exp(-i*Im(z)) * [exp(-2*Re(z)) if Re(z)>0]
    if scaling == Scaling::Exponential {
        let (sin_a, cos_a) = (-z.im).sin_cos();
        k_val = k_val * Complex::new(cos_a, sin_a);
        if z.re > T::zero() {
            let scale = (-two * z.re).exp();
            k_val = k_val * scale;
        }
    }

    Ok(i_val + k_val * (two / pi * sin_nu_pi))
}

#[inline]
fn besselk_internal<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
) -> Result<Complex<T>, BesselError> {
    // K_{-ν}(z) = K_ν(z) (DLMF 10.27.3) — K is even in ν
    let abs_nu = nu.abs();
    let zero = T::zero();
    let mut y = [Complex::new(zero, zero)];
    besk::zbesk(z, abs_nu, scaling, &mut y)?;
    Ok(y[0])
}

#[inline]
fn hankel_internal<T: BesselFloat>(
    kind: HankelKind,
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
) -> Result<Complex<T>, BesselError> {
    let zero = T::zero();
    if nu >= zero {
        let mut y = [Complex::new(zero, zero)];
        besh::zbesh(z, nu, kind, scaling, &mut y)?;
        return Ok(y[0]);
    }

    // Negative order (DLMF 10.4.6, 10.4.7):
    //   H^(1)_{-ν}(z) = exp(νπi) * H^(1)_ν(z)
    //   H^(2)_{-ν}(z) = exp(-νπi) * H^(2)_ν(z)
    let abs_nu = nu.abs();
    let mut y = [Complex::new(zero, zero)];
    besh::zbesh(z, abs_nu, kind, scaling, &mut y)?;
    let h_val = y[0];

    let cos_nu_pi = utils::cospi(abs_nu);
    let sin_nu_pi = utils::sinpi(abs_nu);

    let rotation = match kind {
        // exp(νπi) = cos(νπ) + i*sin(νπ)
        HankelKind::First => Complex::new(cos_nu_pi, sin_nu_pi),
        // exp(-νπi) = cos(νπ) - i*sin(νπ)
        HankelKind::Second => Complex::new(cos_nu_pi, -sin_nu_pi),
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
#[inline]
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
#[inline]
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
#[inline]
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
#[inline]
pub fn besselk<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    besselk_internal(nu, z, Scaling::Unscaled)
}

/// Hankel function of the first kind, H_ν^(1)(z).
///
/// Computes a single value of H_ν^(1)(z) for complex z and real order ν
/// (any real value, including negative).
///
/// For negative ν, the DLMF 10.4.6 reflection formula is applied:
/// `H^(1)_{-ν}(z) = exp(νπi) H^(1)_ν(z)`.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails (overflow, z = 0, etc.).
#[inline]
pub fn hankel1<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    hankel_internal(HankelKind::First, nu, z, Scaling::Unscaled)
}

/// Hankel function of the second kind, H_ν^(2)(z).
///
/// Computes a single value of H_ν^(2)(z) for complex z and real order ν
/// (any real value, including negative).
///
/// For negative ν, the DLMF 10.4.7 reflection formula is applied:
/// `H^(2)_{-ν}(z) = exp(-νπi) H^(2)_ν(z)`.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails (overflow, z = 0, etc.).
#[inline]
pub fn hankel2<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    hankel_internal(HankelKind::Second, nu, z, Scaling::Unscaled)
}

/// Airy function Ai(z).
///
/// Computes the Airy function of the first kind for complex z.
/// Ai(z) is a solution to the differential equation `w'' - z·w = 0`
/// that decays exponentially for large positive real z.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
#[inline]
pub fn airy<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, BesselError> {
    let (result, _nz) = airy::zairy(z, AiryDerivative::Value, Scaling::Unscaled)?;
    Ok(result)
}

/// Derivative of the Airy function, Ai'(z).
///
/// Computes the derivative of the Airy function of the first kind for complex z.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
#[inline]
pub fn airyprime<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, BesselError> {
    let (result, _nz) = airy::zairy(z, AiryDerivative::Derivative, Scaling::Unscaled)?;
    Ok(result)
}

/// Airy function of the second kind, Bi(z).
///
/// Computes the Airy function Bi(z) for complex z.
/// Bi(z) is the solution to `w'' - z·w = 0` that grows exponentially
/// for large positive real z.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
#[inline]
pub fn biry<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, BesselError> {
    airy::zbiry(z, AiryDerivative::Value, Scaling::Unscaled)
}

/// Derivative of the Airy function of the second kind, Bi'(z).
///
/// Computes the derivative of Bi(z) for complex z.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
#[inline]
pub fn biryprime<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, BesselError> {
    airy::zbiry(z, AiryDerivative::Derivative, Scaling::Unscaled)
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
#[inline]
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
#[inline]
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
#[inline]
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
#[inline]
pub fn besselk_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    besselk_internal(nu, z, Scaling::Exponential)
}

/// Scaled Hankel function of the first kind: `exp(-iz) · H_ν^(1)(z)`.
///
/// H^(1) grows exponentially in the lower half-plane;
/// the scaling factor removes this growth, preventing overflow.
///
/// Supports negative ν via the same reflection formula as [`hankel1`].
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
#[inline]
pub fn hankel1_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    hankel_internal(HankelKind::First, nu, z, Scaling::Exponential)
}

/// Scaled Hankel function of the second kind: `exp(iz) · H_ν^(2)(z)`.
///
/// H^(2) grows exponentially in the upper half-plane;
/// the scaling factor removes this growth, preventing overflow.
///
/// Supports negative ν via the same reflection formula as [`hankel2`].
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
#[inline]
pub fn hankel2_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    hankel_internal(HankelKind::Second, nu, z, Scaling::Exponential)
}

/// Scaled Airy function: `exp(ζ) · Ai(z)`, where ζ = (2/3) z√z.
///
/// Ai(z) decays super-exponentially for large positive real z.
/// The scaling factor `exp(ζ)` keeps the result representable.
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
#[inline]
pub fn airy_scaled<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, BesselError> {
    let (result, _nz) = airy::zairy(z, AiryDerivative::Value, Scaling::Exponential)?;
    Ok(result)
}

/// Scaled derivative of the Airy function: `exp(ζ) · Ai'(z)`, where ζ = (2/3) z√z.
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
#[inline]
pub fn airyprime_scaled<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, BesselError> {
    let (result, _nz) = airy::zairy(z, AiryDerivative::Derivative, Scaling::Exponential)?;
    Ok(result)
}

/// Scaled Airy function of the second kind: `exp(-|Re(ζ)|) · Bi(z)`,
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
#[inline]
pub fn biry_scaled<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, BesselError> {
    airy::zbiry(z, AiryDerivative::Value, Scaling::Exponential)
}

/// Scaled derivative of the Airy function of the second kind:
/// `exp(-|Re(ζ)|) · Bi'(z)`, where ζ = (2/3) z√z.
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Errors
///
/// Returns [`BesselError`] if the computation fails.
#[inline]
pub fn biryprime_scaled<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, BesselError> {
    airy::zbiry(z, AiryDerivative::Derivative, Scaling::Exponential)
}

// ── Sequence functions with scaling option (require alloc) ──

#[cfg(feature = "alloc")]
extern crate alloc as alloc_crate;

#[cfg(feature = "alloc")]
fn seq_helper<T: BesselFloat>(
    n: usize,
    f: impl FnOnce(&mut [Complex<T>]) -> Result<(usize, BesselStatus), BesselError>,
) -> Result<BesselResult<T>, BesselError> {
    let zero = T::zero();
    let mut values = alloc_crate::vec![Complex::new(zero, zero); n];
    let (underflow_count, status) = f(&mut values)?;
    Ok(BesselResult {
        values,
        underflow_count,
        status,
    })
}

#[cfg(feature = "alloc")]
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
    seq_helper(n, |y| besj::zbesj(z, nu, scaling, y))
}

#[cfg(feature = "alloc")]
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
    seq_helper(n, |y| besy::zbesy(z, nu, scaling, y))
}

#[cfg(feature = "alloc")]
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
    seq_helper(n, |y| besi::zbesi(z, nu, scaling, y))
}

#[cfg(feature = "alloc")]
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
    seq_helper(n, |y| besk::zbesk(z, nu, scaling, y))
}

#[cfg(feature = "alloc")]
/// Compute H_{ν+j}^(1)(z) for j = 0, 1, …, n−1 in a single call.
///
/// Returns a [`BesselResult`] containing `n` values and a [`BesselStatus`]:
/// - [`BesselStatus::Normal`] — full precision (~14 digits for f64)
/// - [`BesselStatus::ReducedPrecision`] — some precision lost (|z| or ν very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Requires ν ≥ 0. Use [`hankel1`] for negative orders.
///
/// See [crate-level docs](crate#consecutive-orders) for more on sequence functions.
///
/// # Errors
///
/// Returns [`BesselError::InvalidInput`] if ν < 0 or n < 1.
pub fn hankel1_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, BesselError> {
    seq_helper(n, |y| besh::zbesh(z, nu, HankelKind::First, scaling, y))
}

#[cfg(feature = "alloc")]
/// Compute H_{ν+j}^(2)(z) for j = 0, 1, …, n−1 in a single call.
///
/// Returns a [`BesselResult`] containing `n` values and a [`BesselStatus`]:
/// - [`BesselStatus::Normal`] — full precision (~14 digits for f64)
/// - [`BesselStatus::ReducedPrecision`] — some precision lost (|z| or ν very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Requires ν ≥ 0. Use [`hankel2`] for negative orders.
///
/// See [crate-level docs](crate#consecutive-orders) for more on sequence functions.
///
/// # Errors
///
/// Returns [`BesselError::InvalidInput`] if ν < 0 or n < 1.
pub fn hankel2_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, BesselError> {
    seq_helper(n, |y| besh::zbesh(z, nu, HankelKind::Second, scaling, y))
}
