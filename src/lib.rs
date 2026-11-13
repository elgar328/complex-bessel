//! Pure Rust implementation of complex Bessel functions based on Amos Algorithm 644 (ACM TOMS 644).
//!
//! # Features
//!
//! - **f32 & f64** — all functions accept `Complex<f64>` or `Complex<f32>`
//! - **Complete function set** — J, Y, I, K, H<sup>(1)</sup>, H<sup>(2)</sup>, Ai, Bi
//! - **Consecutive orders** — `_seq` variants return ν, ν+1, …, ν+n−1 in one call;
//!   `_seq_into` variants write them into a caller-provided slice
//! - **Exponential scaling** — `_scaled` variants prevent overflow/underflow
//! - **Negative orders** — supports ν < 0 via DLMF reflection formulas (not in Amos)
//! - **`no_std` support** — everything but the `Vec`-returning `_seq` variants (which need `std`) works without an allocator
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
//! # Functions
//!
//! | Function | Description | Scaled variant returns |
//! |----------|-------------|------------------------|
//! | [`besselj`]`(ν, z)`<br>[`besselj_scaled`]`(ν, z)`<br>[`besselj_seq`]`(ν, z, n, scaling)`<br>[`besselj_seq_into`]`(ν, z, scaling, out)` | J<sub>ν</sub>(z), Bessel first kind | exp(−\|Im(z)\|) · J<sub>ν</sub>(z) |
//! | [`bessely`]`(ν, z)`<br>[`bessely_scaled`]`(ν, z)`<br>[`bessely_seq`]`(ν, z, n, scaling)`<br>[`bessely_seq_into`]`(ν, z, scaling, out)` | Y<sub>ν</sub>(z), Bessel second kind | exp(−\|Im(z)\|) · Y<sub>ν</sub>(z) |
//! | [`besseli`]`(ν, z)`<br>[`besseli_scaled`]`(ν, z)`<br>[`besseli_seq`]`(ν, z, n, scaling)`<br>[`besseli_seq_into`]`(ν, z, scaling, out)` | I<sub>ν</sub>(z), modified first kind | exp(−\|Re(z)\|) · I<sub>ν</sub>(z) |
//! | [`besselk`]`(ν, z)`<br>[`besselk_scaled`]`(ν, z)`<br>[`besselk_seq`]`(ν, z, n, scaling)`<br>[`besselk_seq_into`]`(ν, z, scaling, out)` | K<sub>ν</sub>(z), modified second kind | exp(z) · K<sub>ν</sub>(z) |
//! | [`hankel1`]`(ν, z)`<br>[`hankel1_scaled`]`(ν, z)`<br>[`hankel1_seq`]`(ν, z, n, scaling)`<br>[`hankel1_seq_into`]`(ν, z, scaling, out)` | H<sub>ν</sub><sup>(1)</sup>(z), Hankel first kind | exp(−iz) · H<sub>ν</sub><sup>(1)</sup>(z) |
//! | [`hankel2`]`(ν, z)`<br>[`hankel2_scaled`]`(ν, z)`<br>[`hankel2_seq`]`(ν, z, n, scaling)`<br>[`hankel2_seq_into`]`(ν, z, scaling, out)` | H<sub>ν</sub><sup>(2)</sup>(z), Hankel second kind | exp(iz) · H<sub>ν</sub><sup>(2)</sup>(z) |
//! | [`airy`]`(z)`<br>[`airy_scaled`]`(z)`<br>[`airy_raw`]`(z, scaling)` | Ai(z), Airy first kind | exp(ζ) · Ai(z) |
//! | [`airyprime`]`(z)`<br>[`airyprime_scaled`]`(z)`<br>[`airyprime_raw`]`(z, scaling)` | Ai′(z), derivative of Airy first kind | exp(ζ) · Ai′(z) |
//! | [`biry`]`(z)`<br>[`biry_scaled`]`(z)`<br>[`biry_raw`]`(z, scaling)` | Bi(z), Airy second kind | exp(−\|Re(ζ)\|) · Bi(z) |
//! | [`biryprime`]`(z)`<br>[`biryprime_scaled`]`(z)`<br>[`biryprime_raw`]`(z, scaling)` | Bi′(z), derivative of Airy second kind | exp(−\|Re(ζ)\|) · Bi′(z) |
//!
//! where ζ = (2/3) z√z.
//!
//! ## Exponential scaling
//!
//! Each `_scaled` variant multiplies the result by an exponential factor that
//! cancels the asymptotic growth (or decay), keeping values in a representable
//! floating-point range. The exact factor is listed in the "Scaled variant
//! returns" column above. Use scaled functions whenever |z| or |Im(z)| is
//! large enough that the unscaled result would overflow or underflow.
//!
//! # Negative orders
//!
//! All functions accept any real order, including negative values.
//! DLMF reflection formulas are applied automatically:
//!
//! - **J**: J<sub>−ν</sub>(z) = cos(νπ) J<sub>ν</sub>(z) − sin(νπ) Y<sub>ν</sub>(z) (DLMF 10.2.3)
//! - **Y**: Y<sub>−ν</sub>(z) = sin(νπ) J<sub>ν</sub>(z) + cos(νπ) Y<sub>ν</sub>(z) (DLMF 10.2.3)
//! - **I**: I<sub>−ν</sub>(z) = I<sub>ν</sub>(z) + (2/π) sin(νπ) K<sub>ν</sub>(z) (DLMF 10.27.2)
//! - **K**: K<sub>−ν</sub>(z) = K<sub>ν</sub>(z) (DLMF 10.27.3)
//! - **H<sup>(1)</sup>**: H<sup>(1)</sup><sub>−ν</sub>(z) = exp(νπi) H<sup>(1)</sup><sub>ν</sub>(z) (DLMF 10.4.6)
//! - **H<sup>(2)</sup>**: H<sup>(2)</sup><sub>−ν</sub>(z) = exp(−νπi) H<sup>(2)</sup><sub>ν</sub>(z) (DLMF 10.4.6)
//!
//! For integer orders, simplified identities are used (e.g., J<sub>−n</sub>(z) = (−1)<sup>n</sup> J<sub>n</sub>(z)).
//!
//! # Function variants
//!
//! ## Consecutive orders
//!
//! The `_seq` variants ([`besselj_seq`], [`besselk_seq`], …) compute values at
//! consecutive orders ν, ν+1, …, ν+n−1 in a single call and return
//! a [`BesselResult`] that includes an [`Accuracy`] field.
//! The `_seq_into` variants ([`besselj_seq_into`], [`besselk_seq_into`], …)
//! write the sequence into a caller-provided slice and return a [`SeqInfo`]
//! with [`Accuracy`]; they never allocate.
//! The `_raw` Airy variants ([`airy_raw`], [`biry_raw`], …) similarly return
//! an [`AiryResult`] with [`Accuracy`].
//! All other functions return only the computed value without [`Accuracy`].
//!
//! **[`Accuracy`]:**
//!
//! | Status | Meaning |
//! |--------|---------|
//! | [`Normal`](Accuracy::Normal) | Full machine precision |
//! | [`Reduced`](Accuracy::Reduced) | More than half of significant digits may be lost |
//!
//! [`Reduced`](Accuracy::Reduced) is extremely rare in practice. SciPy's Bessel wrappers also silently
//! discard the equivalent Amos IERR=3 flag by default.
//!
//! To check accuracy status, use a `_seq`/`_seq_into` function (Bessel) or a
//! `_raw` function (Airy):
//!
//! ```
//! use complex_bessel::*;
//! use num_complex::Complex;
//!
//! let z = Complex::new(1.0, 2.0);
//! let mut out = [Complex::new(0.0, 0.0); 1];
//! let info = besselk_seq_into(0.0, z, Scaling::Unscaled, &mut out).unwrap();
//! assert!(matches!(info.status, Accuracy::Normal));
//!
//! let result = airy_raw(z, Scaling::Unscaled).unwrap();
//! assert!(matches!(result.status, Accuracy::Normal));
//! ```
//!
//! # Error handling
//!
//! All functions return `Result<_, Error>`. The four error variants are:
//!
//! | Variant | Cause |
//! |---------|-------|
//! | [`InvalidInput`](Error::InvalidInput) | non-finite z or ν, z = 0 for K/Y/H, n < 1 (empty `out`) |
//! | [`Overflow`](Error::Overflow) | \|z\| or \|ν\| too large (or \|z\| too small) for finite result |
//! | [`TotalPrecisionLoss`](Error::TotalPrecisionLoss) | \|z\| or \|ν\| too large for meaningful computation |
//! | [`ConvergenceFailure`](Error::ConvergenceFailure) | Internal algorithm did not converge |
//!
//! # `no_std` support
//!
//! | Cargo features | Available API |
//! |---------------|---------------|
//! | `default-features = false` | 24 single-value functions + 6 `_seq_into` variants |
//! | `features = ["std"]` (default) | full API (+ 6 `_seq` variants), with platform-native math |
//!
//! The 24 single-value functions include 12 Bessel (J/Y/I/K/H<sup>(1)</sup>/H<sup>(2)</sup> × unscaled/scaled)
//! and 12 Airy (Ai/Ai'/Bi/Bi' × unscaled/scaled/raw). No allocator is needed;
//! without `std`, math functions come from `libm`.
//!
//! ```toml
//! complex-bessel = { version = "0.4", default-features = false }
//! ```

#![warn(missing_docs)]
// Without `std`, the docs still link to the `std`-only `_seq` functions.
#![cfg_attr(not(feature = "std"), allow(rustdoc::broken_intra_doc_links))]
#![no_std]

#[cfg(feature = "std")]
extern crate alloc;
#[cfg(feature = "std")]
extern crate std;

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
#[cfg(feature = "std")]
pub use types::BesselResult;
pub use types::{Accuracy, AiryResult, Error, Scaling, SeqInfo};

use num_complex::Complex;
use types::{AiryDerivative, HankelKind};
use utils::{Workspace, with_workspace};

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

// ── Element-wise reflection helpers (shared by single-value and _seq) ──

/// J_{-ν}(z) = cos(νπ)·J_ν(z) − sin(νπ)·Y_ν(z)  (DLMF 10.2.3)
#[inline]
fn reflect_j_element<T: BesselFloat>(order: T, j: Complex<T>, y: Complex<T>) -> Complex<T> {
    j * utils::cospi(order) - y * utils::sinpi(order)
}

/// Y_{-ν}(z) = sin(νπ)·J_ν(z) + cos(νπ)·Y_ν(z)  (DLMF 10.2.3)
#[inline]
fn reflect_y_element<T: BesselFloat>(order: T, j: Complex<T>, y: Complex<T>) -> Complex<T> {
    j * utils::sinpi(order) + y * utils::cospi(order)
}

/// I_{-ν}(z) = I_ν(z) + (2/π)·sin(νπ)·K_ν(z)  (DLMF 10.27.2)
/// `k_val` must already have scaling correction applied by the caller.
#[inline]
fn reflect_i_element<T: BesselFloat>(order: T, i_val: Complex<T>, k_val: Complex<T>) -> Complex<T> {
    let pi = T::from_f64(core::f64::consts::PI);
    let two = T::from_f64(2.0);
    utils::mul_add_scalar(k_val, two / pi * utils::sinpi(order), i_val)
}

/// H^(m)_{-ν}(z) = exp(±νπi)·H^(m)_ν(z)  (DLMF 10.4.6/7)
#[inline]
fn reflect_h_element<T: BesselFloat>(order: T, kind: HankelKind, h: Complex<T>) -> Complex<T> {
    let cos_nu_pi = utils::cospi(order);
    let sin_nu_pi = utils::sinpi(order);
    let rotation = match kind {
        HankelKind::First => Complex::new(cos_nu_pi, sin_nu_pi),
        HankelKind::Second => Complex::new(cos_nu_pi, -sin_nu_pi),
    };
    h * rotation
}

/// Convert K scaling to I scaling for the reflection formula.
/// I_scaled = exp(-|Re(z)|)·I, K_scaled = exp(z)·K.
/// factor = exp(-|Re(z)|) / exp(z) = exp(-i·Im(z)) · [exp(-2·Re(z)) if Re(z)>0]
#[inline]
fn k_to_i_scaling_correction<T: BesselFloat>(z: Complex<T>, k_val: Complex<T>) -> Complex<T> {
    let two = T::from_f64(2.0);
    let (sin_a, cos_a) = (-z.im).sin_cos();
    let mut result = k_val * Complex::new(cos_a, sin_a);
    if z.re > T::zero() {
        result = result * (-two * z.re).exp();
    }
    result
}

/// (-1)^n sign factor for integer order reflection.
#[inline]
fn integer_sign<T: BesselFloat>(n: i64) -> T {
    if n % 2 == 0 { T::one() } else { -T::one() }
}

// ── Internal: compute with given scaling for negative order support ──

#[inline]
fn besselj_internal<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
) -> Result<Complex<T>, Error> {
    let zero = T::zero();
    let czero = Complex::new(zero, zero);
    if nu >= zero {
        let mut y = [czero];
        besj::zbesj(z, nu, scaling, &mut y)?;
        return Ok(y[0]);
    }

    let abs_nu = nu.abs();

    // Integer shortcut: J_{-n}(z) = (-1)^n * J_n(z)
    if let Some(n) = as_integer(abs_nu) {
        let mut y = [czero];
        besj::zbesj(z, abs_nu, scaling, &mut y)?;
        return Ok(y[0] * integer_sign::<T>(n));
    }

    // General case: need both J and Y at positive |ν|
    let mut j_buf = [czero];
    let mut y_buf = [czero];
    besj::zbesj(z, abs_nu, scaling, &mut j_buf)?;
    besy::zbesy(z, abs_nu, scaling, &mut y_buf)?;

    Ok(reflect_j_element(abs_nu, j_buf[0], y_buf[0]))
}

#[inline]
fn bessely_internal<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
) -> Result<Complex<T>, Error> {
    let zero = T::zero();
    let czero = Complex::new(zero, zero);
    if nu >= zero {
        let mut y = [czero];
        besy::zbesy(z, nu, scaling, &mut y)?;
        return Ok(y[0]);
    }

    let abs_nu = nu.abs();

    // Integer shortcut: Y_{-n}(z) = (-1)^n * Y_n(z)
    if let Some(n) = as_integer(abs_nu) {
        let mut y = [czero];
        besy::zbesy(z, abs_nu, scaling, &mut y)?;
        return Ok(y[0] * integer_sign::<T>(n));
    }

    // General case: need both J and Y at positive |ν|
    let mut j_buf = [czero];
    let mut y_buf = [czero];
    besj::zbesj(z, abs_nu, scaling, &mut j_buf)?;
    besy::zbesy(z, abs_nu, scaling, &mut y_buf)?;

    Ok(reflect_y_element(abs_nu, j_buf[0], y_buf[0]))
}

#[inline]
fn besseli_internal<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
) -> Result<Complex<T>, Error> {
    let zero = T::zero();
    let czero = Complex::new(zero, zero);
    if nu >= zero {
        let mut y = [czero];
        besi::zbesi(z, nu, scaling, &mut y)?;
        return Ok(y[0]);
    }

    let abs_nu = nu.abs();

    // Integer shortcut: I_{-n}(z) = I_n(z)
    if as_integer(abs_nu).is_some() {
        let mut y = [czero];
        besi::zbesi(z, abs_nu, scaling, &mut y)?;
        return Ok(y[0]);
    }

    // General case: need both I and K at positive |ν|
    let mut i_buf = [czero];
    let mut k_buf = [czero];
    besi::zbesi(z, abs_nu, scaling, &mut i_buf)?;
    besk::zbesk(z, abs_nu, scaling, &mut k_buf)?;

    let mut k_val = k_buf[0];
    if scaling == Scaling::Exponential {
        k_val = k_to_i_scaling_correction(z, k_val);
    }

    Ok(reflect_i_element(abs_nu, i_buf[0], k_val))
}

#[inline]
fn besselk_internal<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
) -> Result<Complex<T>, Error> {
    // K_{-ν}(z) = K_ν(z) (DLMF 10.27.3)
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
) -> Result<Complex<T>, Error> {
    let zero = T::zero();
    if nu >= zero {
        let mut y = [Complex::new(zero, zero)];
        besh::zbesh(z, nu, kind, scaling, &mut y)?;
        return Ok(y[0]);
    }

    let abs_nu = nu.abs();
    let mut y = [Complex::new(zero, zero)];
    besh::zbesh(z, abs_nu, kind, scaling, &mut y)?;

    Ok(reflect_h_element(abs_nu, kind, y[0]))
}

// ── Single-value convenience functions ──

/// Bessel function of the first kind, J<sub>ν</sub>(z).
///
/// Computes a single value of J_ν(z) for complex z and real order ν
/// (any real value, including negative).
///
/// For negative ν, the DLMF 10.2.3 reflection formula is applied:
/// `J_{-ν}(z) = cos(νπ) J_ν(z) - sin(νπ) Y_ν(z)`.
/// For negative integer n, `J_{−n}(z) = (−1)^n J_n(z)`.
///
/// # Example
///
/// ```
/// use complex_bessel::besselj;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
/// let j0 = besselj(0.0, z).unwrap();
/// assert!((j0.re - 0.7652).abs() < 1e-3); // J_0(1) ≈ 0.7652
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z or ν is not finite (NaN or infinite), or if z = 0 and ν is a negative non-integer.
/// - [`Error::Overflow`] if |z| or |ν| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn besselj<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    besselj_internal(nu, z, Scaling::Unscaled)
}

/// Bessel function of the second kind, Y<sub>ν</sub>(z).
///
/// Computes a single value of Y_ν(z) for complex z and real order ν
/// (any real value, including negative).
///
/// For negative ν, the DLMF 10.2.3 reflection formula is applied:
/// `Y_{-ν}(z) = sin(νπ) J_ν(z) + cos(νπ) Y_ν(z)`.
/// For negative integer n, `Y_{−n}(z) = (−1)^n Y_n(z)`.
///
/// # Example
///
/// ```
/// use complex_bessel::bessely;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
/// let y0 = bessely(0.0, z).unwrap();
/// assert!((y0.re - 0.0883).abs() < 1e-3); // Y_0(1) ≈ 0.0883
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z or ν is not finite (NaN or infinite), or if z = 0.
/// - [`Error::Overflow`] if |z| is too small or too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn bessely<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    bessely_internal(nu, z, Scaling::Unscaled)
}

/// Modified Bessel function of the first kind, I<sub>ν</sub>(z).
///
/// Computes a single value of I_ν(z) for complex z and real order ν
/// (any real value, including negative).
///
/// For negative ν, the DLMF 10.27.2 reflection formula is applied:
/// `I_{-ν}(z) = I_ν(z) + (2/π) sin(νπ) K_ν(z)`.
/// For negative integer n, `I_{−n}(z) = I_n(z)`.
///
/// # Example
///
/// ```
/// use complex_bessel::besseli;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
/// let i0 = besseli(0.0, z).unwrap();
/// assert!((i0.re - 1.2661).abs() < 1e-3); // I_0(1) ≈ 1.2661
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z or ν is not finite (NaN or infinite), or if z = 0 and ν is a negative non-integer.
/// - [`Error::Overflow`] if |z| or |ν| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn besseli<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    besseli_internal(nu, z, Scaling::Unscaled)
}

/// Modified Bessel function of the second kind, K<sub>ν</sub>(z).
///
/// Computes a single value of K_ν(z) for complex z and real order ν
/// (any real value, including negative).
///
/// For negative ν, K_{-ν}(z) = K_ν(z) (DLMF 10.27.3).
///
/// # Example
///
/// ```
/// use complex_bessel::besselk;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
/// let k0 = besselk(0.0, z).unwrap();
/// assert!((k0.re - 0.4211).abs() < 1e-3); // K_0(1) ≈ 0.4211
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z or ν is not finite (NaN or infinite), or if z = 0.
/// - [`Error::Overflow`] if |z| is too small or too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn besselk<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    besselk_internal(nu, z, Scaling::Unscaled)
}

/// Hankel function of the first kind, H<sub>ν</sub><sup>(1)</sup>(z).
///
/// Computes a single value of H_ν^(1)(z) for complex z and real order ν
/// (any real value, including negative).
///
/// For negative ν, the DLMF 10.4.6 reflection formula is applied:
/// `H^(1)_{-ν}(z) = exp(νπi) H^(1)_ν(z)`.
/// For negative integer n, `H^(1)_{−n}(z) = (−1)^n H^(1)_n(z)`.
///
/// # Example
///
/// ```
/// use complex_bessel::hankel1;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
/// let h = hankel1(0.0, z).unwrap();
/// assert!((h.re - 0.7652).abs() < 1e-3); // Re = J_0(1) ≈ 0.7652
/// assert!((h.im - 0.0883).abs() < 1e-3); // Im = Y_0(1) ≈ 0.0883
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z or ν is not finite (NaN or infinite), or if z = 0.
/// - [`Error::Overflow`] if |z| is too small or too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn hankel1<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    hankel_internal(HankelKind::First, nu, z, Scaling::Unscaled)
}

/// Hankel function of the second kind, H<sub>ν</sub><sup>(2)</sup>(z).
///
/// Computes a single value of H_ν^(2)(z) for complex z and real order ν
/// (any real value, including negative).
///
/// For negative ν, the DLMF 10.4.6 reflection formula is applied:
/// `H^(2)_{-ν}(z) = exp(-νπi) H^(2)_ν(z)`.
/// For negative integer n, `H^(2)_{−n}(z) = (−1)^n H^(2)_n(z)`.
///
/// # Example
///
/// ```
/// use complex_bessel::hankel2;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
/// let h = hankel2(0.0, z).unwrap();
/// assert!((h.re - 0.7652).abs() < 1e-3); // Re = J_0(1) ≈ 0.7652
/// assert!((h.im + 0.0883).abs() < 1e-3); // Im = -Y_0(1) ≈ -0.0883
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z or ν is not finite (NaN or infinite), or if z = 0.
/// - [`Error::Overflow`] if |z| is too small or too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn hankel2<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    hankel_internal(HankelKind::Second, nu, z, Scaling::Unscaled)
}

/// Airy function Ai(z).
///
/// Computes the Airy function Ai(z) for complex z.
/// Ai(z) is a solution to the differential equation w′′ − z·w = 0
/// that decays super-exponentially for large positive real z.
///
/// # Example
///
/// ```
/// use complex_bessel::airy;
/// use num_complex::Complex;
///
/// let z = Complex::new(0.0_f64, 0.0);
/// let ai = airy(z).unwrap();
/// assert!((ai.re - 0.3550).abs() < 1e-3); // Ai(0) ≈ 0.3550
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z is not finite (NaN or infinite).
/// - [`Error::Overflow`] if |z| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| is too large for any significant digits (roughly > 10⁶ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn airy<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, Error> {
    let (result, _nz, _status) = airy::zairy(z, AiryDerivative::Value, Scaling::Unscaled)?;
    Ok(result)
}

/// Derivative of the Airy function, Ai'(z).
///
/// Computes the derivative of Ai(z) for complex z.
/// Ai'(z) decays super-exponentially for large positive real z, just as Ai(z) does.
/// Satisfies the differential equation Ai′′(z) = z · Ai(z).
///
/// # Example
///
/// ```
/// use complex_bessel::airyprime;
/// use num_complex::Complex;
///
/// let z = Complex::new(0.0_f64, 0.0);
/// let ai_prime = airyprime(z).unwrap();
/// assert!((ai_prime.re + 0.2588).abs() < 1e-3); // Ai'(0) ≈ -0.2588
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z is not finite (NaN or infinite).
/// - [`Error::Overflow`] if |z| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| is too large for any significant digits (roughly > 10⁶ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn airyprime<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, Error> {
    let (result, _nz, _status) = airy::zairy(z, AiryDerivative::Derivative, Scaling::Unscaled)?;
    Ok(result)
}

/// Airy function of the second kind, Bi(z).
///
/// Computes the Airy function Bi(z) for complex z.
/// Bi(z) is the solution to w′′ − z·w = 0 that grows super-exponentially
/// for large positive real z.
///
/// # Example
///
/// ```
/// use complex_bessel::biry;
/// use num_complex::Complex;
///
/// let z = Complex::new(0.0_f64, 0.0);
/// let bi = biry(z).unwrap();
/// assert!((bi.re - 0.6149).abs() < 1e-3); // Bi(0) ≈ 0.6149
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z is not finite (NaN or infinite).
/// - [`Error::Overflow`] if |z| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| is too large for any significant digits (roughly > 10⁶ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn biry<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, Error> {
    let (result, _status) = airy::zbiry(z, AiryDerivative::Value, Scaling::Unscaled)?;
    Ok(result)
}

/// Derivative of the Airy function of the second kind, Bi'(z).
///
/// Computes the derivative of Bi(z) for complex z.
/// Bi'(z) grows super-exponentially for large positive real z, just as Bi(z) does.
/// Satisfies the differential equation Bi′′(z) = z · Bi(z).
///
/// # Example
///
/// ```
/// use complex_bessel::biryprime;
/// use num_complex::Complex;
///
/// let z = Complex::new(0.0_f64, 0.0);
/// let bi_prime = biryprime(z).unwrap();
/// assert!((bi_prime.re - 0.4483).abs() < 1e-3); // Bi'(0) ≈ 0.4483
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z is not finite (NaN or infinite).
/// - [`Error::Overflow`] if |z| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| is too large for any significant digits (roughly > 10⁶ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn biryprime<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, Error> {
    let (result, _status) = airy::zbiry(z, AiryDerivative::Derivative, Scaling::Unscaled)?;
    Ok(result)
}

// ── Scaled single-value functions ──

/// Scaled Bessel function of the first kind, exp(−|Im(z)|) · J<sub>ν</sub>(z).
///
/// The exponential factor cancels the asymptotic growth of J for large imaginary
/// arguments, keeping results in a representable floating-point range.
///
/// Supports negative ν via the same reflection formula as [`besselj`].
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Example
///
/// ```
/// use complex_bessel::besselj_scaled;
/// use num_complex::Complex;
///
/// let z = Complex::new(0.0_f64, 10.0);
/// let j_s = besselj_scaled(0.0, z).unwrap();
/// assert!((j_s.re - 0.1278).abs() < 1e-3); // exp(-|Im(z)|) * J_0(10i) ≈ 0.1278
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z or ν is not finite (NaN or infinite), or if z = 0 and ν is a negative non-integer.
/// - [`Error::Overflow`] if |z| or |ν| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn besselj_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    besselj_internal(nu, z, Scaling::Exponential)
}

/// Scaled Bessel function of the second kind, exp(−|Im(z)|) · Y<sub>ν</sub>(z).
///
/// The exponential factor cancels the asymptotic growth of Y for large imaginary
/// arguments, keeping results in a representable floating-point range.
///
/// Supports negative ν via the same reflection formula as [`bessely`].
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Example
///
/// ```
/// use complex_bessel::bessely_scaled;
/// use num_complex::Complex;
///
/// let z = Complex::new(0.0_f64, 10.0);
/// let y_s = bessely_scaled(0.0, z).unwrap();
/// assert!((y_s.im - 0.1278).abs() < 1e-3); // exp(-|Im(z)|) * Y_0(10i), Im ≈ 0.1278
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z or ν is not finite (NaN or infinite), or if z = 0.
/// - [`Error::Overflow`] if |z| is too small or too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn bessely_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    bessely_internal(nu, z, Scaling::Exponential)
}

/// Scaled modified Bessel function of the first kind, exp(−|Re(z)|) · I<sub>ν</sub>(z).
///
/// I_ν(z) grows exponentially for large |Re(z)|, so the unscaled value can
/// easily overflow. The scaling factor exp(−|Re(z)|) keeps the result finite.
///
/// Supports negative ν via the same reflection formula as [`besseli`].
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Example
///
/// ```
/// use complex_bessel::besseli_scaled;
/// use num_complex::Complex;
///
/// let z = Complex::new(10.0_f64, 0.0);
/// let i_s = besseli_scaled(0.0, z).unwrap();
/// assert!((i_s.re - 0.1278).abs() < 1e-3); // exp(-10) * I_0(10) ≈ 0.1278
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z or ν is not finite (NaN or infinite), or if z = 0 and ν is a negative non-integer.
/// - [`Error::Overflow`] if |z| or |ν| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn besseli_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    besseli_internal(nu, z, Scaling::Exponential)
}

/// Scaled modified Bessel function of the second kind, exp(z) · K<sub>ν</sub>(z).
///
/// K_ν(z) decays exponentially for large Re(z), so unscaled values can underflow
/// to zero. The scaling factor exp(z) keeps the result in a normal range.
///
/// Supports negative ν via the same reflection formula as [`besselk`].
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Example
///
/// ```
/// use complex_bessel::besselk_scaled;
/// use num_complex::Complex;
///
/// let z = Complex::new(500.0_f64, 0.0);
///
/// // K_0(500) would underflow to 0, but the scaled version stays finite:
/// let k_s = besselk_scaled(0.0, z).unwrap();
/// assert!((k_s.re - 0.0560).abs() < 1e-3); // exp(500) * K_0(500) ≈ 0.0560
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z or ν is not finite (NaN or infinite), or if z = 0.
/// - [`Error::Overflow`] if |z| is too small or too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn besselk_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    besselk_internal(nu, z, Scaling::Exponential)
}

/// Scaled Hankel function of the first kind, exp(−iz) · H<sub>ν</sub><sup>(1)</sup>(z).
///
/// H^(1) grows exponentially in the lower half-plane;
/// the scaling factor removes this growth, preventing overflow.
///
/// Supports negative ν via the same reflection formula as [`hankel1`].
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Example
///
/// ```
/// use complex_bessel::hankel1_scaled;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
/// let h_s = hankel1_scaled(0.0, z).unwrap();
/// assert!((h_s.re - 0.4877).abs() < 1e-3); // exp(-i) * H^(1)_0(1), Re ≈ 0.4877
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z or ν is not finite (NaN or infinite), or if z = 0.
/// - [`Error::Overflow`] if |z| is too small or too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn hankel1_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    hankel_internal(HankelKind::First, nu, z, Scaling::Exponential)
}

/// Scaled Hankel function of the second kind, exp(iz) · H<sub>ν</sub><sup>(2)</sup>(z).
///
/// H^(2) grows exponentially in the upper half-plane;
/// the scaling factor removes this growth, preventing overflow.
///
/// Supports negative ν via the same reflection formula as [`hankel2`].
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Example
///
/// ```
/// use complex_bessel::hankel2_scaled;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
/// let h_s = hankel2_scaled(0.0, z).unwrap();
/// assert!((h_s.re - 0.4877).abs() < 1e-3); // exp(i) * H^(2)_0(1), Re ≈ 0.4877
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z or ν is not finite (NaN or infinite), or if z = 0.
/// - [`Error::Overflow`] if |z| is too small or too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn hankel2_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    hankel_internal(HankelKind::Second, nu, z, Scaling::Exponential)
}

/// Scaled Airy function, exp(ζ) · Ai(z), where ζ = (2/3) z√z.
///
/// Ai(z) decays super-exponentially for large positive real z.
/// The scaling factor exp(ζ) keeps the result representable.
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Example
///
/// ```
/// use complex_bessel::airy_scaled;
/// use num_complex::Complex;
///
/// let z = Complex::new(5.0_f64, 0.0);
/// let ai_s = airy_scaled(z).unwrap();
/// assert!((ai_s.re - 0.1870).abs() < 1e-3); // exp(ζ) * Ai(5) ≈ 0.1870
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z is not finite (NaN or infinite).
/// - [`Error::Overflow`] if |z| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| is too large for any significant digits (roughly > 10⁶ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn airy_scaled<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, Error> {
    let (result, _nz, _status) = airy::zairy(z, AiryDerivative::Value, Scaling::Exponential)?;
    Ok(result)
}

/// Scaled derivative of the Airy function, exp(ζ) · Ai'(z), where ζ = (2/3) z√z.
///
/// Ai'(z) decays super-exponentially for large positive real z, just as Ai(z) does.
/// The scaling factor exp(ζ) keeps the result representable.
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Example
///
/// ```
/// use complex_bessel::airyprime_scaled;
/// use num_complex::Complex;
///
/// let z = Complex::new(5.0_f64, 0.0);
/// let ai_s = airyprime_scaled(z).unwrap();
/// assert!((ai_s.re + 0.4270).abs() < 1e-3); // exp(ζ) * Ai'(5) ≈ -0.4270
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z is not finite (NaN or infinite).
/// - [`Error::Overflow`] if |z| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| is too large for any significant digits (roughly > 10⁶ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn airyprime_scaled<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, Error> {
    let (result, _nz, _status) = airy::zairy(z, AiryDerivative::Derivative, Scaling::Exponential)?;
    Ok(result)
}

/// Scaled Airy function of the second kind, exp(−|Re(ζ)|) · Bi(z), where ζ = (2/3) z√z.
///
/// Bi(z) grows super-exponentially for large positive real z.
/// The scaling factor exp(−|Re(ζ)|) keeps the result representable.
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Example
///
/// ```
/// use complex_bessel::biry_scaled;
/// use num_complex::Complex;
///
/// let z = Complex::new(5.0_f64, 0.0);
/// let bi_s = biry_scaled(z).unwrap();
/// assert!((bi_s.re - 0.3811).abs() < 1e-3); // exp(-|Re(ζ)|) * Bi(5) ≈ 0.3811
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z is not finite (NaN or infinite).
/// - [`Error::Overflow`] if |z| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| is too large for any significant digits (roughly > 10⁶ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn biry_scaled<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, Error> {
    let (result, _status) = airy::zbiry(z, AiryDerivative::Value, Scaling::Exponential)?;
    Ok(result)
}

/// Scaled derivative of the Airy function of the second kind, exp(−|Re(ζ)|) · Bi'(z), where ζ = (2/3) z√z.
///
/// Bi'(z) grows super-exponentially for large positive real z, just as Bi(z) does.
/// The scaling factor exp(−|Re(ζ)|) keeps the result representable.
///
/// See [crate-level docs](crate#exponential-scaling) for the full scaling table.
///
/// # Example
///
/// ```
/// use complex_bessel::biryprime_scaled;
/// use num_complex::Complex;
///
/// let z = Complex::new(5.0_f64, 0.0);
/// let bi_s = biryprime_scaled(z).unwrap();
/// assert!((bi_s.re - 0.8319).abs() < 1e-3); // exp(-|Re(ζ)|) * Bi'(5) ≈ 0.8319
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z is not finite (NaN or infinite).
/// - [`Error::Overflow`] if |z| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| is too large for any significant digits (roughly > 10⁶ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn biryprime_scaled<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, Error> {
    let (result, _status) = airy::zbiry(z, AiryDerivative::Derivative, Scaling::Exponential)?;
    Ok(result)
}

// ── Airy _raw functions (expose Accuracy) ──

/// Airy function Ai(z) with accuracy status.
///
/// Like [`airy`], but returns an [`AiryResult`] that includes an [`Accuracy`] status:
/// - [`Accuracy::Normal`] — no significant precision loss
/// - [`Accuracy::Reduced`] — more than half of significant digits may be lost (|z| very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// # Example
///
/// ```
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(0.0_f64, 0.0);
/// let result = airy_raw(z, Scaling::Unscaled).unwrap();
/// assert!((result.value.re - 0.3550).abs() < 1e-3); // Ai(0) ≈ 0.3550
/// assert!(matches!(result.status, Accuracy::Normal));
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z is not finite (NaN or infinite).
/// - [`Error::Overflow`] if |z| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| is too large for any significant digits (roughly > 10⁶ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn airy_raw<T: BesselFloat>(z: Complex<T>, scaling: Scaling) -> Result<AiryResult<T>, Error> {
    let (value, _nz, status) = airy::zairy(z, AiryDerivative::Value, scaling)?;
    Ok(AiryResult { value, status })
}

/// Derivative of the Airy function Ai'(z) with accuracy status.
///
/// Like [`airyprime`], but returns an [`AiryResult`] that includes an [`Accuracy`] status:
/// - [`Accuracy::Normal`] — no significant precision loss
/// - [`Accuracy::Reduced`] — more than half of significant digits may be lost (|z| very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// # Example
///
/// ```
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(0.0_f64, 0.0);
/// let result = airyprime_raw(z, Scaling::Unscaled).unwrap();
/// assert!((result.value.re + 0.2588).abs() < 1e-3); // Ai'(0) ≈ -0.2588
/// assert!(matches!(result.status, Accuracy::Normal));
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z is not finite (NaN or infinite).
/// - [`Error::Overflow`] if |z| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| is too large for any significant digits (roughly > 10⁶ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn airyprime_raw<T: BesselFloat>(
    z: Complex<T>,
    scaling: Scaling,
) -> Result<AiryResult<T>, Error> {
    let (value, _nz, status) = airy::zairy(z, AiryDerivative::Derivative, scaling)?;
    Ok(AiryResult { value, status })
}

/// Airy function of the second kind Bi(z) with accuracy status.
///
/// Like [`biry`], but returns an [`AiryResult`] that includes an [`Accuracy`] status:
/// - [`Accuracy::Normal`] — no significant precision loss
/// - [`Accuracy::Reduced`] — more than half of significant digits may be lost (|z| very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// # Example
///
/// ```
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(0.0_f64, 0.0);
/// let result = biry_raw(z, Scaling::Unscaled).unwrap();
/// assert!((result.value.re - 0.6149).abs() < 1e-3); // Bi(0) ≈ 0.6149
/// assert!(matches!(result.status, Accuracy::Normal));
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z is not finite (NaN or infinite).
/// - [`Error::Overflow`] if |z| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| is too large for any significant digits (roughly > 10⁶ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn biry_raw<T: BesselFloat>(z: Complex<T>, scaling: Scaling) -> Result<AiryResult<T>, Error> {
    let (value, status) = airy::zbiry(z, AiryDerivative::Value, scaling)?;
    Ok(AiryResult { value, status })
}

/// Derivative of the Airy function of the second kind Bi'(z) with accuracy status.
///
/// Like [`biryprime`], but returns an [`AiryResult`] that includes an [`Accuracy`] status:
/// - [`Accuracy::Normal`] — no significant precision loss
/// - [`Accuracy::Reduced`] — more than half of significant digits may be lost (|z| very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// # Example
///
/// ```
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(0.0_f64, 0.0);
/// let result = biryprime_raw(z, Scaling::Unscaled).unwrap();
/// assert!((result.value.re - 0.4483).abs() < 1e-3); // Bi'(0) ≈ 0.4483
/// assert!(matches!(result.status, Accuracy::Normal));
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if z is not finite (NaN or infinite).
/// - [`Error::Overflow`] if |z| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| is too large for any significant digits (roughly > 10⁶ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
#[inline]
pub fn biryprime_raw<T: BesselFloat>(
    z: Complex<T>,
    scaling: Scaling,
) -> Result<AiryResult<T>, Error> {
    let (value, status) = airy::zbiry(z, AiryDerivative::Derivative, scaling)?;
    Ok(AiryResult { value, status })
}

// ── Sequence functions ──

/// Take the worse of two statuses (Reduced > Normal).
#[inline]
fn worse_status(a: Accuracy, b: Accuracy) -> Accuracy {
    match (a, b) {
        (Accuracy::Reduced, _) | (_, Accuracy::Reduced) => Accuracy::Reduced,
        _ => Accuracy::Normal,
    }
}

/// Underflow count of a reflected sequence: its zero elements. Underflow
/// can be at either end (large |orders| below zero, high orders above), so
/// all of them are counted (Amos also counts the exact zeros at z = 0).
#[inline]
fn recount_underflow<T: BesselFloat>(values: &[Complex<T>]) -> usize {
    let zero = Complex::new(T::zero(), T::zero());
    values.iter().filter(|v| **v == zero).count()
}

/// Compute the number of negative-order elements (before crossing zero).
/// For orders ν, ν+1, …, ν+n−1, this is the count where ν+j < 0.
#[inline]
fn neg_count<T: BesselFloat>(nu: T, n: usize) -> usize {
    if nu >= T::zero() {
        return 0;
    }
    // Number of j in [0,n) with nu + j < 0
    let abs_nu = nu.abs();
    // ceil(|nu|) gives the first j where nu+j >= 0
    let c = abs_nu.ceil();
    let cn = if let Some(v) = c.to_usize() { v } else { n };
    cn.min(n)
}

/// Upper-interface signature shared by the sequence helpers below:
/// `f(start, buf)` fills `buf` with the function at orders start, start+1, ….
type SeqFn<'a, T> = dyn FnMut(T, &mut [Complex<T>]) -> Result<(usize, Accuracy), Error> + 'a;

// ── Negative-order sequences ──
//
// Orders ν, ν+1, …, ν+n−1 with ν < 0 are computed at |order| and reflected.
// Everything happens in `out`; only the secondary function of the non-integer
// reflection formulas (and Y's K values) needs a workspace: the full length on
// the heap for `_seq`, or a fixed stack buffer filled in chunks for `_seq_into`.

/// Set `out[j]` = f(|ν+j|) when the |orders| form a single lattice: all
/// orders negative (nc == n), or ν an integer whose orders cross zero.
///
/// For a crossing, the orders 0..=max(|ν|, ν+n−1) always fit in `out`
/// (n = nc + number of orders ≥ 0), so one call computes them in place and
/// they are mirrored into position.
#[allow(clippy::needless_range_loop)]
fn fill_abs_orders<T: BesselFloat>(
    out: &mut [Complex<T>],
    nu: T,
    nc: usize,
    f: &mut SeqFn<'_, T>,
) -> Result<Accuracy, Error> {
    let n = out.len();
    if nc == n {
        // |orders| = |ν|-(n-1), …, |ν| ascending, then reversed into place
        let start = nu.abs() - T::from_f64((n - 1) as f64);
        let (_, status) = f(start, out)?;
        out.reverse();
        return Ok(status);
    }
    // Integer crossing, e.g. ν=-2, n=5 → orders -2,-1,0,1,2; nc = |ν|
    let pos_n = n - nc; // orders 0, 1, …, pos_n-1
    if pos_n > nc {
        // Largest |order| is positive: compute 0..pos_n into the positive part,
        // then mirror f(nc-j) into out[j]
        let (_, status) = f(T::zero(), &mut out[nc..])?;
        for j in 0..nc {
            out[j] = out[2 * nc - j];
        }
        Ok(status)
    } else {
        // Largest |order| is |ν|: compute 0..=nc, reverse so out[j] = f(nc-j),
        // then mirror f(m) into the positive tail
        let (_, status) = f(T::zero(), &mut out[..=nc])?;
        out[..=nc].reverse();
        for m in 1..pos_n {
            out[nc + m] = out[nc - m];
        }
        Ok(status)
    }
}

/// Set `out[..nc]` = f(|ν+j|) and `out[nc..]` = f(ν+j) for a non-integer ν
/// whose orders cross zero: two lattices, two calls.
fn fill_non_int_crossing<T: BesselFloat>(
    out: &mut [Complex<T>],
    nu: T,
    nc: usize,
    f: &mut SeqFn<'_, T>,
) -> Result<Accuracy, Error> {
    let abs_nu = nu.abs();
    // Negative part: |orders| frac, frac+1, …, |ν| ascending, then reversed
    let frac_neg = abs_nu - abs_nu.floor();
    let (_, status1) = f(frac_neg, &mut out[..nc])?;
    out[..nc].reverse();
    // Positive part: orders ν+nc, …, ν+n-1
    let pos_start = nu + T::from_f64(nc as f64);
    let (_, status2) = f(pos_start, &mut out[nc..])?;
    Ok(worse_status(status1, status2))
}

/// Apply a two-function reflection formula to the negative non-integer orders
/// −(|ν| − j), j = 0, …, len−1, stored in `seg`.
///
/// Their |orders| form the lattice lo, lo+1, …, |ν| with lo = |ν| − (len−1),
/// computed in ascending order and then reversed. `whole` computes one
/// function for the whole lattice into `seg`; `chunked` computes the other
/// into `buf`, `buf.len()` orders at a time. With `chunked_first`, the first
/// chunk is computed before `whole` (to keep the call order of the formula).
/// `combine(|order|, whole, chunked)` returns the reflected value.
fn reflect_lattice<T: BesselFloat>(
    seg: &mut [Complex<T>],
    abs_nu: T,
    chunked_first: bool,
    whole: &mut SeqFn<'_, T>,
    chunked: &mut SeqFn<'_, T>,
    combine: &dyn Fn(T, Complex<T>, Complex<T>) -> Complex<T>,
    buf: &mut [Complex<T>],
) -> Result<Accuracy, Error> {
    let len = seg.len();
    let chunk = buf.len();
    assert!(chunk > 0, "empty workspace");
    // |ν| minus an integer below it stays on the ulp grid of |ν|: exact
    let lo = abs_nu - T::from_f64((len - 1) as f64);
    let mut status = Accuracy::Normal;
    let mut whole_done = false;
    if !chunked_first {
        let (_, s) = whole(lo, seg)?;
        status = s;
        whole_done = true;
    }
    let mut c0 = 0;
    while c0 < len {
        let c1 = (c0 + chunk).min(len);
        // Exact for the same reason as `lo`
        let start = lo + T::from_f64(c0 as f64);
        let (_, s) = chunked(start, &mut buf[..c1 - c0])?;
        status = worse_status(status, s);
        if !whole_done {
            let (_, s) = whole(lo, seg)?;
            status = worse_status(status, s);
            whole_done = true;
        }
        for k in c0..c1 {
            let order = abs_nu - T::from_f64((len - 1 - k) as f64);
            seg[k] = combine(order, seg[k], buf[k - c0]);
        }
        c0 = c1;
    }
    seg.reverse();
    Ok(status)
}

/// Multiply the negative integer orders by (−1)^|order| (J and Y).
fn apply_integer_sign<T: BesselFloat>(out: &mut [Complex<T>], nu: T) {
    for (j, v) in out.iter_mut().enumerate() {
        let order = nu + T::from_f64(j as f64);
        if order < T::zero() {
            let int_order = order.abs().to_i64().unwrap_or(0);
            *v = *v * integer_sign::<T>(int_order);
        }
    }
}

/// `n = 0` would underflow the `n - 1` order arithmetic of the helpers above.
#[inline]
fn check_nonempty<T>(out: &[T]) -> Result<(), Error> {
    if out.is_empty() {
        Err(Error::InvalidInput)
    } else {
        Ok(())
    }
}

/// K with ν < 0: K_{−ν} = K_ν (even function).
fn besselk_seq_neg<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    out: &mut [Complex<T>],
) -> Result<SeqInfo, Error> {
    check_nonempty(out)?;
    let nc = neg_count(nu, out.len());
    let f = &mut |start, buf: &mut [Complex<T>]| besk::zbesk(z, start, scaling, buf);
    let status = if nc == out.len() || as_integer(nu.abs()).is_some() {
        fill_abs_orders(out, nu, nc, f)?
    } else {
        fill_non_int_crossing(out, nu, nc, f)?
    };
    Ok(SeqInfo {
        underflow_count: recount_underflow(out),
        status,
    })
}

/// H with ν < 0: H^(1)_{−ν} = exp(νπi)·H^(1)_ν, H^(2)_{−ν} = exp(−νπi)·H^(2)_ν.
fn hankel_seq_neg<T: BesselFloat>(
    kind: HankelKind,
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    out: &mut [Complex<T>],
) -> Result<SeqInfo, Error> {
    check_nonempty(out)?;
    let abs_nu = nu.abs();
    let nc = neg_count(nu, out.len());
    let f = &mut |start, buf: &mut [Complex<T>]| besh::zbesh(z, start, kind, scaling, buf);
    let status = if nc == out.len() || as_integer(abs_nu).is_some() {
        fill_abs_orders(out, nu, nc, f)?
    } else {
        fill_non_int_crossing(out, nu, nc, f)?
    };
    // Rotate the negative orders; |ν+j| = |ν| - j exactly
    for (j, v) in out.iter_mut().take(nc).enumerate() {
        let order = abs_nu - T::from_f64(j as f64);
        *v = reflect_h_element(order, kind, *v);
    }
    Ok(SeqInfo {
        underflow_count: recount_underflow(out),
        status,
    })
}

/// J with ν < 0: J_{−ν} = cos(νπ)·J_ν − sin(νπ)·Y_ν, J_{−n} = (−1)^n·J_n.
fn besselj_seq_neg<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    out: &mut [Complex<T>],
    ws: Workspace<'_, T>,
) -> Result<SeqInfo, Error> {
    check_nonempty(out)?;
    let n = out.len();
    let nc = neg_count(nu, n);
    let status = if as_integer(nu.abs()).is_some() {
        let f = &mut |start, buf: &mut [Complex<T>]| besj::zbesj(z, start, scaling, buf);
        let status = fill_abs_orders(out, nu, nc, f)?;
        apply_integer_sign(out, nu);
        status
    } else {
        // J into `out`, Y into one half of the workspace and Y's K values
        // into the other (J first, as before)
        with_workspace(ws, 2 * n, |ws| {
            let (neg, pos) = out.split_at_mut(nc);
            let (y_buf, k_buf) = ws.split_at_mut(ws.len() / 2);
            let status1 = reflect_lattice(
                neg,
                nu.abs(),
                false,
                &mut |start, buf| besj::zbesj(z, start, scaling, buf),
                &mut |start, buf| {
                    besy::zbesy_with(z, start, scaling, buf, Workspace::Slice(&mut *k_buf))
                },
                &|order, j, y| reflect_j_element(order, j, y),
                y_buf,
            )?;
            if pos.is_empty() {
                return Ok(status1);
            }
            let (_, status2) = besj::zbesj(z, nu + T::from_f64(nc as f64), scaling, pos)?;
            Ok(worse_status(status1, status2))
        })?
    };
    Ok(SeqInfo {
        underflow_count: recount_underflow(out),
        status,
    })
}

/// Y with ν < 0: Y_{−ν} = sin(νπ)·J_ν + cos(νπ)·Y_ν, Y_{−n} = (−1)^n·Y_n.
fn bessely_seq_neg<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    out: &mut [Complex<T>],
    ws: Workspace<'_, T>,
) -> Result<SeqInfo, Error> {
    check_nonempty(out)?;
    let n = out.len();
    let nc = neg_count(nu, n);
    let is_int = as_integer(nu.abs()).is_some();
    // Heap: Y's K values only, or also the Y values of the reflection formula
    let ws_len = if is_int { n } else { 2 * n };
    let status = with_workspace(ws, ws_len, |ws| {
        if is_int {
            let f = &mut |start, buf: &mut [Complex<T>]| {
                besy::zbesy_with(z, start, scaling, buf, Workspace::Slice(&mut *ws))
            };
            let status = fill_abs_orders(out, nu, nc, f)?;
            apply_integer_sign(out, nu);
            return Ok(status);
        }
        // J into `out`, Y into one half of the workspace and Y's K values
        // into the other (J first, as before)
        let (neg, pos) = out.split_at_mut(nc);
        let status1 = {
            let (y_buf, k_buf) = ws.split_at_mut(ws.len() / 2);
            reflect_lattice(
                neg,
                nu.abs(),
                false,
                &mut |start, buf| besj::zbesj(z, start, scaling, buf),
                &mut |start, buf| {
                    besy::zbesy_with(z, start, scaling, buf, Workspace::Slice(&mut *k_buf))
                },
                &|order, j, y| reflect_y_element(order, j, y),
                y_buf,
            )?
        };
        if pos.is_empty() {
            return Ok(status1);
        }
        let pos_start = nu + T::from_f64(nc as f64);
        let (_, status2) = besy::zbesy_with(z, pos_start, scaling, pos, Workspace::Slice(ws))?;
        Ok(worse_status(status1, status2))
    })?;
    Ok(SeqInfo {
        underflow_count: recount_underflow(out),
        status,
    })
}

/// I with ν < 0: I_{−ν} = I_ν + (2/π)·sin(νπ)·K_ν, I_{−n} = I_n.
fn besseli_seq_neg<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    out: &mut [Complex<T>],
    ws: Workspace<'_, T>,
) -> Result<SeqInfo, Error> {
    check_nonempty(out)?;
    let n = out.len();
    let nc = neg_count(nu, n);
    let i_fn = &mut |start, buf: &mut [Complex<T>]| besi::zbesi(z, start, scaling, buf);
    let status = if as_integer(nu.abs()).is_some() {
        fill_abs_orders(out, nu, nc, i_fn)?
    } else {
        // K into `out` in one call, I into the workspace: Amos's K overflow
        // test uses the last order and length of the call, so splitting K
        // could fail where the whole call succeeds, while I's tests only get
        // more lenient in later chunks. The first I chunk goes first to keep
        // the I-then-K order.
        with_workspace(ws, n, |ws| {
            let (neg, pos) = out.split_at_mut(nc);
            let status1 = reflect_lattice(
                neg,
                nu.abs(),
                true,
                &mut |start, buf| besk::zbesk(z, start, scaling, buf),
                i_fn,
                &|order, k, i| {
                    let k = if scaling == Scaling::Exponential {
                        k_to_i_scaling_correction(z, k)
                    } else {
                        k
                    };
                    reflect_i_element(order, i, k)
                },
                ws,
            )?;
            if pos.is_empty() {
                return Ok(status1);
            }
            let (_, status2) = besi::zbesi(z, nu + T::from_f64(nc as f64), scaling, pos)?;
            Ok(worse_status(status1, status2))
        })?
    };
    Ok(SeqInfo {
        underflow_count: recount_underflow(out),
        status,
    })
}

// ── Shared by `_seq_into` (stack workspace) and `_seq` (heap workspace) ──

fn besselj_seq_impl<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    out: &mut [Complex<T>],
    ws: Workspace<'_, T>,
) -> Result<SeqInfo, Error> {
    if nu < T::zero() {
        return besselj_seq_neg(nu, z, scaling, out, ws);
    }
    let (underflow_count, status) = besj::zbesj(z, nu, scaling, out)?;
    Ok(SeqInfo {
        underflow_count,
        status,
    })
}

fn bessely_seq_impl<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    out: &mut [Complex<T>],
    ws: Workspace<'_, T>,
) -> Result<SeqInfo, Error> {
    if nu < T::zero() {
        return bessely_seq_neg(nu, z, scaling, out, ws);
    }
    let (underflow_count, status) = besy::zbesy_with(z, nu, scaling, out, ws)?;
    Ok(SeqInfo {
        underflow_count,
        status,
    })
}

fn besseli_seq_impl<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    out: &mut [Complex<T>],
    ws: Workspace<'_, T>,
) -> Result<SeqInfo, Error> {
    if nu < T::zero() {
        return besseli_seq_neg(nu, z, scaling, out, ws);
    }
    let (underflow_count, status) = besi::zbesi(z, nu, scaling, out)?;
    Ok(SeqInfo {
        underflow_count,
        status,
    })
}

fn hankel_seq_impl<T: BesselFloat>(
    kind: HankelKind,
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    out: &mut [Complex<T>],
) -> Result<SeqInfo, Error> {
    if nu < T::zero() {
        return hankel_seq_neg(kind, nu, z, scaling, out);
    }
    let (underflow_count, status) = besh::zbesh(z, nu, kind, scaling, out)?;
    Ok(SeqInfo {
        underflow_count,
        status,
    })
}

/// Compute J<sub>ν+j</sub>(z) for j = 0, 1, …, n−1 into `out`, where n = `out.len()`.
///
/// Like [`besselj_seq`], but writes into a caller-provided slice instead of
/// allocating, so it also works without an allocator (`no_std`). Returns a
/// [`SeqInfo`] with the underflow count and [`Accuracy`]. If an error is
/// returned, the contents of `out` are unspecified.
///
/// For negative non-integer ν, the Y values of the reflection formula need
/// workspace. It is a 1 KB stack buffer (f64), so long sequences are computed
/// in chunks. That is somewhat slower than [`besselj_seq`], and the last digits
/// can differ (within the accuracy of the single-value functions); close to
/// the overflow limit, one of the two can even return an error where the
/// other succeeds.
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Supports negative ν via the same reflection formula as [`besselj`].
///
/// # Example
///
/// ```
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
///
/// // J_0(z), J_1(z), J_2(z) in one call
/// let mut out = [Complex::new(0.0, 0.0); 3];
/// let info = besselj_seq_into(0.0, z, Scaling::Unscaled, &mut out).unwrap();
/// assert!((out[0].re - 0.7652).abs() < 1e-3); // J_0(1) ≈ 0.7652
/// assert_eq!(info.status, Accuracy::Normal);
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if `out` is empty, if z or ν is not finite (NaN or infinite), or if z = 0 and ν is a negative non-integer.
/// - [`Error::Overflow`] if |z| or |ν| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
pub fn besselj_seq_into<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    out: &mut [Complex<T>],
) -> Result<SeqInfo, Error> {
    besselj_seq_impl(nu, z, scaling, out, Workspace::Stack)
}

/// Compute Y<sub>ν+j</sub>(z) for j = 0, 1, …, n−1 into `out`, where n = `out.len()`.
///
/// Like [`bessely_seq`], but writes into a caller-provided slice instead of
/// allocating, so it also works without an allocator (`no_std`). Returns a
/// [`SeqInfo`] with the underflow count and [`Accuracy`]. If an error is
/// returned, the contents of `out` are unspecified.
///
/// The K values that Y is computed from (and, for negative non-integer ν, the
/// second function of the reflection formula) need workspace. It is a 1 KB
/// stack buffer (f64), so long sequences are computed in chunks. That is
/// somewhat slower than [`bessely_seq`], and the last digits can differ
/// (within the accuracy of the single-value functions); close to the overflow
/// limit, one of the two can even return an error where the other succeeds.
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Supports negative ν via the same reflection formula as [`bessely`].
///
/// # Example
///
/// ```
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
///
/// // Y_0(z), Y_1(z), Y_2(z) in one call
/// let mut out = [Complex::new(0.0, 0.0); 3];
/// let info = bessely_seq_into(0.0, z, Scaling::Unscaled, &mut out).unwrap();
/// assert!((out[0].re - 0.0883).abs() < 1e-3); // Y_0(1) ≈ 0.0883
/// assert_eq!(info.status, Accuracy::Normal);
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if `out` is empty, if z or ν is not finite (NaN or infinite), or if z = 0.
/// - [`Error::Overflow`] if |z| is too small or too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
pub fn bessely_seq_into<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    out: &mut [Complex<T>],
) -> Result<SeqInfo, Error> {
    bessely_seq_impl(nu, z, scaling, out, Workspace::Stack)
}

/// Compute I<sub>ν+j</sub>(z) for j = 0, 1, …, n−1 into `out`, where n = `out.len()`.
///
/// Like [`besseli_seq`], but writes into a caller-provided slice instead of
/// allocating, so it also works without an allocator (`no_std`). Returns a
/// [`SeqInfo`] with the underflow count and [`Accuracy`]. If an error is
/// returned, the contents of `out` are unspecified.
///
/// For negative non-integer ν, the I values of the reflection formula need
/// workspace. It is a 1 KB stack buffer (f64), so long sequences are computed
/// in chunks. That is somewhat slower than [`besseli_seq`], and the last digits
/// can differ (within the accuracy of the single-value functions); close to
/// the overflow limit, one of the two can even return an error where the
/// other succeeds.
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Supports negative ν via the same reflection formula as [`besseli`].
///
/// # Example
///
/// ```
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
///
/// // I_0(z), I_1(z), I_2(z) in one call
/// let mut out = [Complex::new(0.0, 0.0); 3];
/// let info = besseli_seq_into(0.0, z, Scaling::Unscaled, &mut out).unwrap();
/// assert!((out[0].re - 1.2661).abs() < 1e-3); // I_0(1) ≈ 1.2661
/// assert_eq!(info.status, Accuracy::Normal);
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if `out` is empty, if z or ν is not finite (NaN or infinite), or if z = 0 and ν is a negative non-integer.
/// - [`Error::Overflow`] if |z| or |ν| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
pub fn besseli_seq_into<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    out: &mut [Complex<T>],
) -> Result<SeqInfo, Error> {
    besseli_seq_impl(nu, z, scaling, out, Workspace::Stack)
}

/// Compute K<sub>ν+j</sub>(z) for j = 0, 1, …, n−1 into `out`, where n = `out.len()`.
///
/// Like [`besselk_seq`], but writes into a caller-provided slice instead of
/// allocating, so it also works without an allocator (`no_std`). Returns a
/// [`SeqInfo`] with the underflow count and [`Accuracy`]. If an error is
/// returned, the contents of `out` are unspecified.
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Supports negative ν via the same reflection formula as [`besselk`].
///
/// # Example
///
/// ```
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
///
/// // K_0(z), K_1(z), K_2(z) in one call
/// let mut out = [Complex::new(0.0, 0.0); 3];
/// let info = besselk_seq_into(0.0, z, Scaling::Unscaled, &mut out).unwrap();
/// assert!((out[0].re - 0.4211).abs() < 1e-3); // K_0(1) ≈ 0.4211
/// assert_eq!(info.status, Accuracy::Normal);
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if `out` is empty, if z or ν is not finite (NaN or infinite), or if z = 0.
/// - [`Error::Overflow`] if |z| is too small or too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
pub fn besselk_seq_into<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    out: &mut [Complex<T>],
) -> Result<SeqInfo, Error> {
    if nu < T::zero() {
        return besselk_seq_neg(nu, z, scaling, out);
    }
    let (underflow_count, status) = besk::zbesk(z, nu, scaling, out)?;
    Ok(SeqInfo {
        underflow_count,
        status,
    })
}

/// Compute H<sub>ν+j</sub><sup>(1)</sup>(z) for j = 0, 1, …, n−1 into `out`, where n = `out.len()`.
///
/// Like [`hankel1_seq`], but writes into a caller-provided slice instead of
/// allocating, so it also works without an allocator (`no_std`). Returns a
/// [`SeqInfo`] with the underflow count and [`Accuracy`]. If an error is
/// returned, the contents of `out` are unspecified.
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Supports negative ν via the same reflection formula as [`hankel1`].
///
/// # Example
///
/// ```
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
///
/// // H^(1)_0(z), H^(1)_1(z) in one call
/// let mut out = [Complex::new(0.0, 0.0); 2];
/// let info = hankel1_seq_into(0.0, z, Scaling::Unscaled, &mut out).unwrap();
/// assert!((out[0].im - 0.0883).abs() < 1e-3); // Im = Y_0(1) ≈ 0.0883
/// assert_eq!(info.status, Accuracy::Normal);
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if `out` is empty, if z or ν is not finite (NaN or infinite), or if z = 0.
/// - [`Error::Overflow`] if |z| is too small or too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
pub fn hankel1_seq_into<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    out: &mut [Complex<T>],
) -> Result<SeqInfo, Error> {
    hankel_seq_impl(HankelKind::First, nu, z, scaling, out)
}

/// Compute H<sub>ν+j</sub><sup>(2)</sup>(z) for j = 0, 1, …, n−1 into `out`, where n = `out.len()`.
///
/// Like [`hankel2_seq`], but writes into a caller-provided slice instead of
/// allocating, so it also works without an allocator (`no_std`). Returns a
/// [`SeqInfo`] with the underflow count and [`Accuracy`]. If an error is
/// returned, the contents of `out` are unspecified.
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Supports negative ν via the same reflection formula as [`hankel2`].
///
/// # Example
///
/// ```
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
///
/// // H^(2)_0(z), H^(2)_1(z) in one call
/// let mut out = [Complex::new(0.0, 0.0); 2];
/// let info = hankel2_seq_into(0.0, z, Scaling::Unscaled, &mut out).unwrap();
/// assert!((out[0].im + 0.0883).abs() < 1e-3); // Im = -Y_0(1) ≈ -0.0883
/// assert_eq!(info.status, Accuracy::Normal);
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if `out` is empty, if z or ν is not finite (NaN or infinite), or if z = 0.
/// - [`Error::Overflow`] if |z| is too small or too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
pub fn hankel2_seq_into<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    out: &mut [Complex<T>],
) -> Result<SeqInfo, Error> {
    hankel_seq_impl(HankelKind::Second, nu, z, scaling, out)
}

// ── Allocating sequence functions ──

#[cfg(feature = "std")]
fn seq_helper<T: BesselFloat>(
    n: usize,
    f: impl FnOnce(&mut [Complex<T>]) -> Result<SeqInfo, Error>,
) -> Result<BesselResult<T>, Error> {
    let zero = T::zero();
    let mut values = alloc::vec![Complex::new(zero, zero); n];
    let info = f(&mut values)?;
    Ok(BesselResult {
        values,
        underflow_count: info.underflow_count,
        status: info.status,
    })
}

#[cfg(feature = "std")]
/// Compute J<sub>ν+j</sub>(z) for j = 0, 1, …, n−1 in a single call.
///
/// Allocating version of [`besselj_seq_into`].
///
/// Returns a [`BesselResult`] containing `n` values and an [`Accuracy`]:
/// - [`Accuracy::Normal`] — no significant precision loss
/// - [`Accuracy::Reduced`] — more than half of significant digits may be lost (|z| or |ν| very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Supports negative ν via the same reflection formula as [`besselj`].
///
/// # Example
///
/// ```
/// # #[cfg(feature = "std")] {
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
///
/// // J_0(z), J_1(z), J_2(z) in one call
/// let result = besselj_seq(0.0, z, 3, Scaling::Unscaled).unwrap();
/// assert_eq!(result.values.len(), 3);
/// assert!((result.values[0].re - 0.7652).abs() < 1e-3); // J_0(1) ≈ 0.7652
/// # }
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if n < 1, if z or ν is not finite (NaN or infinite), or if z = 0 and ν is a negative non-integer.
/// - [`Error::Overflow`] if |z| or |ν| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
pub fn besselj_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    seq_helper(n, |out| {
        besselj_seq_impl(nu, z, scaling, out, Workspace::Heap)
    })
}

#[cfg(feature = "std")]
/// Compute Y<sub>ν+j</sub>(z) for j = 0, 1, …, n−1 in a single call.
///
/// Allocating version of [`bessely_seq_into`].
///
/// Returns a [`BesselResult`] containing `n` values and an [`Accuracy`]:
/// - [`Accuracy::Normal`] — no significant precision loss
/// - [`Accuracy::Reduced`] — more than half of significant digits may be lost (|z| or |ν| very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Supports negative ν via the same reflection formula as [`bessely`].
///
/// # Example
///
/// ```
/// # #[cfg(feature = "std")] {
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
///
/// // Y_0(z), Y_1(z), Y_2(z) in one call
/// let result = bessely_seq(0.0, z, 3, Scaling::Unscaled).unwrap();
/// assert_eq!(result.values.len(), 3);
/// assert!((result.values[0].re - 0.0883).abs() < 1e-3); // Y_0(1) ≈ 0.0883
/// # }
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if n < 1, if z or ν is not finite (NaN or infinite), or if z = 0.
/// - [`Error::Overflow`] if |z| is too small or too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
pub fn bessely_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    seq_helper(n, |out| {
        bessely_seq_impl(nu, z, scaling, out, Workspace::Heap)
    })
}

#[cfg(feature = "std")]
/// Compute I<sub>ν+j</sub>(z) for j = 0, 1, …, n−1 in a single call.
///
/// Allocating version of [`besseli_seq_into`].
///
/// Returns a [`BesselResult`] containing `n` values and an [`Accuracy`]:
/// - [`Accuracy::Normal`] — no significant precision loss
/// - [`Accuracy::Reduced`] — more than half of significant digits may be lost (|z| or |ν| very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Supports negative ν via the same reflection formula as [`besseli`].
///
/// # Example
///
/// ```
/// # #[cfg(feature = "std")] {
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
///
/// // I_0(z), I_1(z), I_2(z) in one call
/// let result = besseli_seq(0.0, z, 3, Scaling::Unscaled).unwrap();
/// assert_eq!(result.values.len(), 3);
/// assert!((result.values[0].re - 1.2661).abs() < 1e-3); // I_0(1) ≈ 1.2661
/// # }
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if n < 1, if z or ν is not finite (NaN or infinite), or if z = 0 and ν is a negative non-integer.
/// - [`Error::Overflow`] if |z| or |ν| is too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
pub fn besseli_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    seq_helper(n, |out| {
        besseli_seq_impl(nu, z, scaling, out, Workspace::Heap)
    })
}

#[cfg(feature = "std")]
/// Compute K<sub>ν+j</sub>(z) for j = 0, 1, …, n−1 in a single call.
///
/// Allocating version of [`besselk_seq_into`].
///
/// Returns a [`BesselResult`] containing `n` values and an [`Accuracy`]:
/// - [`Accuracy::Normal`] — no significant precision loss
/// - [`Accuracy::Reduced`] — more than half of significant digits may be lost (|z| or |ν| very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Supports negative ν via the same reflection formula as [`besselk`].
///
/// # Example
///
/// ```
/// # #[cfg(feature = "std")] {
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
///
/// // K_0(z), K_1(z), K_2(z) in one call
/// let result = besselk_seq(0.0, z, 3, Scaling::Unscaled).unwrap();
/// assert_eq!(result.values.len(), 3);
/// assert!((result.values[0].re - 0.4211).abs() < 1e-3); // K_0(1) ≈ 0.4211
/// # }
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if n < 1, if z or ν is not finite (NaN or infinite), or if z = 0.
/// - [`Error::Overflow`] if |z| is too small or too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
pub fn besselk_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    seq_helper(n, |out| besselk_seq_into(nu, z, scaling, out))
}

#[cfg(feature = "std")]
/// Compute H<sub>ν+j</sub><sup>(1)</sup>(z) for j = 0, 1, …, n−1 in a single call.
///
/// Allocating version of [`hankel1_seq_into`].
///
/// Returns a [`BesselResult`] containing `n` values and an [`Accuracy`]:
/// - [`Accuracy::Normal`] — no significant precision loss
/// - [`Accuracy::Reduced`] — more than half of significant digits may be lost (|z| or |ν| very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Supports negative ν via the same reflection formula as [`hankel1`].
///
/// # Example
///
/// ```
/// # #[cfg(feature = "std")] {
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
///
/// // H^(1)_0(z), H^(1)_1(z) in one call
/// let result = hankel1_seq(0.0, z, 2, Scaling::Unscaled).unwrap();
/// assert_eq!(result.values.len(), 2);
/// assert!((result.values[0].im - 0.0883).abs() < 1e-3); // Im = Y_0(1) ≈ 0.0883
/// # }
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if n < 1, if z or ν is not finite (NaN or infinite), or if z = 0.
/// - [`Error::Overflow`] if |z| is too small or too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
pub fn hankel1_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    seq_helper(n, |out| hankel1_seq_into(nu, z, scaling, out))
}

#[cfg(feature = "std")]
/// Compute H<sub>ν+j</sub><sup>(2)</sup>(z) for j = 0, 1, …, n−1 in a single call.
///
/// Allocating version of [`hankel2_seq_into`].
///
/// Returns a [`BesselResult`] containing `n` values and an [`Accuracy`]:
/// - [`Accuracy::Normal`] — no significant precision loss
/// - [`Accuracy::Reduced`] — more than half of significant digits may be lost (|z| or |ν| very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Supports negative ν via the same reflection formula as [`hankel2`].
///
/// # Example
///
/// ```
/// # #[cfg(feature = "std")] {
/// use complex_bessel::*;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
///
/// // H^(2)_0(z), H^(2)_1(z) in one call
/// let result = hankel2_seq(0.0, z, 2, Scaling::Unscaled).unwrap();
/// assert_eq!(result.values.len(), 2);
/// assert!((result.values[0].im + 0.0883).abs() < 1e-3); // Im = -Y_0(1) ≈ -0.0883
/// # }
/// ```
///
/// # Errors
///
/// - [`Error::InvalidInput`] if n < 1, if z or ν is not finite (NaN or infinite), or if z = 0.
/// - [`Error::Overflow`] if |z| is too small or too large for a finite result.
/// - [`Error::TotalPrecisionLoss`] if |z| or |ν| is too large for any significant digits (roughly > 10⁹ for f64).
/// - [`Error::ConvergenceFailure`] if an internal series or recurrence does not converge (rare).
pub fn hankel2_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    seq_helper(n, |out| hankel2_seq_into(nu, z, scaling, out))
}

// ── Tests for non-finite input rejection at the public API level ──
// Deliberately gated on test only (not std), so they run in every
// feature configuration.

#[cfg(test)]
mod non_finite_tests {
    use super::*;
    use num_complex::{Complex32, Complex64};

    #[test]
    fn bessel_non_finite_z_returns_invalid_input() {
        let nan = f64::NAN;
        for z in [
            Complex64::new(nan, 0.0),
            Complex64::new(0.0, nan),
            Complex64::new(nan, nan),
            Complex64::new(f64::INFINITY, 0.0),
        ] {
            assert!(matches!(besselj(0.0, z), Err(Error::InvalidInput)));
            assert!(matches!(besselj_scaled(1.0, z), Err(Error::InvalidInput)));
            assert!(matches!(bessely(0.0, z), Err(Error::InvalidInput)));
            assert!(matches!(besseli(0.0, z), Err(Error::InvalidInput)));
            assert!(matches!(besselk(0.0, z), Err(Error::InvalidInput)));
            assert!(matches!(hankel1(0.0, z), Err(Error::InvalidInput)));
            assert!(matches!(hankel2(0.0, z), Err(Error::InvalidInput)));
            // Negative order exercises the reflection path in the wrappers.
            assert!(matches!(besselj(-2.5, z), Err(Error::InvalidInput)));
            assert!(matches!(besselk(-2.5, z), Err(Error::InvalidInput)));
        }
    }

    #[test]
    fn bessel_non_finite_nu_returns_invalid_input() {
        let z = Complex64::new(1.0, 1.0);
        for nu in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(besselj(nu, z), Err(Error::InvalidInput)));
            assert!(matches!(bessely(nu, z), Err(Error::InvalidInput)));
            assert!(matches!(besseli(nu, z), Err(Error::InvalidInput)));
            assert!(matches!(besselk(nu, z), Err(Error::InvalidInput)));
            assert!(matches!(hankel1(nu, z), Err(Error::InvalidInput)));
            assert!(matches!(hankel2(nu, z), Err(Error::InvalidInput)));
        }
    }

    #[test]
    fn bessel_non_finite_f32_returns_invalid_input() {
        let z = Complex32::new(1.0, 1.0);
        assert!(matches!(
            besselj(0.0f32, Complex32::new(f32::NAN, 0.0)),
            Err(Error::InvalidInput)
        ));
        assert!(matches!(besselj(f32::NAN, z), Err(Error::InvalidInput)));
    }

    #[test]
    fn airy_non_finite_z_returns_invalid_input() {
        let nan = f64::NAN;
        for z in [
            Complex64::new(nan, 0.0),
            Complex64::new(0.0, nan),
            Complex64::new(f64::INFINITY, 0.0),
        ] {
            assert!(matches!(airy(z), Err(Error::InvalidInput)));
            assert!(matches!(airyprime(z), Err(Error::InvalidInput)));
            assert!(matches!(biry(z), Err(Error::InvalidInput)));
            assert!(matches!(biryprime(z), Err(Error::InvalidInput)));
            assert!(matches!(
                airy_raw(z, Scaling::Unscaled),
                Err(Error::InvalidInput)
            ));
            assert!(matches!(
                biry_raw(z, Scaling::Unscaled),
                Err(Error::InvalidInput)
            ));
        }
    }

    #[test]
    fn finite_overflowing_modulus_keeps_total_precision_loss() {
        // |z| overflows to infinity while both components are finite: this
        // must stay TotalPrecisionLoss (range check), not InvalidInput.
        let z = Complex64::new(1.7e308, 1.7e308);
        assert!(matches!(besselj(0.0, z), Err(Error::TotalPrecisionLoss)));
    }

    #[test]
    fn seq_non_finite_inputs_return_invalid_input() {
        let z_nan = Complex64::new(f64::NAN, 0.0);
        let z = Complex64::new(1.0, 1.0);
        let mut out = [Complex64::new(0.0, 0.0); 3];
        assert!(matches!(
            besselj_seq_into(f64::NAN, z, Scaling::Unscaled, &mut out[..2]),
            Err(Error::InvalidInput)
        ));
        assert!(matches!(
            besselj_seq_into(0.0, z_nan, Scaling::Unscaled, &mut out[..2]),
            Err(Error::InvalidInput)
        ));
        assert!(matches!(
            besselk_seq_into(-2.5, z_nan, Scaling::Unscaled, &mut out),
            Err(Error::InvalidInput)
        ));
    }
}

// ── Scaled vs unscaled consistency ──
//
// Every `_scaled` function is the unscaled one times an elementary factor, so
// the two must agree wherever both are representable. The grid reaches the
// regions where the scaled paths switch to underflow/overflow handling
// (2|Im z| above ELIM for Y, large orders), which a value comparison at
// moderate arguments does not.

#[cfg(test)]
mod scaling_tests {
    use super::*;
    use num_complex::Complex64;

    type F<T> = fn(T, Complex<T>) -> Result<Complex<T>, Error>;

    fn fns<T: BesselFloat>() -> [(&'static str, F<T>, F<T>); 6] {
        [
            ("J", besselj, besselj_scaled),
            ("Y", bessely, bessely_scaled),
            ("I", besseli, besseli_scaled),
            ("K", besselk, besselk_scaled),
            ("H1", hankel1, hankel1_scaled),
            ("H2", hankel2, hankel2_scaled),
        ]
    }

    /// Log of the factor in scaled = exp(c)·unscaled.
    fn log_factor(name: &str, z: Complex64) -> Complex64 {
        let i = Complex64::new(0.0, 1.0);
        match name {
            "J" | "Y" => Complex64::new(-z.im.abs(), 0.0),
            "I" => Complex64::new(-z.re.abs(), 0.0),
            "K" => z,
            "H1" => -i * z,
            _ => i * z,
        }
    }

    /// Compare scaled with exp(c)·unscaled, formed in log space so neither
    /// the factor nor the product overflows. Pairs where either value fails,
    /// is non-finite, or lies below `floor` (underflow handling) are skipped.
    fn check<T: BesselFloat>(nus: &[f64], radii: &[f64], floor: f64, tol: f64) {
        let to64 = |v: Complex<T>| Complex64::new(v.re.to_f64().unwrap(), v.im.to_f64().unwrap());
        let mut compared = 0;
        for (name, unscaled, scaled) in fns::<T>() {
            for &nu in nus {
                for &r in radii {
                    for k in 0..12 {
                        // 12 angles over the full circle, including both axes
                        let th =
                            -core::f64::consts::PI + core::f64::consts::FRAC_PI_6 * f64::from(k);
                        let z = Complex::from_polar(T::from_f64(r), T::from_f64(th));
                        let nu_t = T::from_f64(nu);
                        let (Ok(u), Ok(s)) = (unscaled(nu_t, z), scaled(nu_t, z)) else {
                            continue;
                        };
                        let (u, s) = (to64(u), to64(s));
                        let ok = |v: Complex64| v.norm().is_finite() && v.norm() > floor;
                        if !ok(u) || !ok(s) {
                            continue;
                        }
                        let expected = (u.ln() + log_factor(name, to64(z))).exp();
                        let err = (s - expected).norm() / expected.norm();
                        assert!(
                            err < tol,
                            "{name}: nu={nu}, z={z:?}: scaled={s:?}, exp(c)*unscaled={expected:?}, \
                             rel_err={err:e}"
                        );
                        compared += 1;
                    }
                }
            }
        }
        assert!(compared > 0);
    }

    #[test]
    fn scaled_matches_unscaled_f64() {
        // Radii past ELIM/2 ≈ 350 with orders up to beyond |z|. Observed up to
        // ~6e-13 (large ν and |z|); a lost or mis-scaled term gives O(1).
        check::<f64>(
            &[0.0, 0.3, 2.5, 45.7, 250.0, 575.0, 800.0, 1200.0],
            &[0.5, 5.0, 40.0, 200.0, 360.0, 450.0, 600.0],
            1e-290,
            1e-10,
        );
    }

    #[test]
    fn scaled_matches_unscaled_f32() {
        // ELIM/2 ≈ 40 in f32. Observed up to ~4e-5.
        check::<f32>(
            &[0.0, 0.3, 2.5, 25.7, 65.0, 90.0, 130.0],
            &[0.5, 5.0, 20.0, 41.0, 50.0, 60.0, 75.0],
            1e-30,
            1e-3,
        );
    }
}

// ── Tests for _seq functions (negative orders, recurrence paths) ──

#[cfg(test)]
mod seq_tests {
    use super::*;
    use crate::utils::SEQ_WORKSPACE;
    use num_complex::Complex64;

    /// Longest sequence used below; buffers are fixed arrays so that these
    /// tests also run without `std`.
    const MAX_N: usize = 128;

    fn zeros<T: BesselFloat>() -> [Complex<T>; MAX_N] {
        [Complex::new(T::zero(), T::zero()); MAX_N]
    }

    /// Power series path (ZSERI) where the top orders underflow after the
    /// near-underflow scaling was switched on: the lower orders retried after
    /// that must keep the scaling consistent (they came out multiplied by TOL,
    /// and so did every order below them through the recurrence).
    /// Whether a top order underflows is a rounding-level ZUCHK decision on
    /// the real axis, so scan several top orders and require at least one hit.
    #[test]
    fn seq_power_series_underflow_retry_keeps_scaling() {
        fn check<T: BesselFloat>(
            tops: core::ops::Range<i32>,
            ns: &[usize],
            z: Complex<T>,
            tol: f64,
        ) {
            for &n in ns {
                let mut hits = 0;
                for top in tops.clone() {
                    let nu = T::from_f64(f64::from(top) - (n - 1) as f64);
                    let mut buf = zeros::<T>();
                    let values = &mut buf[..n];
                    let info = besselj_seq_into(nu, z, Scaling::Unscaled, values).unwrap();
                    if info.underflow_count == 0 || values[0].norm() == T::zero() {
                        continue;
                    }
                    hits += 1;
                    for (j, &got) in values.iter().enumerate() {
                        if got.norm() == T::zero() {
                            continue;
                        }
                        let order = nu + T::from_f64(j as f64);
                        let single = besselj(order, z).unwrap();
                        let err = ((got - single).norm() / single.norm()).to_f64().unwrap();
                        // Loose bound: the scaled path near underflow differs from the
                        // single value by ~1e-13 (as in Fortran); the bug was a factor of TOL.
                        assert!(
                            err < tol,
                            "n={n} order={:?}: {got:?} vs {single:?}",
                            order.to_f64()
                        );
                    }
                }
                assert!(hits > 0, "n={n}: no case with a partial underflow");
            }
        }
        check(202..219, &[3, 40], Complex64::new(6.5, 0.0), 1e-10);
        check(36..48, &[3, 20], Complex::new(6.5f32, 0.0), 1e-4);
    }

    /// With ν < 0, underflow can be at either end of the sequence: at large
    /// |orders| below zero or at the high orders above it. `underflow_count`
    /// must count every zero, as Amos does for ν ≥ 0 (including z = 0).
    #[test]
    fn seq_neg_underflow_count_counts_all_zeros() {
        fn check<T: BesselFloat>(radii: &[f64]) {
            let mut hits = 0;
            for nu in [-0.5, -3.0, -40.3, -100.0, -127.0] {
                for &r in radii {
                    for deg in [0.0, 60.0, 135.0, -100.0] {
                        let a = f64::to_radians(deg);
                        let z = Complex::new(T::from_f64(r * a.cos()), T::from_f64(r * a.sin()));
                        for s in [Scaling::Unscaled, Scaling::Exponential] {
                            for (name, seq, _, _) in seq_fns::<T>() {
                                let mut out = zeros::<T>();
                                let Ok(info) = seq(T::from_f64(nu), z, s, &mut out) else {
                                    continue;
                                };
                                let count = out.iter().filter(|v| v.norm() == T::zero()).count();
                                assert_eq!(
                                    info.underflow_count, count,
                                    "{name} ν={nu} r={r} {deg}°"
                                );
                                hits += usize::from(count > 0);
                            }
                        }
                    }
                }
            }
            assert!(hits > 0, "no case with underflow");
        }
        check::<f64>(&[0.1, 6.5, 40.0, 700.0]);
        check::<f32>(&[0.1, 6.5, 40.0, 90.0]);

        // Trailing underflow after the orders cross zero (was reported as 0)
        let z = Complex64::new(0.1, 0.0);
        let mut out = zeros::<f64>();
        let info = besselj_seq_into(-0.5, z, Scaling::Unscaled, &mut out).unwrap();
        assert!(out[0].norm() > 0.0 && out[MAX_N - 1].norm() == 0.0);
        assert!(info.underflow_count > 0);

        // z = 0: every order but 0 is zero, counted as for ν ≥ 0
        let z = Complex64::new(0.0, 0.0);
        for nu in [-2.0, 2.0] {
            for seq in [besselj_seq_into::<f64>, besseli_seq_into::<f64>] {
                let mut out = [Complex64::new(0.0, 0.0); 5];
                let info = seq(nu, z, Scaling::Unscaled, &mut out).unwrap();
                let expected = if nu < 0.0 { 4 } else { 5 };
                assert_eq!(info.underflow_count, expected, "ν={nu}");
            }
        }
    }

    /// Above the half-precision threshold of |z| (~32767 for f64, ~2048 for
    /// f32), every path must report `Reduced`: positive ν, the negative
    /// lattices, crossings, and Y's chunked K values (n > SEQ_WORKSPACE).
    #[test]
    fn seq_reports_reduced_accuracy() {
        fn check<T: BesselFloat>(r: f64) {
            let z = Complex::new(T::from_f64(r * 0.6), T::from_f64(r * 0.8));
            for nu in [-80.5, -3.0, -2.5, 0.5] {
                for n in [5, 70] {
                    for (name, seq, _, _) in seq_fns::<T>() {
                        let mut out = zeros::<T>();
                        let info =
                            seq(T::from_f64(nu), z, Scaling::Exponential, &mut out[..n]).unwrap();
                        assert_eq!(info.status, Accuracy::Reduced, "{name} ν={nu} n={n}");
                    }
                }
            }
        }
        check::<f64>(4.0e4);
        check::<f32>(3.0e3);
    }

    // ── Edge: empty output ──

    #[test]
    fn seq_into_empty_output_is_invalid_input() {
        let z = Complex64::new(1.5, 0.5);
        let s = Scaling::Unscaled;
        // Non-integer and integer ν take different code paths
        for nu in [-1.5, -2.0, 1.5] {
            for (name, seq, _, _) in seq_fns::<f64>() {
                assert_eq!(seq(nu, z, s, &mut []), Err(Error::InvalidInput), "{name}");
            }
        }
    }

    // ── Seq/single consistency over a grid of regions ──
    //
    // Every element of a `_seq` result must match the single-value function at
    // the same order. The grid spans negative and positive orders, small to
    // near-overflow |z|, all four quadrants, and long sequences, so it reaches
    // the recurrence, rescaling, and analytic-continuation paths that only run
    // for n > 1.

    type SeqIntoFn<T> = fn(T, Complex<T>, Scaling, &mut [Complex<T>]) -> Result<SeqInfo, Error>;
    type SingleFn<T> = fn(T, Complex<T>) -> Result<Complex<T>, Error>;
    /// (name, seq, single unscaled, single scaled)
    type SeqCase<T> = (&'static str, SeqIntoFn<T>, SingleFn<T>, SingleFn<T>);

    fn seq_fns<T: BesselFloat>() -> [SeqCase<T>; 6] {
        [
            ("J", besselj_seq_into, besselj, besselj_scaled),
            ("Y", bessely_seq_into, bessely, bessely_scaled),
            ("I", besseli_seq_into, besseli, besseli_scaled),
            ("K", besselk_seq_into, besselk, besselk_scaled),
            ("H1", hankel1_seq_into, hankel1, hankel1_scaled),
            ("H2", hankel2_seq_into, hankel2, hankel2_scaled),
        ]
    }

    /// Compare seq[j] with single(ν+j) for every (ν, |z|, arg z) on the grid.
    /// A sequence may fail as a whole only if some of its orders fail on their
    /// own (e.g. orders below the requested range must not cause an error),
    /// with one exception kept from the Fortran: K tests order max(last
    /// order, n), up to one order above the last for orders starting below 1,
    /// so such K, H, and Y sequences can report `Overflow` just below the
    /// limit (e.g. K from ν = 0.5 at z = 0.05 with n = 97). The converse does
    /// not hold either: Amos's overflow tests depend on the orders and length
    /// of the call, so a single value can report `Overflow` where a sequence
    /// that contains it succeeds (as in the Fortran, e.g. f32 K at ν = 10.5
    /// alone vs. from ν = 7.5 with n = 9, z = 90·e^{-2.6i}). Elements whose
    /// single value is an error, non-finite, or below 1e3 times the smallest
    /// normal number are skipped.
    fn check_seq_grid<T: BesselFloat>(
        nus: &[f64],
        radii: &[f64],
        args: &[f64],
        n: usize,
        tol: f64,
    ) {
        let floor = T::MACH_TINY * T::from_f64(1e3);
        for (name, seq, single, single_scaled) in seq_fns::<T>() {
            for scaling in [Scaling::Unscaled, Scaling::Exponential] {
                let single = if scaling == Scaling::Unscaled {
                    single
                } else {
                    single_scaled
                };
                for &nu in nus {
                    for &r in radii {
                        for &th in args {
                            let nu_t = T::from_f64(nu);
                            let z = Complex::from_polar(T::from_f64(r), T::from_f64(th));
                            let mut buf = zeros::<T>();
                            let values = &mut buf[..n];
                            match seq(nu_t, z, scaling, values) {
                                Ok(_) => {}
                                Err(e) => {
                                    let some_order_fails = (0..n)
                                        .any(|j| single(nu_t + T::from_f64(j as f64), z).is_err());
                                    assert!(
                                        some_order_fails,
                                        "{name} {scaling:?}: nu={nu}, z={z:?}, n={n}: seq \
                                         failed with {e:?} but every single value succeeds"
                                    );
                                    continue;
                                }
                            }
                            for (j, &w) in values.iter().enumerate() {
                                let order = nu_t + T::from_f64(j as f64);
                                let Ok(v) = single(order, z) else {
                                    continue;
                                };
                                if !(v.norm() > floor && v.norm().is_finite()) {
                                    continue;
                                }
                                let err = ((w - v).norm() / v.norm()).to_f64().unwrap();
                                assert!(
                                    err < tol,
                                    "{name} {scaling:?}: nu={nu}, z={z:?}, n={n}, j={j}: \
                                     seq={w:?}, single={v:?}, rel_err={err:e}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    /// Angles covering all four quadrants, the real and imaginary axes, and
    /// both sides of the |Im z| = √3·|Re z| region boundary.
    const ARGS: [f64; 8] = [
        0.0,
        0.5,
        1.2,
        core::f64::consts::FRAC_PI_2,
        2.2,
        3.0,
        -1.0,
        -2.6,
    ];

    /// Recurrence paths in f64. The seq values come from recurrence while the
    /// single values come from direct expansions, so their ~1e-14 errors do not
    /// cancel (observed up to ~6e-12 for n = 100); a misplaced or mis-scaled
    /// element gives O(1e-3) error or worse.
    const TOL_GRID_F64: f64 = 1e-10;

    #[test]
    fn seq_matches_single_grid_f64() {
        // Orders: all-negative / crossing (integer and not), small, and beyond
        // FNUL (≈ 86) so the uniform expansions recur. Radii include values near
        // the ELIM/ALIM scaling thresholds (≈ 665–700) and beyond, where low
        // orders overflow but ν = -200 does not. ν = 575 at |z| = 360 reaches
        // scaled Y with 2|Im z| > ELIM, where the K term still dominates. -0.0
        // is an edge case of the sign dispatch; n = 1 is the one-element
        // lattice edge case. For n = 9, ν = -3 / -4 / -7 cover an integer
        // crossing whose largest |order| is positive / both / negative.
        for n in [1, 9] {
            check_seq_grid::<f64>(
                &[
                    -200.0, -20.0, -12.3, -7.0, -5.3, -4.0, -3.0, -2.5, -0.0, 0.0, 0.4, 13.5, 46.2,
                    81.7, 121.3, 575.0,
                ],
                &[
                    0.3, 3.0, 17.0, 45.0, 120.0, 330.0, 360.0, 666.0, 698.0, 705.0, 709.0,
                ],
                &ARGS,
                n,
                TOL_GRID_F64,
            );
        }
    }

    #[test]
    fn seq_matches_single_grid_f64_long() {
        // n = 100 with small ν: ν + n − 1 exceeds FNUL while ν itself does not.
        // Negative non-integer ν spans several stack chunks of the secondary
        // function (all negative, and crossing zero with more negative orders
        // than a chunk holds).
        check_seq_grid::<f64>(
            &[-80.3, -40.3, 0.5, 3.3],
            &[2.2, 25.0, 90.0],
            &ARGS,
            100,
            TOL_GRID_F64,
        );
    }

    /// f32 hits the scaling thresholds at |values| ≈ 1e28, so moderate |z|
    /// exercises the rescaling branches. Observed up to ~1.4e-4 on this grid;
    /// misplaced or mis-scaled elements give O(1) error or worse. (Off this grid,
    /// near the turning point ν ≈ |z| ≈ FNUL, Amos's single-precision algorithm
    /// itself loses up to ~2e-2, so random points would need a looser bound.)
    const TOL_GRID_F32: f64 = 1e-3;

    #[test]
    fn seq_matches_single_grid_f32() {
        // At |z| = 90, low orders overflow but ν = -60 does not. ν = 70 at
        // |z| = 45 reaches scaled Y with 2|Im z| > ELIM (K term dominant).
        check_seq_grid::<f32>(
            &[
                -60.0, -9.0, -4.6, -0.0, 0.0, 0.4, 7.5, 21.2, 29.9, 46.3, 70.0,
            ],
            &[0.3, 3.0, 17.0, 45.0, 70.0, 90.0],
            &ARGS,
            9,
            TOL_GRID_F32,
        );
    }

    // ── Three-term recurrence across ν = 0 ──
    //
    // Negative orders come from reflection formulas shared by the single-value
    // and `_seq` functions, so comparing the two cannot catch an error in a
    // formula. The recurrence can: across ν = 0 it links reflected values to
    // directly computed ones. Every scaling factor is independent of the order,
    // so scaled values satisfy the same recurrence.

    /// Check C(ν-1) + s·C(ν+1) = t·(2ν/z)·C(ν) on consecutive `_seq` elements,
    /// with (s, t) = (1, 1) for J, Y, H1, H2, (-1, 1) for I, and (-1, -1) for K
    /// (DLMF 10.6.1, 10.29.1). The residual is relative to the largest term.
    fn check_recurrence_grid<T: BesselFloat>(
        nus: &[f64],
        radii: &[f64],
        args: &[f64],
        n: usize,
        tol: f64,
    ) {
        for (name, seq, _, _) in seq_fns::<T>() {
            let (s, t) = match name {
                "I" => (-1.0, 1.0),
                "K" => (-1.0, -1.0),
                _ => (1.0, 1.0),
            };
            let (s, t) = (T::from_f64(s), T::from_f64(t));
            for scaling in [Scaling::Unscaled, Scaling::Exponential] {
                for &nu in nus {
                    for &r in radii {
                        for &th in args {
                            let nu_t = T::from_f64(nu);
                            let z = Complex::from_polar(T::from_f64(r), T::from_f64(th));
                            let mut buf = zeros::<T>();
                            let values = &mut buf[..n];
                            if seq(nu_t, z, scaling, values).is_err() {
                                continue;
                            }
                            for (j, w) in values.windows(3).enumerate() {
                                // Underflowed (zeroed) or overflowed elements break the recurrence
                                if w.iter()
                                    .any(|v| v.norm() == T::zero() || !v.norm().is_finite())
                                {
                                    continue;
                                }
                                let order = nu_t + T::from_f64((j + 1) as f64);
                                let lhs = w[0] + w[2] * s;
                                let rhs = w[1] * (z.inv() * (T::from_f64(2.0) * order * t));
                                let mag = w[0].norm().max(w[2].norm()).max(rhs.norm());
                                let err = ((lhs - rhs).norm() / mag).to_f64().unwrap();
                                assert!(
                                    err < tol,
                                    "{name} {scaling:?}: nu={nu}, z={z:?}, order={:?}: \
                                     recurrence residual {err:e}",
                                    order.to_f64()
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn seq_recurrence_across_zero_f64() {
        // Non-integer (generic and half-integer) and integer orders crossing
        // zero, plus all-negative sequences.
        check_recurrence_grid::<f64>(
            &[-12.3, -5.3, -4.0, -3.7, -2.5, -0.3],
            &[0.3, 3.0, 17.0, 45.0],
            &ARGS,
            9,
            1e-11, // observed up to ~1.3e-13
        );
    }

    #[test]
    fn seq_recurrence_across_zero_f32() {
        check_recurrence_grid::<f32>(
            &[-12.3, -5.3, -4.0, -3.7, -2.5, -0.3],
            &[0.3, 3.0, 17.0, 45.0],
            &ARGS,
            9,
            1e-4, // observed up to ~1e-5
        );
    }

    // ── Chunked workspace ──
    //
    // Y's K values and the secondary function of the negative non-integer
    // reflection formulas go to a workspace, in chunks of its length (Y, and
    // I with non-integer ν < 0) or of half of it (J and Y with non-integer
    // ν < 0). A workspace that covers the whole
    // sequence (as `_seq` allocates) gives one Amos call per function, so
    // compare short workspaces against it. Amos's tests depend on the orders
    // and length of each call, so a chunked sequence may fail only where some
    // single value fails, as above.

    type ImplFn<T> = for<'a> fn(
        T,
        Complex<T>,
        Scaling,
        &mut [Complex<T>],
        Workspace<'a, T>,
    ) -> Result<SeqInfo, Error>;

    /// (name, implementation, single unscaled, single scaled)
    type ImplCase<T> = (&'static str, ImplFn<T>, SingleFn<T>, SingleFn<T>);

    fn impl_fns<T: BesselFloat>() -> [ImplCase<T>; 3] {
        [
            ("J", besselj_seq_impl, besselj, besselj_scaled),
            ("Y", bessely_seq_impl, bessely, bessely_scaled),
            ("I", besseli_seq_impl, besseli, besseli_scaled),
        ]
    }

    /// Bitwise equality (exact through the lossless conversion to f64).
    fn same_bits<T: BesselFloat>(x: Complex<T>, y: Complex<T>) -> bool {
        let bits = |v: T| v.to_f64().unwrap().to_bits();
        bits(x.re) == bits(y.re) && bits(x.im) == bits(y.im)
    }

    fn check_chunk_invariance<T: BesselFloat>(
        nus: &[f64],
        radii: &[f64],
        ns: &[usize],
        ws_len: usize,
        tol: f64,
    ) {
        let floor = T::MACH_TINY * T::from_f64(1e3);
        let czero = Complex::new(T::zero(), T::zero());
        for (name, f, single, single_scaled) in impl_fns::<T>() {
            for scaling in [Scaling::Unscaled, Scaling::Exponential] {
                let single = if scaling == Scaling::Unscaled {
                    single
                } else {
                    single_scaled
                };
                for &nu in nus {
                    for &r in radii {
                        for &th in &ARGS {
                            for &n in ns {
                                let nu_t = T::from_f64(nu);
                                let z = Complex::from_polar(T::from_f64(r), T::from_f64(th));
                                let (mut a, mut b) = (zeros::<T>(), zeros::<T>());
                                let (a, b) = (&mut a[..n], &mut b[..n]);
                                let mut ws_whole = [czero; 2 * MAX_N];
                                let mut ws_short = [czero; 2 * MAX_N];
                                let whole = f(nu_t, z, scaling, a, Workspace::Slice(&mut ws_whole));
                                let chunked = f(
                                    nu_t,
                                    z,
                                    scaling,
                                    b,
                                    Workspace::Slice(&mut ws_short[..ws_len]),
                                );
                                let some_single_fails = || {
                                    (0..n).any(|j| single(nu_t + T::from_f64(j as f64), z).is_err())
                                };
                                let (info_a, info_b) = match (whole, chunked) {
                                    (Ok(ia), Ok(ib)) => (ia, ib),
                                    (Ok(_), Err(e)) => {
                                        assert!(
                                            some_single_fails(),
                                            "{name} {scaling:?} ws={ws_len}: nu={nu}, z={z:?}, \
                                             n={n}: chunked failed with {e:?}, whole call and \
                                             every single value succeed"
                                        );
                                        continue;
                                    }
                                    // Chunks can pass Amos's tests where the whole call
                                    // fails (e.g. K's test of max(last order, n)); the
                                    // values must then match the single values.
                                    (Err(_), Ok(_)) => {
                                        for (j, &y) in b.iter().enumerate() {
                                            let Ok(v) = single(nu_t + T::from_f64(j as f64), z)
                                            else {
                                                continue;
                                            };
                                            if !(v.norm() > floor && v.norm().is_finite()) {
                                                continue;
                                            }
                                            let err = ((y - v).norm() / v.norm()).to_f64().unwrap();
                                            assert!(
                                                err < tol,
                                                "{name} {scaling:?} ws={ws_len}: nu={nu}, \
                                                 z={z:?}, n={n}, j={j}: chunked={y:?}, \
                                                 single={v:?}"
                                            );
                                        }
                                        continue;
                                    }
                                    (Err(_), Err(_)) => continue,
                                };
                                if 2 * n <= ws_len {
                                    // One chunk: the same calls as the whole sequence
                                    assert_eq!(
                                        info_a, info_b,
                                        "{name} {scaling:?} ws={ws_len}: nu={nu}, n={n}"
                                    );
                                    for (&x, &y) in a.iter().zip(b.iter()) {
                                        assert!(
                                            same_bits(x, y),
                                            "{name} {scaling:?} ws={ws_len}: nu={nu}, z={z:?}, \
                                             n={n}: {x:?} vs {y:?}"
                                        );
                                    }
                                    continue;
                                }
                                // The underflow count is decided per call, so only the
                                // status must match
                                assert_eq!(
                                    info_a.status, info_b.status,
                                    "{name} {scaling:?} ws={ws_len}: nu={nu}, n={n}"
                                );
                                for (j, (&x, &y)) in a.iter().zip(b.iter()).enumerate() {
                                    if !(x.norm() > floor && x.norm().is_finite()) {
                                        continue;
                                    }
                                    let err = ((x - y).norm() / x.norm()).to_f64().unwrap();
                                    assert!(
                                        err < tol,
                                        "{name} {scaling:?} ws={ws_len}: nu={nu}, z={z:?}, \
                                         n={n}, j={j}: whole={x:?}, chunked={y:?}, \
                                         rel_err={err:e}"
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// Workspace lengths: chunks of 1, 2, 7 orders for J and Y with non-integer
    /// ν < 0 (half the workspace) and of 2, 4, 14 otherwise, and the stack
    /// workspace.
    const WS_LENS: [usize; 4] = [2, 4, 14, SEQ_WORKSPACE];

    #[test]
    fn seq_chunk_size_invariance_f64() {
        // Negative non-integer ν (all negative, and crossing zero) chunk the
        // secondary function; positive ν and negative integer ν (both mirror
        // branches of a crossing: n = 100 from −40 and from −70) chunk the K
        // values inside Y. Radii reach the ELIM/ALIM thresholds. With the
        // stack workspace, n = 32 and 33 end at a chunk boundary or one past it
        // for J and Y with non-integer ν < 0, and n = 9 is a single chunk.
        let nus = [-80.3, -70.0, -40.3, -40.0, -12.3, -0.7, 0.4, 3.3, 81.7];
        let radii = [0.3, 3.0, 17.0, 45.0, 120.0, 666.0, 698.0, 705.0];
        let ns = [9, 32, 33, 65, 100];
        for ws_len in WS_LENS {
            check_chunk_invariance::<f64>(&nus, &radii, &ns, ws_len, TOL_CHUNK_F64);
        }
    }

    #[test]
    fn seq_chunk_size_invariance_f32() {
        // Thresholds sit at |z| ≈ 34–100 in f32; z = 90·e^{-2.6i} is where K's
        // overflow test depends on the call (see `check_seq_grid`).
        let nus = [-70.0, -40.3, -40.0, -15.5, -4.6, 0.4, 7.5, 29.9];
        let radii = [3.0, 17.0, 34.0, 45.0, 70.0, 90.0, 100.0];
        let ns = [9, 32, 33, 65, 100];
        for ws_len in WS_LENS {
            check_chunk_invariance::<f32>(&nus, &radii, &ns, ws_len, TOL_CHUNK_F32);
        }
    }

    /// Chunks restart Amos's algorithms, so their values differ from the whole
    /// call as single values do, plus where Amos's own error is largest: on
    /// these grids up to ~1.4e-9 in f64 (Y near |z| = 705 on the real axis,
    /// where both differ from mpmath by ~3e-9) and ~3.9e-3 in f32 (near the
    /// turning point). Off the grids, cancellation in the reflection formulas
    /// can reach ~1.3e-8 (J, ν = −250.7, z = 698, where |J| ≈ 5e-7 comes from
    /// J and Y terms ≈ 3e-2); the whole call itself shifts by ~3e-9 there with
    /// n. A misplaced or mis-scaled chunk gives O(1) error.
    const TOL_CHUNK_F64: f64 = 1e-8;
    const TOL_CHUNK_F32: f64 = 2e-2;

    // ── Output buffer contents ──

    #[test]
    fn seq_into_ignores_initial_out_contents() {
        // Every path must overwrite `out`; NaN left behind would show up
        let z_list = [Complex64::new(1.5, 0.5), Complex64::new(-30.0, 20.0)];
        for (name, seq, _, _) in seq_fns::<f64>() {
            for nu in [-40.3, -7.0, -2.5, 0.4, 3.0] {
                for z in z_list {
                    for n in [1, 9, 40] {
                        let mut a = [Complex64::new(0.0, 0.0); 40];
                        let mut b = [Complex64::new(f64::NAN, f64::NAN); 40];
                        for scaling in [Scaling::Unscaled, Scaling::Exponential] {
                            let ra = seq(nu, z, scaling, &mut a[..n]);
                            let rb = seq(nu, z, scaling, &mut b[..n]);
                            assert_eq!(ra, rb, "{name}: nu={nu}, z={z:?}, n={n}");
                            if ra.is_err() {
                                continue;
                            }
                            for (&x, &y) in a[..n].iter().zip(b[..n].iter()) {
                                assert!(
                                    same_bits(x, y),
                                    "{name}: nu={nu}, z={z:?}, n={n}: {x:?} vs {y:?}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}

// ── Allocating `_seq` wrappers ──

#[cfg(all(test, feature = "std"))]
mod seq_vec_tests {
    use super::*;
    use num_complex::Complex64;

    /// `_seq` gives every function a heap workspace of the full length, so each
    /// is one Amos call (as in 0.4), unlike the stack chunks of `_seq_into`.
    /// These long sequences span several stack chunks.
    #[test]
    fn seq_computes_each_function_in_one_call() {
        type Seq = fn(f64, Complex64, usize, Scaling) -> Result<BesselResult<f64>, Error>;
        type Impl = for<'a> fn(
            f64,
            Complex64,
            Scaling,
            &mut [Complex64],
            Workspace<'a, f64>,
        ) -> Result<SeqInfo, Error>;
        let cases: [(Seq, Impl, f64, Complex64); 4] = [
            (
                bessely_seq,
                bessely_seq_impl,
                0.4,
                Complex64::from_polar(25.0, 0.5),
            ),
            (
                bessely_seq,
                bessely_seq_impl,
                -40.3,
                Complex64::from_polar(17.0, -1.0),
            ),
            (
                besselj_seq,
                besselj_seq_impl,
                -80.3,
                Complex64::from_polar(45.0, 3.0),
            ),
            (
                besseli_seq,
                besseli_seq_impl,
                -40.3,
                Complex64::from_polar(17.0, 2.2),
            ),
        ];
        const N: usize = 100;
        for (seq, imp, nu, z) in cases {
            let r = seq(nu, z, N, Scaling::Unscaled).unwrap();
            let mut out = [Complex64::new(0.0, 0.0); N];
            let mut ws = [Complex64::new(0.0, 0.0); 2 * N];
            let info = imp(
                nu,
                z,
                Scaling::Unscaled,
                &mut out,
                Workspace::Slice(&mut ws),
            )
            .unwrap();
            assert_eq!(
                (r.underflow_count, r.status),
                (info.underflow_count, info.status)
            );
            for (j, (x, y)) in r.values.iter().zip(out.iter()).enumerate() {
                assert!(
                    x.re.to_bits() == y.re.to_bits() && x.im.to_bits() == y.im.to_bits(),
                    "nu={nu}, z={z:?}, j={j}: {x:?} vs {y:?}"
                );
            }
        }
    }

    #[test]
    fn seq_matches_seq_into() {
        type Pair = (
            fn(f64, Complex64, usize, Scaling) -> Result<BesselResult<f64>, Error>,
            fn(f64, Complex64, Scaling, &mut [Complex64]) -> Result<SeqInfo, Error>,
        );
        let pairs: [Pair; 6] = [
            (besselj_seq, besselj_seq_into),
            (bessely_seq, bessely_seq_into),
            (besseli_seq, besseli_seq_into),
            (besselk_seq, besselk_seq_into),
            (hankel1_seq, hankel1_seq_into),
            (hankel2_seq, hankel2_seq_into),
        ];
        let z = Complex64::new(1.5, -0.5);
        for (seq, seq_into) in pairs {
            for nu in [-5.3, -2.0, 0.0, 1.7] {
                for n in [0, 1, 6] {
                    let mut out = [Complex64::new(0.0, 0.0); 6];
                    let expected = seq_into(nu, z, Scaling::Unscaled, &mut out[..n]);
                    match seq(nu, z, n, Scaling::Unscaled) {
                        Ok(r) => {
                            let info = expected.unwrap();
                            assert_eq!(r.values, &out[..n]);
                            assert_eq!(r.underflow_count, info.underflow_count);
                            assert_eq!(r.status, info.status);
                        }
                        Err(e) => assert_eq!(Err(e), expected),
                    }
                }
            }
        }
    }
}
