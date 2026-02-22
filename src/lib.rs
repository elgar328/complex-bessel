//! Pure Rust implementation of complex Bessel functions based on Amos Algorithm 644 (ACM TOMS 644).
//!
//! # Features
//!
//! - **f32 & f64** — all functions accept `Complex<f64>` or `Complex<f32>`
//! - **Complete function set** — J, Y, I, K, H<sup>(1)</sup>, H<sup>(2)</sup>, Ai, Bi
//! - **Consecutive orders** — `_seq` variants return ν, ν+1, …, ν+n−1 in one call
//! - **Exponential scaling** — `_scaled` variants prevent overflow/underflow
//! - **Negative orders** — supports ν < 0 via DLMF reflection formulas (not in Amos)
//! - **`no_std` support** — 3-tier: bare `no_std` (no allocator), `alloc`, `std` (default)
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
//! | [`besselj`]`(ν, z)`<br>[`besselj_scaled`]`(ν, z)`<br>[`besselj_seq`]`(ν, z, n, scaling)` | J<sub>ν</sub>(z), Bessel first kind | exp(−\|Im(z)\|) · J<sub>ν</sub>(z) |
//! | [`bessely`]`(ν, z)`<br>[`bessely_scaled`]`(ν, z)`<br>[`bessely_seq`]`(ν, z, n, scaling)` | Y<sub>ν</sub>(z), Bessel second kind | exp(−\|Im(z)\|) · Y<sub>ν</sub>(z) |
//! | [`besseli`]`(ν, z)`<br>[`besseli_scaled`]`(ν, z)`<br>[`besseli_seq`]`(ν, z, n, scaling)` | I<sub>ν</sub>(z), modified first kind | exp(−\|Re(z)\|) · I<sub>ν</sub>(z) |
//! | [`besselk`]`(ν, z)`<br>[`besselk_scaled`]`(ν, z)`<br>[`besselk_seq`]`(ν, z, n, scaling)` | K<sub>ν</sub>(z), modified second kind | exp(z) · K<sub>ν</sub>(z) |
//! | [`hankel1`]`(ν, z)`<br>[`hankel1_scaled`]`(ν, z)`<br>[`hankel1_seq`]`(ν, z, n, scaling)` | H<sub>ν</sub><sup>(1)</sup>(z), Hankel first kind | exp(−iz) · H<sub>ν</sub><sup>(1)</sup>(z) |
//! | [`hankel2`]`(ν, z)`<br>[`hankel2_scaled`]`(ν, z)`<br>[`hankel2_seq`]`(ν, z, n, scaling)` | H<sub>ν</sub><sup>(2)</sup>(z), Hankel second kind | exp(iz) · H<sub>ν</sub><sup>(2)</sup>(z) |
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
//! All functions (single-value and `_seq` variants) accept any real order, including negative values.
//! DLMF reflection formulas are applied automatically:
//!
//! - **J**: J<sub>−ν</sub>(z) = cos(νπ) J<sub>ν</sub>(z) − sin(νπ) Y<sub>ν</sub>(z) (DLMF 10.2.3)
//! - **Y**: Y<sub>−ν</sub>(z) = sin(νπ) J<sub>ν</sub>(z) + cos(νπ) Y<sub>ν</sub>(z) (DLMF 10.2.3)
//! - **I**: I<sub>−ν</sub>(z) = I<sub>ν</sub>(z) + (2/π) sin(νπ) K<sub>ν</sub>(z) (DLMF 10.27.2)
//! - **K**: K<sub>−ν</sub>(z) = K<sub>ν</sub>(z) (even in ν, DLMF 10.27.3)
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
//! The `_raw` Airy variants ([`airy_raw`], [`biry_raw`], …) similarly return
//! an [`AiryResult`] with [`Accuracy`].
//! All other functions return only the computed value without [`Accuracy`].
//!
//! **[`Accuracy`]:**
//!
//! | Status | Meaning |
//! |--------|---------|
//! | [`Normal`](Accuracy::Normal) | Full machine precision |
//! | [`Reduced`](Accuracy::Reduced) | More than half of significant digits may be lost.<br>Occurs only when \|z\| or ν exceeds ~32767 |
//!
//! [`Reduced`](Accuracy::Reduced) is extremely rare in practice. SciPy's Bessel wrappers also silently
//! discard the equivalent Amos IERR=3 flag by default.
//!
//! To check accuracy status, use a `_seq` function (Bessel) or a `_raw` function (Airy):
//!
//! ```
//! # #[cfg(feature = "alloc")] {
//! use complex_bessel::*;
//! use num_complex::Complex;
//!
//! let z = Complex::new(1.0, 2.0);
//! let result = besselk_seq(0.0, z, 1, Scaling::Unscaled).unwrap();
//! assert!(matches!(result.status, Accuracy::Normal));
//!
//! let result = airy_raw(z, Scaling::Unscaled).unwrap();
//! assert!(matches!(result.status, Accuracy::Normal));
//! # }
//! ```
//!
//! # Error handling
//!
//! All functions return `Result<_, Error>`. The four error variants are:
//!
//! | Variant | Cause |
//! |---------|-------|
//! | [`InvalidInput`](Error::InvalidInput) | z = 0 for K/Y/H, n < 1 |
//! | [`Overflow`](Error::Overflow) | \|z\| or ν too large (or too small) for finite result |
//! | [`TotalPrecisionLoss`](Error::TotalPrecisionLoss) | Complete loss of significant digits; \|z\| or ν too large |
//! | [`ConvergenceFailure`](Error::ConvergenceFailure) | Internal algorithm did not converge |
//!
//! [`Error`] implements `Display` always and `std::error::Error` with the `std` feature.
//!
//! # `no_std` support
//!
//! | Cargo features | Available API |
//! |---------------|---------------|
//! | `default-features = false` | 24 single-value functions |
//! | `features = ["alloc"]` | + 6 `_seq` variants + [`BesselResult`] |
//! | `features = ["std"]` (default) | + `impl Error for Error` |
//!
//! The 24 single-value functions include 12 Bessel (J/Y/I/K/H<sup>(1)</sup>/H<sup>(2)</sup> × unscaled/scaled)
//! and 12 Airy (Ai/Ai'/Bi/Bi' × unscaled/scaled/raw).
//!
//! ```toml
//! # Bare no_std — no allocator needed:
//! complex-bessel = { version = "0.1", default-features = false }
//!
//! # no_std with alloc (adds _seq functions and BesselResult):
//! complex-bessel = { version = "0.1", default-features = false, features = ["alloc"] }
//! ```

#![warn(missing_docs)]
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
pub use types::{Accuracy, AiryResult, Error, Scaling};

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
    i_val + k_val * (two / pi * utils::sinpi(order))
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

/// Bessel function of the first kind, J_ν(z).
///
/// Computes a single value of the Bessel function J_ν(z) for complex z
/// and real order ν (any real value, including negative).
///
/// For negative ν, the DLMF 10.2.3 reflection formula is applied:
/// `J_{-ν}(z) = cos(νπ) J_ν(z) - sin(νπ) Y_ν(z)`.
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
/// Returns [`Error`] if the computation fails (overflow, precision loss, etc.).
#[inline]
pub fn besselj<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    besselj_internal(nu, z, Scaling::Unscaled)
}

/// Bessel function of the second kind, Y_ν(z).
///
/// Computes a single value of the Bessel function Y_ν(z) for complex z
/// and real order ν (any real value, including negative).
///
/// For negative ν, the DLMF 10.2.3 reflection formula is applied:
/// `Y_{-ν}(z) = sin(νπ) J_ν(z) + cos(νπ) Y_ν(z)`.
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
/// Returns [`Error`] if the computation fails (overflow, z = 0, etc.).
#[inline]
pub fn bessely<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
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
/// Returns [`Error`] if the computation fails (overflow, precision loss, etc.).
#[inline]
pub fn besseli<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    besseli_internal(nu, z, Scaling::Unscaled)
}

/// Modified Bessel function of the second kind, K_ν(z).
///
/// Computes a single value of K_ν(z) for complex z and real order ν
/// (any real value, including negative).
///
/// For negative ν, K_{-ν}(z) = K_ν(z) (K is even in ν, DLMF 10.27.3).
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
/// Returns [`Error`] if the computation fails (overflow, z = 0, etc.).
#[inline]
pub fn besselk<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
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
/// Returns [`Error`] if the computation fails (overflow, z = 0, etc.).
#[inline]
pub fn hankel1<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    hankel_internal(HankelKind::First, nu, z, Scaling::Unscaled)
}

/// Hankel function of the second kind, H_ν^(2)(z).
///
/// Computes a single value of H_ν^(2)(z) for complex z and real order ν
/// (any real value, including negative).
///
/// For negative ν, the DLMF 10.4.6 reflection formula is applied:
/// `H^(2)_{-ν}(z) = exp(-νπi) H^(2)_ν(z)`.
///
/// # Example
///
/// ```
/// use complex_bessel::hankel2;
/// use num_complex::Complex;
///
/// let z = Complex::new(1.0_f64, 0.0);
/// let h = hankel2(0.0, z).unwrap();
/// assert!((h.re - 0.7652).abs() < 1e-3);  // Re = J_0(1) ≈ 0.7652
/// assert!((h.im + 0.0883).abs() < 1e-3);  // Im = -Y_0(1) ≈ -0.0883
/// ```
///
/// # Errors
///
/// Returns [`Error`] if the computation fails (overflow, z = 0, etc.).
#[inline]
pub fn hankel2<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    hankel_internal(HankelKind::Second, nu, z, Scaling::Unscaled)
}

/// Airy function Ai(z).
///
/// Computes the Airy function of the first kind for complex z.
/// Ai(z) is a solution to the differential equation `w'' - z·w = 0`
/// that decays exponentially for large positive real z.
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn airy<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, Error> {
    let (result, _nz, _status) = airy::zairy(z, AiryDerivative::Value, Scaling::Unscaled)?;
    Ok(result)
}

/// Derivative of the Airy function, Ai'(z).
///
/// Computes the derivative of the Airy function of the first kind for complex z.
/// Satisfies the differential equation `Ai''(z) = z · Ai(z)`.
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn airyprime<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, Error> {
    let (result, _nz, _status) = airy::zairy(z, AiryDerivative::Derivative, Scaling::Unscaled)?;
    Ok(result)
}

/// Airy function of the second kind, Bi(z).
///
/// Computes the Airy function Bi(z) for complex z.
/// Bi(z) is the solution to `w'' - z·w = 0` that grows exponentially
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn biry<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, Error> {
    let (result, _status) = airy::zbiry(z, AiryDerivative::Value, Scaling::Unscaled)?;
    Ok(result)
}

/// Derivative of the Airy function of the second kind, Bi'(z).
///
/// Computes the derivative of Bi(z) for complex z.
/// Satisfies the differential equation `Bi''(z) = z · Bi(z)`.
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn biryprime<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, Error> {
    let (result, _status) = airy::zbiry(z, AiryDerivative::Derivative, Scaling::Unscaled)?;
    Ok(result)
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn besselj_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn bessely_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn besseli_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
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
/// let z = Complex::new(500.0_f64, 0.0);
///
/// // K_0(500) would underflow to 0, but the scaled version stays finite:
/// let k_s = besselk_scaled(0.0, z).unwrap();
/// assert!((k_s.re - 0.0560).abs() < 1e-3); // exp(500) * K_0(500) ≈ 0.0560
/// ```
///
/// # Errors
///
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn besselk_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn hankel1_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn hankel2_scaled<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, Error> {
    hankel_internal(HankelKind::Second, nu, z, Scaling::Exponential)
}

/// Scaled Airy function: `exp(ζ) · Ai(z)`, where ζ = (2/3) z√z.
///
/// Ai(z) decays super-exponentially for large positive real z.
/// The scaling factor `exp(ζ)` keeps the result representable.
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn airy_scaled<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, Error> {
    let (result, _nz, _status) = airy::zairy(z, AiryDerivative::Value, Scaling::Exponential)?;
    Ok(result)
}

/// Scaled derivative of the Airy function: `exp(ζ) · Ai'(z)`, where ζ = (2/3) z√z.
///
/// Ai'(z) decays super-exponentially for large positive real z, just as Ai(z) does.
/// Satisfies the differential equation `Ai''(z) = z · Ai(z)`.
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn airyprime_scaled<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, Error> {
    let (result, _nz, _status) = airy::zairy(z, AiryDerivative::Derivative, Scaling::Exponential)?;
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn biry_scaled<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, Error> {
    let (result, _status) = airy::zbiry(z, AiryDerivative::Value, Scaling::Exponential)?;
    Ok(result)
}

/// Scaled derivative of the Airy function of the second kind:
/// `exp(-|Re(ζ)|) · Bi'(z)`, where ζ = (2/3) z√z.
///
/// Bi'(z) grows super-exponentially for large positive real z, just as Bi(z) does.
/// Satisfies the differential equation `Bi''(z) = z · Bi(z)`.
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn biryprime_scaled<T: BesselFloat>(z: Complex<T>) -> Result<Complex<T>, Error> {
    let (result, _status) = airy::zbiry(z, AiryDerivative::Derivative, Scaling::Exponential)?;
    Ok(result)
}

// ── Airy _raw functions (expose Accuracy) ──

/// Airy function Ai(z) with precision status.
///
/// Like [`airy`], but returns an [`AiryResult`] that includes
/// [`Accuracy`] for detecting precision loss at large |z|.
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn airy_raw<T: BesselFloat>(z: Complex<T>, scaling: Scaling) -> Result<AiryResult<T>, Error> {
    let (value, _nz, status) = airy::zairy(z, AiryDerivative::Value, scaling)?;
    Ok(AiryResult { value, status })
}

/// Derivative of the Airy function Ai'(z) with precision status.
///
/// Like [`airyprime`], but returns an [`AiryResult`] that includes
/// [`Accuracy`] for detecting precision loss at large |z|.
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn airyprime_raw<T: BesselFloat>(
    z: Complex<T>,
    scaling: Scaling,
) -> Result<AiryResult<T>, Error> {
    let (value, _nz, status) = airy::zairy(z, AiryDerivative::Derivative, scaling)?;
    Ok(AiryResult { value, status })
}

/// Airy function of the second kind Bi(z) with precision status.
///
/// Like [`biry`], but returns an [`AiryResult`] that includes
/// [`Accuracy`] for detecting precision loss at large |z|.
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn biry_raw<T: BesselFloat>(z: Complex<T>, scaling: Scaling) -> Result<AiryResult<T>, Error> {
    let (value, status) = airy::zbiry(z, AiryDerivative::Value, scaling)?;
    Ok(AiryResult { value, status })
}

/// Derivative of the Airy function of the second kind Bi'(z) with precision status.
///
/// Like [`biryprime`], but returns an [`AiryResult`] that includes
/// [`Accuracy`] for detecting precision loss at large |z|.
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
/// Returns [`Error`] if the computation fails.
#[inline]
pub fn biryprime_raw<T: BesselFloat>(
    z: Complex<T>,
    scaling: Scaling,
) -> Result<AiryResult<T>, Error> {
    let (value, status) = airy::zbiry(z, AiryDerivative::Derivative, scaling)?;
    Ok(AiryResult { value, status })
}

// ── Sequence functions with scaling option (require alloc) ──

#[cfg(feature = "alloc")]
extern crate alloc as alloc_crate;

#[cfg(feature = "alloc")]
fn seq_helper<T: BesselFloat>(
    n: usize,
    f: impl FnOnce(&mut [Complex<T>]) -> Result<(usize, Accuracy), Error>,
) -> Result<BesselResult<T>, Error> {
    let zero = T::zero();
    let mut values = alloc_crate::vec![Complex::new(zero, zero); n];
    let (underflow_count, status) = f(&mut values)?;
    Ok(BesselResult {
        values,
        underflow_count,
        status,
    })
}

/// Take the worse of two statuses (Reduced > Normal).
#[cfg(feature = "alloc")]
#[inline]
fn worse_status(a: Accuracy, b: Accuracy) -> Accuracy {
    match (a, b) {
        (Accuracy::Reduced, _) | (_, Accuracy::Reduced) => Accuracy::Reduced,
        _ => Accuracy::Normal,
    }
}

/// Recount underflow: number of leading zero elements in the output.
#[cfg(feature = "alloc")]
#[inline]
fn recount_underflow<T: BesselFloat>(values: &[Complex<T>]) -> usize {
    let zero = Complex::new(T::zero(), T::zero());
    values.iter().take_while(|v| **v == zero).count()
}

/// Compute the number of negative-order elements (before crossing zero).
/// For orders ν, ν+1, …, ν+n−1, this is the count where ν+j < 0.
#[cfg(feature = "alloc")]
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

// ── Negative-order _seq implementations ──

/// K_seq with negative order support: K_{-ν} = K_ν (even function).
///
/// Orders: ν, ν+1, …, ν+n−1. Negative orders map to |order|.
/// - All negative: single call zbesk(|ν|-(n-1), n), reverse.
/// - All positive: single call zbesk(ν, n).
/// - Integer crossing: single call covering max(|ν|, ν+n-1), mirror.
/// - Non-integer crossing: two calls (neg lattice + pos lattice).
#[cfg(feature = "alloc")]
#[allow(clippy::needless_range_loop, clippy::manual_memcpy)]
fn besselk_seq_neg<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    let zero = T::zero();
    let czero = Complex::new(zero, zero);
    let abs_nu = nu.abs();
    let nc = neg_count(nu, n);

    // All negative: absolute values are |ν|, |ν|-1, …, |ν|-(nc-1) (descending)
    // = ascending from |ν|-(nc-1) to |ν|
    if nc == n {
        // All orders negative, e.g., ν=-3.7, n=3 → orders -3.7,-2.7,-1.7
        // |orders| = 3.7, 2.7, 1.7 → start at 1.7, compute 3 ascending → reverse
        let start = abs_nu - T::from_f64((n - 1) as f64);
        let mut buf = alloc_crate::vec![czero; n];
        let (_, status) = besk::zbesk(z, start, scaling, &mut buf)?;
        buf.reverse();
        let uc = recount_underflow(&buf);
        return Ok(BesselResult {
            values: buf,
            underflow_count: uc,
            status,
        });
    }

    // Check if integer crossing (fractional part is 0 for integer nu)
    let is_int = as_integer(abs_nu).is_some();

    if is_int {
        // Integer crossing: e.g., ν=-2, n=5 → orders -2,-1,0,1,2
        // All have same fractional part (0). One call covers all |orders|.
        // Max absolute order = max(|ν|, ν+n-1)
        let last_order = nu + T::from_f64((n - 1) as f64);
        let max_abs = if abs_nu > last_order {
            abs_nu
        } else {
            last_order
        };
        let buf_n = max_abs.to_usize().unwrap_or(0) + 1; // orders 0..=max_abs
        let mut buf = alloc_crate::vec![czero; buf_n];
        let (_, status) = besk::zbesk(z, zero, scaling, &mut buf)?;

        let mut values = alloc_crate::vec![czero; n];
        for j in 0..n {
            let order = nu + T::from_f64(j as f64);
            let idx = order.abs().to_usize().unwrap_or(0);
            values[j] = buf[idx];
        }
        let uc = recount_underflow(&values);
        return Ok(BesselResult {
            values,
            underflow_count: uc,
            status,
        });
    }

    // Non-integer crossing: two separate lattices
    // Negative part: orders ν, ν+1, …, ν+(nc-1) → |orders| descending
    // frac_neg = fractional part of |ν|
    let frac_neg = abs_nu - abs_nu.floor();
    let neg_start = frac_neg; // smallest |order| in negative lattice
    let mut neg_buf = alloc_crate::vec![czero; nc];
    let (_, status1) = besk::zbesk(z, neg_start, scaling, &mut neg_buf)?;

    // Positive part: orders ν+nc, ν+nc+1, …, ν+n-1
    let pos_start = nu + T::from_f64(nc as f64);
    let pos_n = n - nc;
    let mut pos_buf = alloc_crate::vec![czero; pos_n];
    let (_, status2) = besk::zbesk(z, pos_start, scaling, &mut pos_buf)?;

    let mut values = alloc_crate::vec![czero; n];
    // Negative part: buf[0]=K(frac_neg), buf[1]=K(frac_neg+1), …
    // Need to map: order ν+j → |ν+j| → index in neg_buf
    // |ν+j| = abs_nu - j, and neg_buf[k] = K(frac_neg + k)
    // So k = |ν+j| - frac_neg = (abs_nu - j) - frac_neg = floor(abs_nu) - j
    for j in 0..nc {
        let k = (abs_nu - T::from_f64(j as f64))
            .floor()
            .to_usize()
            .unwrap_or(0);
        values[j] = neg_buf[k];
    }
    // Positive part: direct copy
    for j in 0..pos_n {
        values[nc + j] = pos_buf[j];
    }

    let status = worse_status(status1, status2);
    let uc = recount_underflow(&values);
    Ok(BesselResult {
        values,
        underflow_count: uc,
        status,
    })
}

/// H_seq with negative order support: H^(m)_{-ν} = exp(±νπi)·H^(m)_ν.
///
/// Same lattice-splitting as K, but negative-order elements get rotation.
#[cfg(feature = "alloc")]
#[allow(clippy::needless_range_loop, clippy::manual_memcpy)]
fn hankel_seq_neg<T: BesselFloat>(
    kind: HankelKind,
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    let zero = T::zero();
    let czero = Complex::new(zero, zero);
    let abs_nu = nu.abs();
    let nc = neg_count(nu, n);
    let is_int = as_integer(abs_nu).is_some();

    if nc == n {
        // All negative: compute at positive |orders| then rotate each
        let start = abs_nu - T::from_f64((n - 1) as f64);
        let mut buf = alloc_crate::vec![czero; n];
        let (_, status) = besh::zbesh(z, start, kind, scaling, &mut buf)?;
        buf.reverse();
        // Apply rotation: element j has |order| = abs_nu - j
        for j in 0..n {
            let order = abs_nu - T::from_f64(j as f64);
            buf[j] = reflect_h_element(order, kind, buf[j]);
        }
        let uc = recount_underflow(&buf);
        return Ok(BesselResult {
            values: buf,
            underflow_count: uc,
            status,
        });
    }

    if is_int {
        // Integer crossing: single call, rotate negative elements
        let last_order = nu + T::from_f64((n - 1) as f64);
        let max_abs = if abs_nu > last_order {
            abs_nu
        } else {
            last_order
        };
        let buf_n = max_abs.to_usize().unwrap_or(0) + 1;
        let mut buf = alloc_crate::vec![czero; buf_n];
        let (_, status) = besh::zbesh(z, zero, kind, scaling, &mut buf)?;

        let mut values = alloc_crate::vec![czero; n];
        for j in 0..n {
            let order = nu + T::from_f64(j as f64);
            let idx = order.abs().to_usize().unwrap_or(0);
            values[j] = buf[idx];
            if order < zero {
                values[j] = reflect_h_element(order.abs(), kind, values[j]);
            }
        }
        let uc = recount_underflow(&values);
        return Ok(BesselResult {
            values,
            underflow_count: uc,
            status,
        });
    }

    // Non-integer crossing: two lattices
    let frac_neg = abs_nu - abs_nu.floor();
    let mut neg_buf = alloc_crate::vec![czero; nc];
    let (_, status1) = besh::zbesh(z, frac_neg, kind, scaling, &mut neg_buf)?;

    let pos_start = nu + T::from_f64(nc as f64);
    let pos_n = n - nc;
    let mut pos_buf = alloc_crate::vec![czero; pos_n];
    let (_, status2) = besh::zbesh(z, pos_start, kind, scaling, &mut pos_buf)?;

    let mut values = alloc_crate::vec![czero; n];
    for j in 0..nc {
        let abs_order = abs_nu - T::from_f64(j as f64);
        let k = abs_order.floor().to_usize().unwrap_or(0);
        values[j] = reflect_h_element(abs_order, kind, neg_buf[k]);
    }
    for j in 0..pos_n {
        values[nc + j] = pos_buf[j];
    }

    let status = worse_status(status1, status2);
    let uc = recount_underflow(&values);
    Ok(BesselResult {
        values,
        underflow_count: uc,
        status,
    })
}

/// J_seq with negative order support.
/// J_{-ν} = cos(νπ)·J_ν − sin(νπ)·Y_ν (non-integer)
/// J_{-n} = (-1)^n · J_n (integer)
#[cfg(feature = "alloc")]
#[allow(clippy::needless_range_loop, clippy::manual_memcpy)]
fn besselj_seq_neg<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    let zero = T::zero();
    let czero = Complex::new(zero, zero);
    let abs_nu = nu.abs();
    let nc = neg_count(nu, n);
    let is_int = as_integer(abs_nu).is_some();

    if is_int {
        // Integer: J_{-n}(z) = (-1)^n J_n(z). Same lattice.
        let last_order = nu + T::from_f64((n - 1) as f64);
        let max_abs = if abs_nu > last_order {
            abs_nu
        } else {
            last_order
        };
        let buf_n = max_abs.to_usize().unwrap_or(0) + 1;
        let mut buf = alloc_crate::vec![czero; buf_n];
        let (_, status) = besj::zbesj(z, zero, scaling, &mut buf)?;

        let mut values = alloc_crate::vec![czero; n];
        for j in 0..n {
            let order = nu + T::from_f64(j as f64);
            let idx = order.abs().to_usize().unwrap_or(0);
            values[j] = buf[idx];
            if order < zero {
                let int_order = order.abs().to_i64().unwrap_or(0);
                values[j] = values[j] * integer_sign::<T>(int_order);
            }
        }
        let uc = recount_underflow(&values);
        return Ok(BesselResult {
            values,
            underflow_count: uc,
            status,
        });
    }

    if nc == n {
        // All negative, non-integer: need J and Y at positive |orders|
        let start = abs_nu - T::from_f64((n - 1) as f64);
        let mut j_buf = alloc_crate::vec![czero; n];
        let mut y_buf = alloc_crate::vec![czero; n];
        let (_, s1) = besj::zbesj(z, start, scaling, &mut j_buf)?;
        let (_, s2) = besy::zbesy(z, start, scaling, &mut y_buf)?;

        let mut values = alloc_crate::vec![czero; n];
        // buf[k] = f(start+k). For output j, |order| = abs_nu - j, k = |order| - start = n-1-j
        for j in 0..n {
            let k = n - 1 - j;
            let order = abs_nu - T::from_f64(j as f64);
            values[j] = reflect_j_element(order, j_buf[k], y_buf[k]);
        }
        let uc = recount_underflow(&values);
        return Ok(BesselResult {
            values,
            underflow_count: uc,
            status: worse_status(s1, s2),
        });
    }

    // Non-integer crossing: two lattices
    let frac_neg = abs_nu - abs_nu.floor();
    let mut neg_j = alloc_crate::vec![czero; nc];
    let mut neg_y = alloc_crate::vec![czero; nc];
    let (_, s1) = besj::zbesj(z, frac_neg, scaling, &mut neg_j)?;
    let (_, s2) = besy::zbesy(z, frac_neg, scaling, &mut neg_y)?;

    let pos_start = nu + T::from_f64(nc as f64);
    let pos_n = n - nc;
    let mut pos_buf = alloc_crate::vec![czero; pos_n];
    let (_, s3) = besj::zbesj(z, pos_start, scaling, &mut pos_buf)?;

    let mut values = alloc_crate::vec![czero; n];
    for j in 0..nc {
        let abs_order = abs_nu - T::from_f64(j as f64);
        let k = abs_order.floor().to_usize().unwrap_or(0);
        values[j] = reflect_j_element(abs_order, neg_j[k], neg_y[k]);
    }
    for j in 0..pos_n {
        values[nc + j] = pos_buf[j];
    }

    let status = worse_status(worse_status(s1, s2), s3);
    let uc = recount_underflow(&values);
    Ok(BesselResult {
        values,
        underflow_count: uc,
        status,
    })
}

/// Y_seq with negative order support.
/// Y_{-ν} = sin(νπ)·J_ν + cos(νπ)·Y_ν (non-integer)
/// Y_{-n} = (-1)^n · Y_n (integer)
#[cfg(feature = "alloc")]
#[allow(clippy::needless_range_loop, clippy::manual_memcpy)]
fn bessely_seq_neg<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    let zero = T::zero();
    let czero = Complex::new(zero, zero);
    let abs_nu = nu.abs();
    let nc = neg_count(nu, n);
    let is_int = as_integer(abs_nu).is_some();

    if is_int {
        let last_order = nu + T::from_f64((n - 1) as f64);
        let max_abs = if abs_nu > last_order {
            abs_nu
        } else {
            last_order
        };
        let buf_n = max_abs.to_usize().unwrap_or(0) + 1;
        let mut buf = alloc_crate::vec![czero; buf_n];
        let (_, status) = besy::zbesy(z, zero, scaling, &mut buf)?;

        let mut values = alloc_crate::vec![czero; n];
        for j in 0..n {
            let order = nu + T::from_f64(j as f64);
            let idx = order.abs().to_usize().unwrap_or(0);
            values[j] = buf[idx];
            if order < zero {
                let int_order = order.abs().to_i64().unwrap_or(0);
                values[j] = values[j] * integer_sign::<T>(int_order);
            }
        }
        let uc = recount_underflow(&values);
        return Ok(BesselResult {
            values,
            underflow_count: uc,
            status,
        });
    }

    if nc == n {
        let start = abs_nu - T::from_f64((n - 1) as f64);
        let mut j_buf = alloc_crate::vec![czero; n];
        let mut y_buf = alloc_crate::vec![czero; n];
        let (_, s1) = besj::zbesj(z, start, scaling, &mut j_buf)?;
        let (_, s2) = besy::zbesy(z, start, scaling, &mut y_buf)?;

        let mut values = alloc_crate::vec![czero; n];
        for j in 0..n {
            let k = n - 1 - j;
            let order = abs_nu - T::from_f64(j as f64);
            values[j] = reflect_y_element(order, j_buf[k], y_buf[k]);
        }
        let uc = recount_underflow(&values);
        return Ok(BesselResult {
            values,
            underflow_count: uc,
            status: worse_status(s1, s2),
        });
    }

    // Non-integer crossing: two lattices
    let frac_neg = abs_nu - abs_nu.floor();
    let mut neg_j = alloc_crate::vec![czero; nc];
    let mut neg_y = alloc_crate::vec![czero; nc];
    let (_, s1) = besj::zbesj(z, frac_neg, scaling, &mut neg_j)?;
    let (_, s2) = besy::zbesy(z, frac_neg, scaling, &mut neg_y)?;

    let pos_start = nu + T::from_f64(nc as f64);
    let pos_n = n - nc;
    let mut pos_buf = alloc_crate::vec![czero; pos_n];
    let (_, s3) = besy::zbesy(z, pos_start, scaling, &mut pos_buf)?;

    let mut values = alloc_crate::vec![czero; n];
    for j in 0..nc {
        let abs_order = abs_nu - T::from_f64(j as f64);
        let k = abs_order.floor().to_usize().unwrap_or(0);
        values[j] = reflect_y_element(abs_order, neg_j[k], neg_y[k]);
    }
    for j in 0..pos_n {
        values[nc + j] = pos_buf[j];
    }

    let status = worse_status(worse_status(s1, s2), s3);
    let uc = recount_underflow(&values);
    Ok(BesselResult {
        values,
        underflow_count: uc,
        status,
    })
}

/// I_seq with negative order support.
/// I_{-ν} = I_ν + (2/π)sin(νπ)K_ν (non-integer)
/// I_{-n} = I_n (integer)
#[cfg(feature = "alloc")]
#[allow(clippy::needless_range_loop, clippy::manual_memcpy)]
fn besseli_seq_neg<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    let zero = T::zero();
    let czero = Complex::new(zero, zero);
    let abs_nu = nu.abs();
    let nc = neg_count(nu, n);
    let is_int = as_integer(abs_nu).is_some();

    if is_int {
        // I_{-n} = I_n. Same lattice, no correction needed.
        let last_order = nu + T::from_f64((n - 1) as f64);
        let max_abs = if abs_nu > last_order {
            abs_nu
        } else {
            last_order
        };
        let buf_n = max_abs.to_usize().unwrap_or(0) + 1;
        let mut buf = alloc_crate::vec![czero; buf_n];
        let (_, status) = besi::zbesi(z, zero, scaling, &mut buf)?;

        let mut values = alloc_crate::vec![czero; n];
        for j in 0..n {
            let order = nu + T::from_f64(j as f64);
            let idx = order.abs().to_usize().unwrap_or(0);
            values[j] = buf[idx];
        }
        let uc = recount_underflow(&values);
        return Ok(BesselResult {
            values,
            underflow_count: uc,
            status,
        });
    }

    if nc == n {
        // All negative, non-integer: need I and K at positive |orders|
        let start = abs_nu - T::from_f64((n - 1) as f64);
        let mut i_buf = alloc_crate::vec![czero; n];
        let mut k_buf = alloc_crate::vec![czero; n];
        let (_, s1) = besi::zbesi(z, start, scaling, &mut i_buf)?;
        let (_, s2) = besk::zbesk(z, start, scaling, &mut k_buf)?;

        let mut values = alloc_crate::vec![czero; n];
        for j in 0..n {
            let k_idx = n - 1 - j;
            let order = abs_nu - T::from_f64(j as f64);
            let mut k_val = k_buf[k_idx];
            if scaling == Scaling::Exponential {
                k_val = k_to_i_scaling_correction(z, k_val);
            }
            values[j] = reflect_i_element(order, i_buf[k_idx], k_val);
        }
        let uc = recount_underflow(&values);
        return Ok(BesselResult {
            values,
            underflow_count: uc,
            status: worse_status(s1, s2),
        });
    }

    // Non-integer crossing: two lattices
    let frac_neg = abs_nu - abs_nu.floor();
    let mut neg_i = alloc_crate::vec![czero; nc];
    let mut neg_k = alloc_crate::vec![czero; nc];
    let (_, s1) = besi::zbesi(z, frac_neg, scaling, &mut neg_i)?;
    let (_, s2) = besk::zbesk(z, frac_neg, scaling, &mut neg_k)?;

    let pos_start = nu + T::from_f64(nc as f64);
    let pos_n = n - nc;
    let mut pos_buf = alloc_crate::vec![czero; pos_n];
    let (_, s3) = besi::zbesi(z, pos_start, scaling, &mut pos_buf)?;

    let mut values = alloc_crate::vec![czero; n];
    for j in 0..nc {
        let abs_order = abs_nu - T::from_f64(j as f64);
        let k_idx = abs_order.floor().to_usize().unwrap_or(0);
        let mut k_val = neg_k[k_idx];
        if scaling == Scaling::Exponential {
            k_val = k_to_i_scaling_correction(z, k_val);
        }
        values[j] = reflect_i_element(abs_order, neg_i[k_idx], k_val);
    }
    for j in 0..pos_n {
        values[nc + j] = pos_buf[j];
    }

    let status = worse_status(worse_status(s1, s2), s3);
    let uc = recount_underflow(&values);
    Ok(BesselResult {
        values,
        underflow_count: uc,
        status,
    })
}

#[cfg(feature = "alloc")]
/// Compute J_{ν+j}(z) for j = 0, 1, …, n−1 in a single call.
///
/// Returns a [`BesselResult`] containing `n` values and an [`Accuracy`]:
/// - [`Accuracy::Normal`] — full machine precision
/// - [`Accuracy::Reduced`] — more than half of significant digits may be lost (|z| or ν very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Negative orders are supported via DLMF reflection formulas:
/// - Non-integer ν: J_{−ν}(z) = cos(νπ) J_ν(z) − sin(νπ) Y_ν(z)
/// - Integer ν: J_{−n}(z) = (−1)^n J_n(z)
///
/// See [crate-level docs](crate#consecutive-orders) for more on sequence functions.
///
/// # Example
///
/// ```
/// # #[cfg(feature = "alloc")] {
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
/// Returns [`Error::InvalidInput`] if n < 1.
pub fn besselj_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    if nu < T::zero() {
        return besselj_seq_neg(nu, z, n, scaling);
    }
    seq_helper(n, |y| besj::zbesj(z, nu, scaling, y))
}

#[cfg(feature = "alloc")]
/// Compute Y_{ν+j}(z) for j = 0, 1, …, n−1 in a single call.
///
/// Returns a [`BesselResult`] containing `n` values and an [`Accuracy`]:
/// - [`Accuracy::Normal`] — full machine precision
/// - [`Accuracy::Reduced`] — more than half of significant digits may be lost (|z| or ν very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Negative orders are supported via DLMF reflection formulas:
/// - Non-integer ν: Y_{−ν}(z) = sin(νπ) J_ν(z) + cos(νπ) Y_ν(z)
/// - Integer ν: Y_{−n}(z) = (−1)^n Y_n(z)
///
/// See [crate-level docs](crate#consecutive-orders) for more on sequence functions.
///
/// # Example
///
/// ```
/// # #[cfg(feature = "alloc")] {
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
/// Returns [`Error::InvalidInput`] if n < 1.
pub fn bessely_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    if nu < T::zero() {
        return bessely_seq_neg(nu, z, n, scaling);
    }
    seq_helper(n, |y| besy::zbesy(z, nu, scaling, y))
}

#[cfg(feature = "alloc")]
/// Compute I_{ν+j}(z) for j = 0, 1, …, n−1 in a single call.
///
/// Returns a [`BesselResult`] containing `n` values and an [`Accuracy`]:
/// - [`Accuracy::Normal`] — full machine precision
/// - [`Accuracy::Reduced`] — more than half of significant digits may be lost (|z| or ν very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Negative orders are supported via DLMF reflection formulas:
/// - Non-integer ν: I_{−ν}(z) = I_ν(z) + (2/π) sin(νπ) K_ν(z)
/// - Integer ν: I_{−n}(z) = I_n(z)
///
/// See [crate-level docs](crate#consecutive-orders) for more on sequence functions.
///
/// # Example
///
/// ```
/// # #[cfg(feature = "alloc")] {
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
/// Returns [`Error::InvalidInput`] if n < 1.
pub fn besseli_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    if nu < T::zero() {
        return besseli_seq_neg(nu, z, n, scaling);
    }
    seq_helper(n, |y| besi::zbesi(z, nu, scaling, y))
}

#[cfg(feature = "alloc")]
/// Compute K_{ν+j}(z) for j = 0, 1, …, n−1 in a single call.
///
/// Returns a [`BesselResult`] containing `n` values and an [`Accuracy`]:
/// - [`Accuracy::Normal`] — full machine precision
/// - [`Accuracy::Reduced`] — more than half of significant digits may be lost (|z| or ν very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Negative orders are supported: K_{−ν}(z) = K_ν(z) (K is even in ν).
///
/// See [crate-level docs](crate#consecutive-orders) for more on sequence functions.
///
/// # Example
///
/// ```
/// # #[cfg(feature = "alloc")] {
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
/// Returns [`Error::InvalidInput`] if n < 1.
pub fn besselk_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    if nu < T::zero() {
        return besselk_seq_neg(nu, z, n, scaling);
    }
    seq_helper(n, |y| besk::zbesk(z, nu, scaling, y))
}

#[cfg(feature = "alloc")]
/// Compute H_{ν+j}^(1)(z) for j = 0, 1, …, n−1 in a single call.
///
/// Returns a [`BesselResult`] containing `n` values and an [`Accuracy`]:
/// - [`Accuracy::Normal`] — full machine precision
/// - [`Accuracy::Reduced`] — more than half of significant digits may be lost (|z| or ν very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Negative orders are supported: H^(1)_{−ν}(z) = exp(νπi) H^(1)_ν(z) (DLMF 10.4.6).
///
/// See [crate-level docs](crate#consecutive-orders) for more on sequence functions.
///
/// # Example
///
/// ```
/// # #[cfg(feature = "alloc")] {
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
/// Returns [`Error::InvalidInput`] if n < 1.
pub fn hankel1_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    if nu < T::zero() {
        return hankel_seq_neg(HankelKind::First, nu, z, n, scaling);
    }
    seq_helper(n, |y| besh::zbesh(z, nu, HankelKind::First, scaling, y))
}

#[cfg(feature = "alloc")]
/// Compute H_{ν+j}^(2)(z) for j = 0, 1, …, n−1 in a single call.
///
/// Returns a [`BesselResult`] containing `n` values and an [`Accuracy`]:
/// - [`Accuracy::Normal`] — full machine precision
/// - [`Accuracy::Reduced`] — more than half of significant digits may be lost (|z| or ν very large)
///
/// The `scaling` parameter selects [`Scaling::Unscaled`] or [`Scaling::Exponential`];
/// see [crate-level docs](crate#exponential-scaling) for details.
///
/// Negative orders are supported: H^(2)_{−ν}(z) = exp(−νπi) H^(2)_ν(z) (DLMF 10.4.6).
///
/// See [crate-level docs](crate#consecutive-orders) for more on sequence functions.
///
/// # Example
///
/// ```
/// # #[cfg(feature = "alloc")] {
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
/// Returns [`Error::InvalidInput`] if n < 1.
pub fn hankel2_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, Error> {
    if nu < T::zero() {
        return hankel_seq_neg(HankelKind::Second, nu, z, n, scaling);
    }
    seq_helper(n, |y| besh::zbesh(z, nu, HankelKind::Second, scaling, y))
}

// ── Tests for negative order _seq functions ──

#[cfg(all(test, feature = "alloc"))]
mod neg_order_seq_tests {
    use super::*;
    use num_complex::Complex64;

    const TOL: f64 = 2e-14;

    fn rel_err(a: Complex64, b: Complex64) -> f64 {
        let diff = (a - b).norm();
        let mag = a.norm().max(b.norm());
        if mag == 0.0 { diff } else { diff / mag }
    }

    /// Verify _seq(ν, z, n)[j] == single_value(ν+j, z) for all j.
    fn check_seq_vs_single<F, G>(
        seq_fn: F,
        single_fn: G,
        nu: f64,
        z: Complex64,
        n: usize,
        scaling: Scaling,
        tol: f64,
        label: &str,
    ) where
        F: FnOnce(f64, Complex64, usize, Scaling) -> Result<BesselResult<f64>, Error>,
        G: Fn(f64, Complex64) -> Result<Complex64, Error>,
    {
        let result = seq_fn(nu, z, n, scaling).unwrap();
        assert_eq!(result.values.len(), n, "{label}: wrong length");
        for j in 0..n {
            let order = nu + j as f64;
            let expected = single_fn(order, z).unwrap();
            let err = rel_err(result.values[j], expected);
            assert!(
                err < tol,
                "{label}: order={order}, seq={:?}, single={expected:?}, rel_err={err}",
                result.values[j]
            );
        }
    }

    // ── K (even function) ──

    #[test]
    fn besselk_seq_all_negative_nonint() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            besselk_seq,
            |nu, z| besselk(nu, z),
            -3.7,
            z,
            3,
            Scaling::Unscaled,
            TOL,
            "K all neg",
        );
    }

    #[test]
    fn besselk_seq_crossing_nonint() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            besselk_seq,
            |nu, z| besselk(nu, z),
            -1.5,
            z,
            5,
            Scaling::Unscaled,
            TOL,
            "K crossing nonint",
        );
    }

    #[test]
    fn besselk_seq_int_negative() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            besselk_seq,
            |nu, z| besselk(nu, z),
            -3.0,
            z,
            4,
            Scaling::Unscaled,
            TOL,
            "K int neg",
        );
    }

    #[test]
    fn besselk_seq_int_crossing() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            besselk_seq,
            |nu, z| besselk(nu, z),
            -2.0,
            z,
            5,
            Scaling::Unscaled,
            TOL,
            "K int cross",
        );
    }

    #[test]
    fn besselk_seq_neg_zero() {
        let z = Complex64::new(1.0, 0.0);
        check_seq_vs_single(
            besselk_seq,
            |nu, z| besselk(nu, z),
            -0.0,
            z,
            1,
            Scaling::Unscaled,
            TOL,
            "K neg zero",
        );
    }

    // ── H^(1) ──

    #[test]
    fn hankel1_seq_all_negative_nonint() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            hankel1_seq,
            |nu, z| hankel1(nu, z),
            -3.7,
            z,
            3,
            Scaling::Unscaled,
            TOL,
            "H1 all neg",
        );
    }

    #[test]
    fn hankel1_seq_crossing_nonint() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            hankel1_seq,
            |nu, z| hankel1(nu, z),
            -1.5,
            z,
            5,
            Scaling::Unscaled,
            TOL,
            "H1 crossing",
        );
    }

    #[test]
    fn hankel1_seq_int_crossing() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            hankel1_seq,
            |nu, z| hankel1(nu, z),
            -2.0,
            z,
            5,
            Scaling::Unscaled,
            TOL,
            "H1 int cross",
        );
    }

    // ── H^(2) ──

    #[test]
    fn hankel2_seq_all_negative_nonint() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            hankel2_seq,
            |nu, z| hankel2(nu, z),
            -3.7,
            z,
            3,
            Scaling::Unscaled,
            TOL,
            "H2 all neg",
        );
    }

    #[test]
    fn hankel2_seq_crossing_nonint() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            hankel2_seq,
            |nu, z| hankel2(nu, z),
            -1.5,
            z,
            5,
            Scaling::Unscaled,
            TOL,
            "H2 crossing",
        );
    }

    #[test]
    fn hankel2_seq_int_crossing() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            hankel2_seq,
            |nu, z| hankel2(nu, z),
            -2.0,
            z,
            5,
            Scaling::Unscaled,
            TOL,
            "H2 int cross",
        );
    }

    // ── J ──

    #[test]
    fn besselj_seq_all_negative_nonint() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            besselj_seq,
            |nu, z| besselj(nu, z),
            -3.7,
            z,
            3,
            Scaling::Unscaled,
            TOL,
            "J all neg",
        );
    }

    #[test]
    fn besselj_seq_crossing_nonint() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            besselj_seq,
            |nu, z| besselj(nu, z),
            -1.5,
            z,
            5,
            Scaling::Unscaled,
            TOL,
            "J crossing",
        );
    }

    #[test]
    fn besselj_seq_int_negative() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            besselj_seq,
            |nu, z| besselj(nu, z),
            -3.0,
            z,
            4,
            Scaling::Unscaled,
            TOL,
            "J int neg",
        );
    }

    #[test]
    fn besselj_seq_int_crossing() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            besselj_seq,
            |nu, z| besselj(nu, z),
            -2.0,
            z,
            5,
            Scaling::Unscaled,
            TOL,
            "J int cross",
        );
    }

    // ── Y ──

    #[test]
    fn bessely_seq_all_negative_nonint() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            bessely_seq,
            |nu, z| bessely(nu, z),
            -3.7,
            z,
            3,
            Scaling::Unscaled,
            TOL,
            "Y all neg",
        );
    }

    #[test]
    fn bessely_seq_crossing_nonint() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            bessely_seq,
            |nu, z| bessely(nu, z),
            -1.5,
            z,
            5,
            Scaling::Unscaled,
            TOL,
            "Y crossing",
        );
    }

    #[test]
    fn bessely_seq_int_negative() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            bessely_seq,
            |nu, z| bessely(nu, z),
            -3.0,
            z,
            4,
            Scaling::Unscaled,
            TOL,
            "Y int neg",
        );
    }

    #[test]
    fn bessely_seq_int_crossing() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            bessely_seq,
            |nu, z| bessely(nu, z),
            -2.0,
            z,
            5,
            Scaling::Unscaled,
            TOL,
            "Y int cross",
        );
    }

    // ── I ──

    #[test]
    fn besseli_seq_all_negative_nonint() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            besseli_seq,
            |nu, z| besseli(nu, z),
            -3.7,
            z,
            3,
            Scaling::Unscaled,
            TOL,
            "I all neg",
        );
    }

    #[test]
    fn besseli_seq_crossing_nonint() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            besseli_seq,
            |nu, z| besseli(nu, z),
            -1.5,
            z,
            5,
            Scaling::Unscaled,
            TOL,
            "I crossing",
        );
    }

    #[test]
    fn besseli_seq_int_negative() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            besseli_seq,
            |nu, z| besseli(nu, z),
            -3.0,
            z,
            4,
            Scaling::Unscaled,
            TOL,
            "I int neg",
        );
    }

    #[test]
    fn besseli_seq_int_crossing() {
        let z = Complex64::new(1.5, 0.5);
        check_seq_vs_single(
            besseli_seq,
            |nu, z| besseli(nu, z),
            -2.0,
            z,
            5,
            Scaling::Unscaled,
            TOL,
            "I int cross",
        );
    }

    // ── Scaled variants ──

    #[test]
    fn besseli_seq_scaled_negative_nonint() {
        let z = Complex64::new(2.0, 1.0);
        check_seq_vs_single(
            besseli_seq,
            |nu, z| besseli_scaled(nu, z),
            -2.3,
            z,
            5,
            Scaling::Exponential,
            TOL,
            "I scaled neg",
        );
    }

    #[test]
    fn besselj_seq_scaled_crossing() {
        let z = Complex64::new(2.0, 1.0);
        check_seq_vs_single(
            besselj_seq,
            |nu, z| besselj_scaled(nu, z),
            -1.5,
            z,
            4,
            Scaling::Exponential,
            TOL,
            "J scaled cross",
        );
    }

    #[test]
    fn besselk_seq_scaled_negative() {
        let z = Complex64::new(2.0, 1.0);
        check_seq_vs_single(
            besselk_seq,
            |nu, z| besselk_scaled(nu, z),
            -2.5,
            z,
            4,
            Scaling::Exponential,
            TOL,
            "K scaled neg",
        );
    }

    #[test]
    fn hankel1_seq_scaled_negative() {
        let z = Complex64::new(2.0, 1.0);
        check_seq_vs_single(
            hankel1_seq,
            |nu, z| hankel1_scaled(nu, z),
            -2.5,
            z,
            4,
            Scaling::Exponential,
            TOL,
            "H1 scaled neg",
        );
    }

    #[test]
    fn bessely_seq_scaled_crossing() {
        let z = Complex64::new(2.0, 1.0);
        check_seq_vs_single(
            bessely_seq,
            |nu, z| bessely_scaled(nu, z),
            -1.5,
            z,
            4,
            Scaling::Exponential,
            TOL,
            "Y scaled cross",
        );
    }

    // ── Edge: -0.0 ──

    #[test]
    fn besselj_seq_neg_zero() {
        let z = Complex64::new(1.0, 0.0);
        check_seq_vs_single(
            besselj_seq,
            |nu, z| besselj(nu, z),
            -0.0,
            z,
            3,
            Scaling::Unscaled,
            TOL,
            "J neg zero",
        );
    }
}
