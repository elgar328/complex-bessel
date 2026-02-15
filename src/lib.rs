//! Pure Rust implementation of complex Bessel functions based on Amos Algorithm 644 (TOMS 644).
//!
//! This crate provides complex-valued Bessel functions of the first kind (J),
//! second kind (Y), modified first kind (I), modified second kind (K),
//! Hankel functions (H), and Airy functions (Ai, Bi).
//!
//! # Status
//!
//! This crate is in early development (alpha). Function signatures are defined
//! but implementations are not yet available.

// TODO: remove this once functions are implemented
#![allow(unused)]
#![cfg_attr(not(feature = "std"), no_std)]

#[cfg(not(feature = "std"))]
extern crate alloc;

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
pub use types::{AiryDerivative, BesselError, BesselResult, HankelKind, Scaling};

use num_complex::Complex;

// ── Single-value convenience functions ──

/// Bessel function of the first kind, J_ν(z).
pub fn besselj<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    let result = besj::zbesj(z, nu, Scaling::Unscaled, 1)?;
    Ok(result.values[0])
}

/// Bessel function of the second kind, Y_ν(z).
pub fn bessely<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    let result = besy::zbesy(z, nu, Scaling::Unscaled, 1)?;
    Ok(result.values[0])
}

/// Modified Bessel function of the first kind, I_ν(z).
pub fn besseli<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    let result = besi::zbesi(z, nu, Scaling::Unscaled, 1)?;
    Ok(result.values[0])
}

/// Modified Bessel function of the second kind, K_ν(z).
pub fn besselk<T: BesselFloat>(nu: T, z: Complex<T>) -> Result<Complex<T>, BesselError> {
    let result = besk::zbesk(z, nu, Scaling::Unscaled, 1)?;
    Ok(result.values[0])
}

/// Hankel function, H_ν^(m)(z).
pub fn hankel<T: BesselFloat>(
    kind: HankelKind,
    nu: T,
    z: Complex<T>,
) -> Result<Complex<T>, BesselError> {
    let result = besh::zbesh(z, nu, kind, Scaling::Unscaled, 1)?;
    Ok(result.values[0])
}

/// Airy function Ai(z) or its derivative Ai'(z).
pub fn airy<T: BesselFloat>(
    z: Complex<T>,
    deriv: AiryDerivative,
) -> Result<Complex<T>, BesselError> {
    todo!()
}

/// Airy function Bi(z) or its derivative Bi'(z).
pub fn biry<T: BesselFloat>(
    z: Complex<T>,
    deriv: AiryDerivative,
) -> Result<Complex<T>, BesselError> {
    todo!()
}

// ── Sequence functions with scaling option ──

/// Compute J_{ν+j}(z) for j = 0, 1, ..., n-1.
pub fn besselj_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, BesselError> {
    besj::zbesj(z, nu, scaling, n)
}

/// Compute Y_{ν+j}(z) for j = 0, 1, ..., n-1.
pub fn bessely_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, BesselError> {
    besy::zbesy(z, nu, scaling, n)
}

/// Compute I_{ν+j}(z) for j = 0, 1, ..., n-1.
pub fn besseli_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, BesselError> {
    besi::zbesi(z, nu, scaling, n)
}

/// Compute K_{ν+j}(z) for j = 0, 1, ..., n-1.
pub fn besselk_seq<T: BesselFloat>(
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, BesselError> {
    besk::zbesk(z, nu, scaling, n)
}

/// Compute H_{ν+j}^(m)(z) for j = 0, 1, ..., n-1.
pub fn hankel_seq<T: BesselFloat>(
    kind: HankelKind,
    nu: T,
    z: Complex<T>,
    n: usize,
    scaling: Scaling,
) -> Result<BesselResult<T>, BesselError> {
    besh::zbesh(z, nu, kind, scaling, n)
}
