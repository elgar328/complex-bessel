# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- Reject non-finite inputs (NaN or infinite z/ν) with `Error::InvalidInput`
  instead of panicking at internal `to_i32().unwrap()` calls (all Bessel
  functions).
- Return `Error::TotalPrecisionLoss` from Y functions for finite orders
  ν ≥ 2³¹ instead of panicking (`zbesy` lacked the range check the other
  upper interfaces have).
- Report `Accuracy::Reduced` from Y functions when |z| or ν exceeds the
  half-precision threshold (~32767 for f64), matching the other functions
  and Fortran ZBESY's IERR=3; previously the status was always `Normal`.

### Changed

- Infinite z or ν now classifies as `Error::InvalidInput` instead of
  `Error::TotalPrecisionLoss` (Bessel and Airy functions).
- Airy functions now return `Error::InvalidInput` for NaN z instead of
  silently returning `Ok(NaN)`.

## [0.2.0] - 2026-02-28

### Changed

- Use unconditional `#![no_std]` with `core::error::Error` (MSRV 1.85+, edition 2024).
- Apply FMA (fused multiply-add) optimization across all arithmetic hot paths
  (`mul_add`, `mul_add_scalar`, scalar `fma` — 38 call sites in 19 files).
- Improve doc comments and README for public API, types, and module headers.

## [0.1.0] - 2026-02-22

### Added

- Implement complex Bessel functions J, Y, I, K, H(1), H(2) and Airy Ai, Ai', Bi, Bi'.
- Support both f32 and f64 precision.
- Accept complex arguments with fractional orders.
- Support negative orders (ν < 0) via DLMF reflection formulas.
- Provide exponentially scaled variants (`_scaled`) for all functions.
- Provide sequence computation (`_seq` variants) for consecutive orders.
- Support `no_std` with 3-tier feature flags (`no_std` / `alloc` / `std`).

[Unreleased]: https://github.com/elgar328/complex-bessel/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/elgar328/complex-bessel/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/elgar328/complex-bessel/releases/tag/v0.1.0
