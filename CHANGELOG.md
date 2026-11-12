# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- `besselj_seq`, `bessely_seq`, and `besseli_seq` compute only the requested
  orders when all of them are negative integers, so cost and memory scale with
  `n` instead of |ν|. Very large |ν| now returns `Error::TotalPrecisionLoss`
  immediately instead of first allocating a buffer of |ν| elements.
- Complex division now uses Smith's algorithm instead of Amos's ZDIV. For f32,
  large-order results near the turning point (ν ≈ |z|) were up to several times
  less accurate than Amos's own single-precision routines; they now match them.
  f64 results shift slightly (relative changes up to ~5e-13 in ill-conditioned
  regions) and are on average marginally closer to high-precision references.

### Fixed

- `besselj_seq`, `bessely_seq`, `besseli_seq`, `hankel1_seq`, and `hankel2_seq`
  returned wrong values for all but the last element when `n ≥ 3` and the
  uniform asymptotic expansion was used (roughly when ν + n − 1 exceeds 86 for
  f64, 34 for f32, and |z| is not small), while reporting `Accuracy::Normal`.
  Single-value functions were not affected.
- `besselj_seq`, `bessely_seq`, `besseli_seq`, `hankel1_seq`, and `hankel2_seq`
  could return elements off by a large factor (about 1/ε) for `n ≥ 3` when |z|
  exceeds about 86 (f64; 34 for f32) and the values are large enough to trigger
  internal rescaling (roughly above 1e289 for f64 and 1e28 for f32, so f32 was
  affected far more often). Single-value functions were not affected.
- `besselk_seq`, `hankel1_seq`, and `hankel2_seq` without scaling returned wrong
  values from the third element on (relative error around 1e-3, sometimes zeros)
  when the K argument has a real part above about 665 for f64 (64 for f32), where
  results approach the underflow limit. Single-value functions were not affected.
- `besselk_seq` returned wrong values when ν is below about 86 (f64; 34 for f32)
  but ν + n − 1 is above it, e.g. `besselk_seq(0.5, 2+i, 100)`, because it
  applied the large-order expansion from the first order.
- `besselk` and `besselk_seq` returned zeros instead of `Error::Overflow` for
  Re z < 0 when the result overflows.
- An internal threshold test for results near the underflow and overflow limits
  used a wrong phase, so some results were zero (or `Error::Overflow`) where
  Amos's algorithm computes them, or the reverse:
  - J, Y, I, K, and Hankel functions (single values and `_seq`, including
    `underflow_count`) for magnitudes of about 1e-305 to 1e-289 for f64
    (1e-35 to 1e-28 for f32).
  - `besselk` and the Hankel functions without scaling, for magnitudes of about
    1e289 to 1e304 for f64 (1e28 to 1e35 for f32) when Re z < 0 for K,
    Im z < 0 for H(1), or Im z > 0 for H(2). Here the double-precision Fortran
    evaluates the test in the wrong half plane and reports `Error::Overflow`
    for many representable results; the crate follows Amos's single-precision
    routines instead.
- `besselj_seq`, `bessely_seq`, and `besseli_seq` no longer return
  `Error::Overflow` for negative integer orders when only orders below the
  requested range would overflow.
- All `_seq` functions now return `Error::InvalidInput` for `n = 0` with negative
  ν, as documented. Previously they panicked in debug builds; in release builds,
  `besselj_seq`, `bessely_seq`, and `besseli_seq` returned an empty result for
  integer ν.

## [0.3.0] - 2026-08-06

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

[Unreleased]: https://github.com/elgar328/complex-bessel/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/elgar328/complex-bessel/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/elgar328/complex-bessel/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/elgar328/complex-bessel/releases/tag/v0.1.0
