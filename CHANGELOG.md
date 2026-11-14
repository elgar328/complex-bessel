# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `besselj_seq_into`, `bessely_seq_into`, `besseli_seq_into`, `besselk_seq_into`,
  `hankel1_seq_into`, and `hankel2_seq_into` compute the same sequences as the
  `_seq` functions into a caller-provided slice and return the new `SeqInfo`
  (underflow count and accuracy). They never allocate, for any n and ν. Where Y
  and the negative-order reflection formulas need workspace, they use a 1 KB
  stack buffer (f64) in chunks, so long sequences are somewhat slower than with
  `_seq` and can differ in the last digits, or close to the overflow limit in
  whether an error is returned.

### Fixed

- `besselj_seq`, `bessely_seq`, `besseli_seq`, `besselk_seq`, `hankel1_seq`,
  and `hankel2_seq` with ν < 0 reported in `underflow_count` only the zeros at
  the start of the sequence, so underflow at the high orders (e.g. J from
  ν = −0.5 at z = 6.5 with n = 230, whose last 11 values are zero) was
  reported as 0. They now count all underflowed values.
- `airy`, `airyprime`, `biry`, and `biryprime` (and their `_scaled` and `_raw`
  variants) returned wrong values on the negative real axis when the imaginary
  part of z was −0.0, e.g. −0.114 + 0.206i for Ai(−2 − 0i) instead of 0.227.
  As in Amos, −0 is now treated like +0.
- `bessely`, `bessely_scaled`, and `bessely_seq` with a negative half-integer
  order returned `Error::Overflow` when Y of the positive order overflows,
  although its term vanishes there (cos(νπ) = 0), e.g. `bessely(-200.5, 1.0)`
  instead of 0. Such orders now use J alone, so at z = 0 they return 0
  instead of `Error::InvalidInput` (for `bessely_seq`, when all orders are
  negative half-integers).
- `besselj`, `bessely`, `besseli`, `besselk`, `hankel1`, and `hankel2` (and
  their `_scaled`, `_seq`, and `_seq_into` variants) returned 0 or wrong values
  for |z| between about 7.2e8 and 1.07e9 (2.8e6 and 4.2e6 for f32) within
  about 30° of the real axis (J, Y, H) or the imaginary axis (I, K): for K, H,
  and Y when |ν| exceeds about 86 (34 for f32), and for J and I when ν exceeds
  about √(2|z|). E.g. `hankel1(86.0, 1e9)` returned 0 instead of
  −2.47e-5 + 5.21e-6i. Amos's uniform asymptotic expansions evaluate an Airy
  function there beyond its range limit and use the undefined result; the
  crate now evaluates it. The results have the accuracy Amos's algorithm has
  at smaller |z| (status `Accuracy::Reduced`; for f32, few correct digits).
- `besselk_scaled` (Re z < 0), `hankel1_scaled` (Im z < 0), `hankel2_scaled`
  (Im z > 0 or on the negative real axis), and `besseli_scaled` with negative
  non-integer ν and Re z < 0 (and their `_seq` and `_seq_into` variants) lost
  relative accuracy of up to about ε·|z| for orders |ν| below about 86 (34 for
  f32), where the analytic continuation's K term is not negligible (|Re z| up
  to about 10 for K and I, |Im z| for H). E.g. `hankel2_scaled(0.5, -1e8)` had
  a relative error of 7.9e-9, and 0.21 for f32 at −4e6; both are now about ε.
  Amos's ZS1S2 rounds the phase of the factor exp(−2z) to an ulp of 2·Im z;
  the crate applies exp(−2i·Im z) separately, deliberately differing from
  Amos. Other results computed through ZS1S2 (scaled Airy functions with
  Re z < 0, and the uniform asymptotic expansions, which also serve some
  unscaled J, Y, I, K, and H values) can change, but only within their
  existing error, which comes from elsewhere.
- `bessely`, `bessely_scaled`, and `bessely_seq` (and `bessely_seq_into`) with
  a negative non-integer order close to a half-integer lost accuracy, e.g. a
  relative error of 4.2e-8 for `bessely(-1.500000001, 1e-5)` and 4.8e-5 for
  `bessely(-20.499999999999, -3+0.5i)`. The reflection formula multiplies
  Y of the positive order by cos(νπ), which was computed with a relative
  error of about ε/|cos(νπ)| there. J and the Hankel functions use the same
  factor, but their results change only in the last digits.
- Values on the real and imaginary axes whose real and imaginary parts are
  real functions of very different size had a small part that could be
  wrong entirely, although the value was accurate relative to its modulus:
  `hankel1` and `hankel2` on the real axis, `bessely` (and `bessely_scaled`)
  on the negative real axis, `besselk` on the imaginary axis, and `hankel2`
  (`hankel2_scaled`) for z = ix and `hankel1` (`hankel1_scaled`) for z = −ix,
  x > 0, with their `_seq` and `_seq_into` variants. E.g. `hankel1(10.0, 1.0)`
  had the real part 8.4e-9 instead of J₁₀(1) = 2.6e-10, and
  `hankel1(1.0, 2e-32)` 1.9e15 instead of 1e-32 (openspecfun#45). Such a part
  is now computed from the two real functions (J and Y, or I and K, at |z|),
  deliberately differing from Amos, which returns it only to normwise
  accuracy. For `bessely` this includes the real part for x > |ν| where
  it is small: for orders close to a half-integer, e.g.
  `bessely(0.500000000001, -5.0)` had the real part 3.263e-13 instead of
  3.180e-13, and for 0 ≤ ν < 18 next to a zero of Y_ν, e.g.
  `bessely(0.5404395744314369, -14.199204647196062)` had 1.7790767e-7
  instead of 1.7790766e-7. Where a part is recomputed, these values take
  up to about 2 times as long (sequences up to about 4 times); elsewhere
  nothing changes. Still accurate only normwise are parts of the other
  functions that are a nearly vanishing cos(νπ) or sin(νπ) times one real
  function, for orders close to an integer or a half-integer, e.g.
  Re J_ν(−x) = cos(νπ)·J_ν(x).
- `besselj` and `bessely` (and their `_scaled`, `_seq`, and `_seq_into`
  variants) on the positive real axis returned a rounding residue as the
  imaginary part, e.g. −5.4e-16 for `besselj(0.3, 5.0)`; they are now
  exactly real.

### Removed

- The `alloc` feature. Without `std`, the whole API except the `Vec`-returning
  `_seq` functions is available and needs no allocator; use the `_seq_into`
  functions instead (with a `Vec` if an allocator is available). The `_seq`
  functions and `BesselResult` now require `std`.

## [0.4.1] - 2026-11-13

### Fixed

- J, I, Bi, and scaled Ai functions were up to about 25% slower in 0.4.0 than
  in 0.3.0 because of the new complex division. Their 0.3.0 speed is restored,
  with unchanged accuracy.

## [0.4.0] - 2026-11-13

### Changed

- `besselj_seq`, `bessely_seq`, and `besseli_seq` compute only the requested
  orders when all of them are negative integers, so cost and memory scale with
  `n` instead of |ν|. Very large |ν| now returns `Error::TotalPrecisionLoss`
  immediately instead of first allocating a buffer of |ν| elements.
- Complex division now uses Smith's algorithm instead of Amos's ZDIV. For f32,
  large-order results near the turning point (ν ≈ |z|) were up to several times
  less accurate than Amos's own single-precision routines; they now match them.
  f64 results shift slightly and are on average marginally closer to
  high-precision references. Relative changes reach ~3e-12 only in
  ill-conditioned regions, where Amos's own routines err by a similar amount.

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
- `bessely_scaled` and `bessely_seq` with `Scaling::Exponential` (and
  `besselj_scaled` / `besselj_seq` for negative non-integer ν, which use them)
  returned wrong values, reported as `Accuracy::Normal`, when |Im z| exceeds
  about 350 for f64 (40 for f32) and ν is larger than about 1.5|Im z|, e.g.
  `bessely_scaled(575.0, 360i)`. Amos's routine drops a term there as an
  underflow that in fact dominates; the crate keeps it, so it deliberately
  differs from the Fortran in this region.
- `besselj_seq` and `besseli_seq` could return every nonzero element too small
  by a factor of about 5e15 (1e7 for f32) for small |z| (power series range),
  when the highest requested orders underflow to zero and the highest nonzero
  one is just above about 1e-289 (1e-28 for f32). Lower orders were affected
  however large they are.
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

[Unreleased]: https://github.com/elgar328/complex-bessel/compare/v0.4.1...HEAD
[0.4.1]: https://github.com/elgar328/complex-bessel/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/elgar328/complex-bessel/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/elgar328/complex-bessel/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/elgar328/complex-bessel/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/elgar328/complex-bessel/releases/tag/v0.1.0
