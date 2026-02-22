# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2026-02-22

### Added

- Implement complex Bessel functions J, Y, I, K, H(1), H(2) and Airy Ai, Ai', Bi, Bi'.
- Support both f32 and f64 precision.
- Accept complex arguments with fractional orders.
- Support negative orders (ν < 0) via DLMF reflection formulas.
- Provide exponentially scaled variants (`_scaled`) for all functions.
- Provide sequence computation (`_seq` variants) for consecutive orders.
- Support `no_std` with 3-tier feature flags (`no_std` / `alloc` / `std`).

[Unreleased]: https://github.com/elgar328/complex-bessel/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/elgar328/complex-bessel/releases/tag/v0.1.0
