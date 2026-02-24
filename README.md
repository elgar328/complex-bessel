# complex-bessel

Pure Rust implementation of complex Bessel functions based on **Amos Algorithm 644** (ACM TOMS 644).

[![crates.io](https://img.shields.io/crates/v/complex-bessel.svg)](https://crates.io/crates/complex-bessel)
[![docs.rs](https://docs.rs/complex-bessel/badge.svg)](https://docs.rs/complex-bessel)
[![license](https://img.shields.io/crates/l/complex-bessel.svg)](https://github.com/elgar328/complex-bessel)
[![CI](https://github.com/elgar328/complex-bessel/actions/workflows/ci.yml/badge.svg)](https://github.com/elgar328/complex-bessel/actions/workflows/ci.yml)

## Features

- **f32 & f64** — all functions accept `Complex<f64>` or `Complex<f32>`
- **Complete function set** — J, Y, I, K, H<sup>(1)</sup>, H<sup>(2)</sup>, Ai, Bi
- **Consecutive orders** — `_seq` variants return ν, ν+1, …, ν+n−1 in one call
- **Exponential scaling** — `_scaled` variants prevent overflow/underflow
- **Negative orders** — supports ν < 0 via DLMF reflection formulas (not in Amos)
- **`no_std` support** — 3-tier: bare `no_std`, `alloc`, `std` (default)

## Quick start

```toml
[dependencies]
complex-bessel = "0.1"
```

```rust
use complex_bessel::*;
use num_complex::Complex;

let z = Complex::new(1.0, 2.0);

let j = besselj(0.5, z).unwrap();
let k = besselk(1.0, z).unwrap();
let h = hankel1(0.0, z).unwrap();

// Scaled versions prevent overflow/underflow
let k_scaled = besselk_scaled(1.0, z).unwrap();

// Airy functions
let ai = airy(z).unwrap();
let ai_prime = airyprime(z).unwrap();
```

## Functions

| Function | Description | Scaled variant returns |
|----------|-------------|------------------------|
| `besselj(ν, z)`<br>`besselj_scaled(ν, z)`<br>`besselj_seq(ν, z, n, scaling)` | J<sub>ν</sub>(z), Bessel first kind | exp(−\|Im(z)\|) · J<sub>ν</sub>(z) |
| `bessely(ν, z)`<br>`bessely_scaled(ν, z)`<br>`bessely_seq(ν, z, n, scaling)` | Y<sub>ν</sub>(z), Bessel second kind | exp(−\|Im(z)\|) · Y<sub>ν</sub>(z) |
| `besseli(ν, z)`<br>`besseli_scaled(ν, z)`<br>`besseli_seq(ν, z, n, scaling)` | I<sub>ν</sub>(z), modified first kind | exp(−\|Re(z)\|) · I<sub>ν</sub>(z) |
| `besselk(ν, z)`<br>`besselk_scaled(ν, z)`<br>`besselk_seq(ν, z, n, scaling)` | K<sub>ν</sub>(z), modified second kind | exp(z) · K<sub>ν</sub>(z) |
| `hankel1(ν, z)`<br>`hankel1_scaled(ν, z)`<br>`hankel1_seq(ν, z, n, scaling)` | H<sub>ν</sub><sup>(1)</sup>(z), Hankel first kind | exp(−iz) · H<sub>ν</sub><sup>(1)</sup>(z) |
| `hankel2(ν, z)`<br>`hankel2_scaled(ν, z)`<br>`hankel2_seq(ν, z, n, scaling)` | H<sub>ν</sub><sup>(2)</sup>(z), Hankel second kind | exp(iz) · H<sub>ν</sub><sup>(2)</sup>(z) |
| `airy(z)`<br>`airy_scaled(z)`<br>`airy_raw(z, scaling)` | Ai(z), Airy first kind | exp(ζ) · Ai(z) |
| `airyprime(z)`<br>`airyprime_scaled(z)`<br>`airyprime_raw(z, scaling)` | Ai′(z), derivative of Airy first kind | exp(ζ) · Ai′(z) |
| `biry(z)`<br>`biry_scaled(z)`<br>`biry_raw(z, scaling)` | Bi(z), Airy second kind | exp(−\|Re(ζ)\|) · Bi(z) |
| `biryprime(z)`<br>`biryprime_scaled(z)`<br>`biryprime_raw(z, scaling)` | Bi′(z), derivative of Airy second kind | exp(−\|Re(ζ)\|) · Bi′(z) |

where ζ = (2/3) z√z.

## Function variants

The `_seq` variants (`besselj_seq`, `besselk_seq`, …) compute values at consecutive orders ν, ν+1, …, ν+n−1 in a single call and return a `BesselResult` that includes an `Accuracy` field.
The `_raw` Airy variants (`airy_raw`, `biry_raw`, …) similarly return an `AiryResult` with `Accuracy`.
All other functions return only the computed value without `Accuracy`.

**`Accuracy`:**

| Status | Meaning |
|--------|---------|
| `Normal` | Full machine precision |
| `Reduced` | More than half of significant digits may be lost |

`Reduced` is extremely rare in practice. SciPy's Bessel wrappers also silently discard the equivalent Amos IERR=3 flag by default.

To check accuracy status, use a `_seq` function (Bessel) or a `_raw` function (Airy):

```rust
use complex_bessel::*;
use num_complex::Complex;

let z = Complex::new(1.0, 2.0);
let result = besselk_seq(0.0, z, 1, Scaling::Unscaled).unwrap();
assert!(matches!(result.status, Accuracy::Normal));

let result = airy_raw(z, Scaling::Unscaled).unwrap();
assert!(matches!(result.status, Accuracy::Normal));
```

## Error handling

All functions return `Result<_, Error>`. The four error variants are:

| Variant | Cause |
|---------|-------|
| `InvalidInput` | z = 0 for K/Y/H, n < 1 |
| `Overflow` | \|z\| or ν too large (or too small) for finite result |
| `TotalPrecisionLoss` | Complete loss of significant digits; \|z\| or ν too large |
| `ConvergenceFailure` | Internal algorithm did not converge |

## `no_std` support

| Cargo features | Available API |
|---------------|---------------|
| `default-features = false` | 24 single-value functions |
| `features = ["alloc"]` | full API (+ 6 `_seq` variants) |
| `features = ["std"]` (default) | full API, with platform-native math |

The 24 single-value functions include 12 Bessel (J/Y/I/K/H<sup>(1)</sup>/H<sup>(2)</sup> × unscaled/scaled) and 12 Airy (Ai/Ai'/Bi/Bi' × unscaled/scaled/raw).

```toml
# Bare no_std — no allocator needed:
complex-bessel = { version = "0.1", default-features = false }

# no_std + alloc — full API:
complex-bessel = { version = "0.1", default-features = false, features = ["alloc"] }
```

## Accuracy & performance

Precision of complex-bessel (Rust), AMOS/TOMS 644 (Fortran), and SciPy (Python) across all functions.

[![Accuracy — relative error vs mpmath](https://raw.githubusercontent.com/elgar328/complex-bessel-test/main/images/accuracy.svg)](https://github.com/elgar328/complex-bessel-test/blob/main/images/accuracy.pdf)

Since this library is a Rust translation of AMOS/TOMS 644, the two should produce nearly identical results.

[![Rust vs Fortran — relative error](https://raw.githubusercontent.com/elgar328/complex-bessel-test/main/images/fidelity.svg)](https://github.com/elgar328/complex-bessel-test/blob/main/images/fidelity.pdf)

Median evaluation time per function call.

[![Performance — median time per call](https://raw.githubusercontent.com/elgar328/complex-bessel-test/main/images/eval_time.svg)](https://github.com/elgar328/complex-bessel-test/blob/main/images/eval_time.pdf)

For methodology, datasets, and full benchmark results, see the [test repository](https://github.com/elgar328/complex-bessel-test).

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
- MIT License ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.
