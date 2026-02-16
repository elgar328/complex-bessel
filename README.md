# complex-bessel

Pure Rust implementation of complex Bessel functions based on **Amos Algorithm 644** (ACM TOMS 644).

[![crates.io](https://img.shields.io/crates/v/complex-bessel.svg)](https://crates.io/crates/complex-bessel)
[![docs.rs](https://docs.rs/complex-bessel/badge.svg)](https://docs.rs/complex-bessel)
[![license](https://img.shields.io/crates/l/complex-bessel.svg)](https://github.com/elgar328/complex-bessel)

## Features

- **Dual precision** — all functions accept `Complex<f64>` or `Complex<f32>`
- **Complete function set** — J, Y, I, K, H<sup>(1)</sup>, H<sup>(2)</sup>, Ai, Bi
- **Consecutive orders** — `_seq` variants return ν, ν+1, …, ν+n−1 in one call
- **Exponential scaling** — `_scaled` variants prevent overflow/underflow
- **Negative orders** — single-value functions accept ν < 0 via reflection formulas
- **`no_std`** — works with `alloc` only

## Quick start

```toml
[dependencies]
complex-bessel = "0.1.0-alpha.1"
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
| `airy(z)`<br>`airy_scaled(z)` | Ai(z), Airy first kind | exp(ζ) · Ai(z) |
| `airyprime(z)`<br>`airyprime_scaled(z)` | Ai′(z), derivative of Airy first kind | exp(ζ) · Ai′(z) |
| `biry(z)`<br>`biry_scaled(z)` | Bi(z), Airy second kind | exp(−\|Re(ζ)\|) · Bi(z) |
| `biryprime(z)`<br>`biryprime_scaled(z)` | Bi′(z), derivative of Airy second kind | exp(−\|Re(ζ)\|) · Bi′(z) |

where ζ = (2/3) z√z.

## Accuracy

Results match the original Fortran TOMS 644 to ~14 significant digits (f64).
The `f32` generic implementation provides ~6–7 digit accuracy.
Comprehensive accuracy analysis is maintained in a [separate repository](https://github.com/elgar328/complex-bessel-accuracy).

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
- MIT License ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.
