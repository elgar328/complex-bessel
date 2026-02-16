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
let h = hankel(HankelKind::First, 0.0, z).unwrap();

// Scaled versions prevent overflow/underflow
let k_scaled = besselk_scaled(1.0, z).unwrap();

// Airy functions
let ai = airy(z, AiryDerivative::Value).unwrap();
```

## Functions

| Function | Description | Scaled variant returns |
|----------|-------------|------------------------|
| [`besselj`][d.besselj]<br>[`besselj_scaled`][d.besselj_scaled]<br>[`besselj_seq`][d.besselj_seq] | J<sub>ν</sub>(z), Bessel first kind | exp(−\|Im(z)\|) · J<sub>ν</sub>(z) |
| [`bessely`][d.bessely]<br>[`bessely_scaled`][d.bessely_scaled]<br>[`bessely_seq`][d.bessely_seq] | Y<sub>ν</sub>(z), Bessel second kind | exp(−\|Im(z)\|) · Y<sub>ν</sub>(z) |
| [`besseli`][d.besseli]<br>[`besseli_scaled`][d.besseli_scaled]<br>[`besseli_seq`][d.besseli_seq] | I<sub>ν</sub>(z), modified first kind | exp(−\|Re(z)\|) · I<sub>ν</sub>(z) |
| [`besselk`][d.besselk]<br>[`besselk_scaled`][d.besselk_scaled]<br>[`besselk_seq`][d.besselk_seq] | K<sub>ν</sub>(z), modified second kind | exp(z) · K<sub>ν</sub>(z) |
| [`hankel`][d.hankel]<br>[`hankel_scaled`][d.hankel_scaled]<br>[`hankel_seq`][d.hankel_seq] | H<sub>ν</sub><sup>(m)</sup>(z), Hankel | exp(∓iz) · H<sub>ν</sub><sup>(m)</sup>(z) |
| [`airy`][d.airy]<br>[`airy_scaled`][d.airy_scaled] | Ai(z), Ai′(z) | exp(ζ) · Ai(z) |
| [`biry`][d.biry]<br>[`biry_scaled`][d.biry_scaled] | Bi(z), Bi′(z) | exp(−\|Re(ζ)\|) · Bi(z) |

[d.besselj]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.besselj.html
[d.besselj_scaled]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.besselj_scaled.html
[d.besselj_seq]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.besselj_seq.html
[d.bessely]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.bessely.html
[d.bessely_scaled]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.bessely_scaled.html
[d.bessely_seq]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.bessely_seq.html
[d.besseli]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.besseli.html
[d.besseli_scaled]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.besseli_scaled.html
[d.besseli_seq]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.besseli_seq.html
[d.besselk]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.besselk.html
[d.besselk_scaled]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.besselk_scaled.html
[d.besselk_seq]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.besselk_seq.html
[d.hankel]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.hankel.html
[d.hankel_scaled]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.hankel_scaled.html
[d.hankel_seq]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.hankel_seq.html
[d.airy]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.airy.html
[d.airy_scaled]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.airy_scaled.html
[d.biry]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.biry.html
[d.biry_scaled]: https://docs.rs/complex-bessel/latest/complex_bessel/fn.biry_scaled.html

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
