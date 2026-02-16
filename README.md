# complex-bessel

Pure Rust implementation of complex Bessel functions based on
**Amos Algorithm 644** (ACM TOMS 644).

Provides Bessel functions J, Y, I, K, Hankel H<sup>(1)</sup>/H<sup>(2)</sup>, and Airy functions
Ai/Bi for complex arguments and real orders, with ~14-digit accuracy matching
the original Fortran implementation.

## Installation

```toml
[dependencies]
complex-bessel = "0.1.0-alpha.1"
```

## Usage

```rust
use complex_bessel::*;
use num_complex::Complex;

let z = Complex::new(1.0, 2.0);

// Single-value functions (support negative order)
let j = besselj(0.5, z).unwrap();
let k = besselk(1.0, z).unwrap();
let h = hankel(HankelKind::First, 0.0, z).unwrap();

// Negative order via DLMF reflection formulas
let j_neg = besselj(-1.5, z).unwrap();
let k_neg = besselk(-1.0, z).unwrap();  // K_{-v} = K_v

// Scaled computation to prevent overflow/underflow
let k_scaled = besselk_scaled(1.0, z).unwrap();  // exp(z) * K_1(z)

// Sequence: compute K_0(z), K_1(z), K_2(z)
let seq = besselk_seq(0.0, z, 3, Scaling::Unscaled).unwrap();
assert_eq!(seq.values.len(), 3);
assert_eq!(seq.status, BesselStatus::Normal);  // precision status

// Airy functions
let ai = airy(z, AiryDerivative::Value).unwrap();
let bi_prime = biry(z, AiryDerivative::Derivative).unwrap();
```

## Functions

| Function | Description | Scaled version returns |
|----------|-------------|------------------------|
| `besselj` / `besselj_scaled` | J<sub>ν</sub>(z), Bessel first kind | exp(−\|Im(z)\|) · J<sub>ν</sub>(z) |
| `bessely` / `bessely_scaled` | Y<sub>ν</sub>(z), Bessel second kind | exp(−\|Im(z)\|) · Y<sub>ν</sub>(z) |
| `besseli` / `besseli_scaled` | I<sub>ν</sub>(z), modified first kind | exp(−\|Re(z)\|) · I<sub>ν</sub>(z) |
| `besselk` / `besselk_scaled` | K<sub>ν</sub>(z), modified second kind | exp(z) · K<sub>ν</sub>(z) |
| `hankel` / `hankel_scaled` | H<sub>ν</sub><sup>(m)</sup>(z), Hankel | exp(∓iz) · H<sub>ν</sub><sup>(m)</sup>(z) |
| `airy` / `airy_scaled` | Ai(z), Ai′(z) | exp(ζ) · Ai(z) |
| `biry` / `biry_scaled` | Bi(z), Bi′(z) | exp(−\|Re(ζ)\|) · Bi(z) |

Sequence variants (`besselj_seq`, `bessely_seq`, `besseli_seq`, `besselk_seq`,
`hankel_seq`) compute values at consecutive orders ν, ν+1, …, ν+n−1 in a
single call. These require ν ≥ 0.

## Negative Order

Single-value functions accept any real order (positive or negative) using
DLMF reflection formulas:

- **K**: even in ν, K<sub>−ν</sub>(z) = K<sub>ν</sub>(z)
- **J**: J<sub>−ν</sub>(z) = cos(νπ) J<sub>ν</sub>(z) − sin(νπ) Y<sub>ν</sub>(z)
- **Y**: Y<sub>−ν</sub>(z) = sin(νπ) J<sub>ν</sub>(z) + cos(νπ) Y<sub>ν</sub>(z)
- **I**: I<sub>−ν</sub>(z) = I<sub>ν</sub>(z) + (2/π) sin(νπ) K<sub>ν</sub>(z)
- **H<sup>(1)</sup>**: H<sup>(1)</sup><sub>−ν</sub>(z) = exp(νπi) H<sup>(1)</sup><sub>ν</sub>(z)
- **H<sup>(2)</sup>**: H<sup>(2)</sup><sub>−ν</sub>(z) = exp(−νπi) H<sup>(2)</sup><sub>ν</sub>(z)

For integer orders, exact shortcuts are used (e.g., J<sub>−n</sub> = (−1)<sup>n</sup> J<sub>n</sub>).

## Accuracy

Results match the original Fortran TOMS 644 to ~14 significant digits (f64).
The `f32` generic implementation provides ~6-7 digit accuracy.
Comprehensive accuracy analysis with Fortran reference comparisons is maintained in a
[separate repository](https://github.com/elgar328/complex-bessel-accuracy).

## `no_std`

Disable the default `std` feature to use in `no_std` environments (requires `alloc`):

```toml
[dependencies]
complex-bessel = { version = "0.1", default-features = false }
```

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
- MIT License ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.
