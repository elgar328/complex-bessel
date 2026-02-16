# complex-bessel

Pure Rust implementation of complex Bessel functions based on
**Amos Algorithm 644** (ACM TOMS 644).

Provides Bessel functions J, Y, I, K, Hankel H^(1)/H^(2), and Airy functions
Ai/Bi for complex arguments and real orders, with ~14-digit accuracy matching
the original Fortran implementation.

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
| `besselj` / `besselj_scaled` | J_v(z), Bessel first kind | exp(-\|Im(z)\|) · J_v(z) |
| `bessely` / `bessely_scaled` | Y_v(z), Bessel second kind | exp(-\|Im(z)\|) · Y_v(z) |
| `besseli` / `besseli_scaled` | I_v(z), modified first kind | exp(-\|Re(z)\|) · I_v(z) |
| `besselk` / `besselk_scaled` | K_v(z), modified second kind | exp(z) · K_v(z) |
| `hankel` / `hankel_scaled` | H_v^(m)(z), Hankel | exp(∓iz) · H_v^(m)(z) |
| `airy` / `airy_scaled` | Ai(z), Ai'(z) | exp(zta) · Ai(z) |
| `biry` / `biry_scaled` | Bi(z), Bi'(z) | exp(-\|Re(zta)\|) · Bi(z) |

Sequence variants (`besselj_seq`, `bessely_seq`, `besseli_seq`, `besselk_seq`,
`hankel_seq`) compute values at consecutive orders v, v+1, ..., v+n-1 in a
single call. These require v >= 0.

## Negative Order

Single-value functions accept any real order (positive or negative) using
DLMF reflection formulas:

- **K**: even in v, K_{-v}(z) = K_v(z)
- **J**: J_{-v}(z) = cos(vpi) J_v(z) - sin(vpi) Y_v(z)
- **Y**: Y_{-v}(z) = sin(vpi) J_v(z) + cos(vpi) Y_v(z)
- **I**: I_{-v}(z) = I_v(z) + (2/pi) sin(vpi) K_v(z)
- **H^(1)**: H^(1)_{-v}(z) = exp(vpi*i) H^(1)_v(z)
- **H^(2)**: H^(2)_{-v}(z) = exp(-vpi*i) H^(2)_v(z)

For integer orders, exact shortcuts are used (e.g., J_{-n} = (-1)^n J_n).

## Accuracy

Results match the original Fortran TOMS 644 to ~14 significant digits (f64).
The `f32` generic implementation provides ~6-7 digit accuracy.
Comprehensive accuracy analysis with Fortran reference comparisons is available at
[complex-bessel-accuracy](https://github.com/elgar328/complex-bessel-accuracy)
(coming soon).

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
