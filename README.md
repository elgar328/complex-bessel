# complex-bessel

Pure Rust implementation of complex Bessel functions based on **Amos Algorithm 644** (ACM TOMS 644).

Provides Bessel functions J, Y, I, K, Hankel H<sup>(1)</sup>/H<sup>(2)</sup>, and Airy functions Ai/Bi for complex arguments and real orders.

## Installation

```toml
[dependencies]
complex-bessel = "0.1.0-alpha.1"
```

For `no_std` environments (requires `alloc`):

```toml
[dependencies]
complex-bessel = { version = "0.1", default-features = false }
```

## Usage

```rust
use complex_bessel::*;
use num_complex::Complex;

let z = Complex::new(1.0, 2.0);

let j = besselj(0.5, z).unwrap();
let k = besselk(1.0, z).unwrap();
let h = hankel(HankelKind::First, 0.0, z).unwrap();
let j_neg = besselj(-1.5, z).unwrap();  // negative orders work directly

// Scaled versions prevent overflow/underflow
let k_scaled = besselk_scaled(1.0, z).unwrap();

// Airy functions
let ai = airy(z, AiryDerivative::Value).unwrap();
```

### Sequence computation

`_seq` variants compute values at consecutive orders ν, ν+1, …, ν+n−1 in a
single call. Internal recurrence is shared, so this is more efficient than
calling the single-value function n times.

```rust
// J_0(z), J_1(z), ..., J_9(z)
let seq = besselj_seq(0.0, z, 10, Scaling::Unscaled).unwrap();

// Partial-wave expansion: Σ (2n+1) · J_n(z)
let sum: Complex<f64> = seq.values.iter().enumerate()
    .map(|(n, &jn)| jn * (2 * n + 1) as f64)
    .sum();
```

Sequence results include a `status` field (`BesselStatus::Normal` or `ReducedPrecision`) that warns when precision loss may have occurred. Single-value functions silently return the best available result.

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

where ζ = (2/3) z√z.

Each function has a `_seq` variant for consecutive-order computation (see [above](#sequence-computation)). Sequence variants require ν ≥ 0.

All single-value functions accept negative orders — DLMF reflection formulas are applied automatically.

## Accuracy

Results match the original Fortran TOMS 644 to ~14 significant digits (f64).
The `f32` generic implementation provides ~6-7 digit accuracy.
Comprehensive accuracy analysis with Fortran reference comparisons is maintained in a [separate repository](https://github.com/elgar328/complex-bessel-accuracy).

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
- MIT License ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.
