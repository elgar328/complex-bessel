# complex-bessel

Pure Rust implementation of complex Bessel functions based on **Amos Algorithm 644** (ACM TOMS 644).

[![crates.io](https://img.shields.io/crates/v/complex-bessel.svg)](https://crates.io/crates/complex-bessel)
[![docs.rs](https://docs.rs/complex-bessel/badge.svg)](https://docs.rs/complex-bessel)
[![license](https://img.shields.io/crates/l/complex-bessel.svg)](https://github.com/elgar328/complex-bessel)
[![CI](https://github.com/elgar328/complex-bessel/actions/workflows/ci.yml/badge.svg)](https://github.com/elgar328/complex-bessel/actions/workflows/ci.yml)

## Features

- **f32 & f64** — all functions accept `Complex<f64>` or `Complex<f32>`
- **Complete function set** — J, Y, I, K, H<sup>(1)</sup>, H<sup>(2)</sup>, Ai, Ai′, Bi, Bi′
- **Consecutive orders** — `_seq` and `_seq_into` variants evaluate at orders ν, ν+1, …, ν+n−1 in one call
- **Exponential scaling** — `_scaled` variants prevent overflow/underflow
- **Negative orders** — supports ν < 0 via DLMF reflection formulas (not in Amos)
- **`no_std`** — no allocator needed; only the `Vec`-returning `_seq` variants require `std`

## Accuracy & performance

Relative error against mpmath for complex-bessel (Rust), Amos TOMS 644 (Fortran), SciPy (Python), and SpecialFunctions.jl (Julia).

[![Accuracy — relative error vs mpmath](https://raw.githubusercontent.com/elgar328/complex-bessel-test/main/images/accuracy.svg)](https://github.com/elgar328/complex-bessel-test/blob/main/images/accuracy.pdf)

Evaluation time per call.

[![Performance — median time per call](https://raw.githubusercontent.com/elgar328/complex-bessel-test/main/images/eval_time.svg)](https://github.com/elgar328/complex-bessel-test/blob/main/images/eval_time.pdf)

The charts are for f64. The [test repository](https://github.com/elgar328/complex-bessel-test) has the methodology, full results, and f32 accuracy.

## Quick start

```toml
[dependencies]
complex-bessel = "0.4"
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
| `besselj(ν, z)`<br>`besselj_scaled(ν, z)`<br>`besselj_seq(ν, z, n, scaling)`<br>`besselj_seq_into(ν, z, scaling, out)` | J<sub>ν</sub>(z), Bessel first kind | exp(−\|Im(z)\|) · J<sub>ν</sub>(z) |
| `bessely(ν, z)`<br>`bessely_scaled(ν, z)`<br>`bessely_seq(ν, z, n, scaling)`<br>`bessely_seq_into(ν, z, scaling, out)` | Y<sub>ν</sub>(z), Bessel second kind | exp(−\|Im(z)\|) · Y<sub>ν</sub>(z) |
| `besseli(ν, z)`<br>`besseli_scaled(ν, z)`<br>`besseli_seq(ν, z, n, scaling)`<br>`besseli_seq_into(ν, z, scaling, out)` | I<sub>ν</sub>(z), modified first kind | exp(−\|Re(z)\|) · I<sub>ν</sub>(z) |
| `besselk(ν, z)`<br>`besselk_scaled(ν, z)`<br>`besselk_seq(ν, z, n, scaling)`<br>`besselk_seq_into(ν, z, scaling, out)` | K<sub>ν</sub>(z), modified second kind | exp(z) · K<sub>ν</sub>(z) |
| `hankel1(ν, z)`<br>`hankel1_scaled(ν, z)`<br>`hankel1_seq(ν, z, n, scaling)`<br>`hankel1_seq_into(ν, z, scaling, out)` | H<sub>ν</sub><sup>(1)</sup>(z), Hankel first kind | exp(−iz) · H<sub>ν</sub><sup>(1)</sup>(z) |
| `hankel2(ν, z)`<br>`hankel2_scaled(ν, z)`<br>`hankel2_seq(ν, z, n, scaling)`<br>`hankel2_seq_into(ν, z, scaling, out)` | H<sub>ν</sub><sup>(2)</sup>(z), Hankel second kind | exp(iz) · H<sub>ν</sub><sup>(2)</sup>(z) |
| `airy(z)`<br>`airy_scaled(z)`<br>`airy_raw(z, scaling)` | Ai(z), Airy first kind | exp(ζ) · Ai(z) |
| `airyprime(z)`<br>`airyprime_scaled(z)`<br>`airyprime_raw(z, scaling)` | Ai′(z), derivative of Airy first kind | exp(ζ) · Ai′(z) |
| `biry(z)`<br>`biry_scaled(z)`<br>`biry_raw(z, scaling)` | Bi(z), Airy second kind | exp(−\|Re(ζ)\|) · Bi(z) |
| `biryprime(z)`<br>`biryprime_scaled(z)`<br>`biryprime_raw(z, scaling)` | Bi′(z), derivative of Airy second kind | exp(−\|Re(ζ)\|) · Bi′(z) |

where ζ = (2/3) z√z.

## Function variants

The `_seq` variants (`besselj_seq`, `besselk_seq`, …) compute values at consecutive orders ν, ν+1, …, ν+n−1 in a single call and return a `BesselResult` that includes an `Accuracy` field.
The `_seq_into` variants (`besselj_seq_into`, `besselk_seq_into`, …) write the sequence into a caller-provided slice and return a `SeqInfo` with `Accuracy`; they never allocate.
The `_raw` Airy variants (`airy_raw`, `biry_raw`, …) similarly return an `AiryResult` with `Accuracy`.
All other functions return only the computed value without `Accuracy`.

**`Accuracy`:**

| Status | Meaning |
|--------|---------|
| `Normal` | Full machine precision |
| `Reduced` | More than half of significant digits may be lost |

`Reduced` is extremely rare in practice. SciPy's Bessel wrappers also silently discard the equivalent Amos IERR=3 flag by default.

To check the accuracy status, use a `_seq`/`_seq_into` function (Bessel) or a `_raw` function (Airy):

```rust
use complex_bessel::*;
use num_complex::Complex;

let z = Complex::new(1.0, 2.0);
let mut out = [Complex::new(0.0, 0.0); 1];
let info = besselk_seq_into(0.0, z, Scaling::Unscaled, &mut out).unwrap();
assert!(matches!(info.status, Accuracy::Normal));

let result = airy_raw(z, Scaling::Unscaled).unwrap();
assert!(matches!(result.status, Accuracy::Normal));
```

## Error handling

All functions return `Result<_, Error>`. The four error variants are:

| Variant | Cause |
|---------|-------|
| `InvalidInput` | non-finite z or ν; z = 0 for K, H, Y (except at negative half-integer orders) and for J, I of negative non-integer order; n < 1 (empty `out`) |
| `Overflow` | \|z\| or \|ν\| too large (or \|z\| too small) for a finite result |
| `TotalPrecisionLoss` | \|z\| or \|ν\| too large for meaningful computation |
| `ConvergenceFailure` | internal algorithm did not converge |

## Using without `std`

Turn off the default `std` feature:

```toml
complex-bessel = { version = "0.4", default-features = false }
```

All functions except the six `_seq` variants (which return a `Vec`) remain
available, and none of them allocate. Math functions then come from `libm`
instead of the platform library.

## Improvements over Amos

- **Scaled Y at large |Im z|** — Amos's ZBESY drops a dominant term as an underflow (large ν, |Im z| > ~350); complex-bessel keeps it.
- **No spurious overflow in K and H** — double-precision ZUOIK reports `Overflow` for representable values near 1e300; complex-bessel follows the single-precision CUOIK test.
- **Complex division** — Smith's algorithm instead of ZDIV, for better f32 accuracy.
- **Airy limit in uniform expansions** — at |z| ≈ 7e8–1e9 and large ν, ZUNI2/ZUNK2 use an undefined result when the Airy argument exceeds its range limit; complex-bessel evaluates the Airy function there anyway.
- **Small components on the axes** — where a value combines two real functions of very different size (e.g. H = J + iY on the real axis), Amos's small component can lose all its digits; complex-bessel recomputes it from the real functions.
- **Phase of scaled K and H** — in the analytic continuation, Amos's ZS1S2 rounds the phase of exp(−2i Im z), a relative error up to ε·|z|; complex-bessel applies this factor separately, keeping the phase exact.

The last two show up in the componentwise relative error (of the real and imaginary parts separately):

[![Componentwise accuracy — on the axes, and of scaled functions near them](https://raw.githubusercontent.com/elgar328/complex-bessel-test/main/images/componentwise_accuracy.svg)](https://github.com/elgar328/complex-bessel-test/blob/main/images/componentwise_accuracy.pdf)

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
- MIT License ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.
