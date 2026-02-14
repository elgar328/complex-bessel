# complex-bessel

Pure Rust implementation of complex Bessel functions based on Amos Algorithm 644 (TOMS 644).

> **Alpha**: This crate is in early development. Function signatures are defined but implementations are not yet available.

## Planned Functions

| Function | Description |
|----------|-------------|
| `besselj` | Bessel function of the first kind, J_ν(z) |
| `bessely` | Bessel function of the second kind, Y_ν(z) |
| `besseli` | Modified Bessel function of the first kind, I_ν(z) |
| `besselk` | Modified Bessel function of the second kind, K_ν(z) |
| `hankel` | Hankel functions, H_ν^(1)(z) and H_ν^(2)(z) |
| `airy` | Airy function Ai(z) and derivative Ai'(z) |
| `biry` | Airy function Bi(z) and derivative Bi'(z) |

## Features

- Generic over `f32` and `f64`
- Sequence computation for consecutive orders ν, ν+1, ..., ν+n-1
- Exponential scaling option to prevent overflow/underflow
- `no_std` compatible (with `alloc`)

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
- MIT License ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.
