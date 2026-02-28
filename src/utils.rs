//! Overflow-safe complex arithmetic utilities.
//!
//! Translations of Amos ZABS and ZDIV from TOMS 644.

use num_complex::Complex;

use crate::machine::BesselFloat;

/// Fused multiply-add for Complex: `s * a + b`.
///
/// Uses hardware FMA (std) or plain arithmetic (no_std) via
/// [`BesselFloat::fma`], avoiding the slow software FMA in libm
/// that `num_complex::Complex::mul_add` would invoke under `no_std`.
#[inline]
pub(crate) fn mul_add<T: BesselFloat>(s: Complex<T>, a: Complex<T>, b: Complex<T>) -> Complex<T> {
    // (s.re + s.im*i)(a.re + a.im*i) + (b.re + b.im*i)
    // re = s.re*a.re - s.im*a.im + b.re
    // im = s.re*a.im + s.im*a.re + b.im
    Complex::new(
        BesselFloat::fma(s.re, a.re, b.re) - s.im * a.im,
        BesselFloat::fma(s.re, a.im, BesselFloat::fma(s.im, a.re, b.im)),
    )
}

/// Fused multiply-add for Complex × scalar: `s * a + b`.
///
/// More efficient than full complex `mul_add` — no cross-term multiplications.
/// Uses hardware FMA (std) or plain arithmetic (no_std) via [`BesselFloat::fma`].
#[inline]
pub(crate) fn mul_add_scalar<T: BesselFloat>(s: Complex<T>, a: T, b: Complex<T>) -> Complex<T> {
    Complex::new(s.re.fma(a, b.re), s.im.fma(a, b.im))
}

/// Multiply a complex number by i: (a+bi)·i = -b+ai.
#[inline]
pub(crate) fn mul_i<T: BesselFloat>(c: Complex<T>) -> Complex<T> {
    Complex::new(-c.im, c.re)
}

/// Multiply a complex number by -i: (a+bi)·(-i) = b-ai.
#[inline]
pub(crate) fn mul_neg_i<T: BesselFloat>(c: Complex<T>) -> Complex<T> {
    Complex::new(c.im, -c.re)
}

/// Overflow-safe complex absolute value.
///
/// Computes `|z| = sqrt(re² + im²)` without intermediate overflow by
/// factoring out the larger component:
///   `max * sqrt(1 + (min/max)²)`
///
/// Equivalent to Fortran ZABS in TOMS 644.
#[inline]
pub(crate) fn zabs<T: BesselFloat>(z: Complex<T>) -> T {
    let u = z.re.abs();
    let v = z.im.abs();
    let s = u + v;
    if s == T::zero() {
        return T::zero();
    }
    if u > v {
        let q = v / u;
        u * (T::one() + q * q).sqrt()
    } else {
        let q = u / v;
        v * (T::one() + q * q).sqrt()
    }
}

/// Overflow-safe complex division `a / b`.
///
/// Uses the Amos algorithm: normalize b by its magnitude first to avoid
/// overflow in intermediate products. Computes:
///   `bm = 1 / |b|`
///   `c = (a.re * b.re + a.im * b.im, a.im * b.re - a.re * b.im) * bm²`
///
/// Equivalent to Fortran ZDIV in TOMS 644.
#[inline]
pub(crate) fn zdiv<T: BesselFloat>(a: Complex<T>, b: Complex<T>) -> Complex<T> {
    let bm = T::one() / zabs(b);
    let cc = b.re * bm;
    let cd = b.im * bm;
    Complex::new((a.re * cc + a.im * cd) * bm, (a.im * cc - a.re * cd) * bm)
}

/// Overflow-safe reciprocal: `rz = 2 / z`.
///
/// Computes `2/z` without intermediate overflow by factoring out the
/// magnitude. Equivalent to the `RZ` computation in Fortran ZBKNU,
/// ZSERI, ZASYI, etc.
#[inline]
pub(crate) fn reciprocal_z<T: BesselFloat>(z: Complex<T>) -> Complex<T> {
    let one = T::one();
    let raz = one / zabs(z);
    let str = z.re * raz;
    let sti = -z.im * raz;
    Complex::new((str + str) * raz, (sti + sti) * raz)
}

/// Compute sin(π·x) with exact values at half-integers.
///
/// Reduces the argument modulo 2 first, so `sinpi(n)` is exactly 0 for
/// any integer `n`, and `sinpi(n + 0.5)` is exactly ±1. This avoids the
/// catastrophic rounding errors of `(x * PI).sin()` when x is a
/// half-integer (e.g. `sin(1.5 * PI)` = −1.837e-16 instead of 0).
///
/// Algorithm follows scipy/xsf: reduce to [0, 0.5], use symmetry.
#[inline]
pub(crate) fn sinpi<T: BesselFloat>(x: T) -> T {
    let zero = T::zero();
    let one = T::one();
    let two = T::from_f64(2.0);
    let half = T::from_f64(0.5);
    let one_half = T::from_f64(1.5);
    let pi = T::from_f64(core::f64::consts::PI);

    // sinpi is odd: sinpi(-x) = -sinpi(x)
    let (ax, sign) = if x < zero { (-x, -one) } else { (x, one) };

    // Reduce to [0, 2): r = ax mod 2
    let r = ax % two;

    // Exact special values
    if r == zero || r == one {
        return zero;
    }
    if r == half {
        return sign;
    }
    if r == one_half {
        return -sign;
    }

    // Use symmetry to reduce to [0, 0.5]
    let s = if r < half {
        (r * pi).sin()
    } else if r < one {
        ((one - r) * pi).sin()
    } else if r < one_half {
        -((r - one) * pi).sin()
    } else {
        -((two - r) * pi).sin()
    };

    sign * s
}

/// Compute cos(π·x) with exact values at integers and half-integers.
///
/// Reduces the argument modulo 2 first, so `cospi(n + 0.5)` is exactly 0
/// for any integer `n`, and `cospi(n)` is exactly ±1. This avoids the
/// catastrophic rounding errors of `(x * PI).cos()` when x is a
/// half-integer (e.g. `cos(1.5 * PI)` = −1.837e-16 instead of 0).
///
/// Algorithm follows scipy/xsf: reduce to [0, 0.5], use symmetry.
#[inline]
pub(crate) fn cospi<T: BesselFloat>(x: T) -> T {
    let zero = T::zero();
    let one = T::one();
    let two = T::from_f64(2.0);
    let half = T::from_f64(0.5);
    let one_half = T::from_f64(1.5);
    let pi = T::from_f64(core::f64::consts::PI);

    // cospi is even: cospi(-x) = cospi(x)
    let ax = x.abs();

    // Reduce to [0, 2): r = ax mod 2
    let r = ax % two;

    // Exact special values
    if r == zero {
        return one;
    }
    if r == half || r == one_half {
        return zero;
    }
    if r == one {
        return -one;
    }

    // Use symmetry to reduce to [0, 0.5]
    if r < half {
        (r * pi).cos()
    } else if r < one {
        -((one - r) * pi).cos()
    } else if r < one_half {
        -((r - one) * pi).cos()
    } else {
        ((two - r) * pi).cos()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    // ── zabs tests ──

    #[test]
    fn zabs_zero() {
        assert_eq!(zabs(Complex64::new(0.0, 0.0)), 0.0);
    }

    #[test]
    fn zabs_real_only() {
        let z = Complex64::new(3.0, 0.0);
        assert!((zabs(z) - 3.0).abs() < 1e-15);

        let z = Complex64::new(-5.0, 0.0);
        assert!((zabs(z) - 5.0).abs() < 1e-15);
    }

    #[test]
    fn zabs_imag_only() {
        let z = Complex64::new(0.0, 4.0);
        assert!((zabs(z) - 4.0).abs() < 1e-15);
    }

    #[test]
    fn zabs_3_4_triangle() {
        let z = Complex64::new(3.0, 4.0);
        assert!((zabs(z) - 5.0).abs() < 1e-15);
    }

    #[test]
    fn zabs_large_values_no_overflow() {
        // Values near sqrt(f64::MAX) would overflow with naive re²+im²
        let big = 1.0e154;
        let z = Complex64::new(big, big);
        let result = zabs(z);
        let expected = big * 2.0_f64.sqrt();
        assert!((result - expected).abs() / expected < 1e-15);
    }

    #[test]
    fn zabs_tiny_values_no_underflow() {
        let tiny = 1.0e-308;
        let z = Complex64::new(tiny, tiny);
        let result = zabs(z);
        assert!(result > 0.0);
        let expected = tiny * 2.0_f64.sqrt();
        assert!((result - expected).abs() / expected < 1e-15);
    }

    #[test]
    fn zabs_asymmetric() {
        // |1 + i| = sqrt(2)
        let z = Complex64::new(1.0, 1.0);
        let expected = 2.0_f64.sqrt();
        assert!((zabs(z) - expected).abs() < 1e-15);
    }

    // ── zdiv tests ──

    #[test]
    fn zdiv_simple() {
        // (1+0i) / (1+0i) = 1+0i
        let a = Complex64::new(1.0, 0.0);
        let b = Complex64::new(1.0, 0.0);
        let c = zdiv(a, b);
        assert!((c.re - 1.0).abs() < 1e-15);
        assert!(c.im.abs() < 1e-15);
    }

    #[test]
    fn zdiv_i_div_i() {
        // i / i = 1
        let a = Complex64::new(0.0, 1.0);
        let b = Complex64::new(0.0, 1.0);
        let c = zdiv(a, b);
        assert!((c.re - 1.0).abs() < 1e-15);
        assert!(c.im.abs() < 1e-15);
    }

    #[test]
    fn zdiv_known_result() {
        // (3+4i) / (1+2i) = (11/5 - 2/5 i)
        let a = Complex64::new(3.0, 4.0);
        let b = Complex64::new(1.0, 2.0);
        let c = zdiv(a, b);
        assert!((c.re - 2.2).abs() < 1e-14);
        assert!((c.im - (-0.4)).abs() < 1e-14);
    }

    #[test]
    fn zdiv_inverse() {
        // a / b * b ≈ a
        let a = Complex64::new(1.23456789, -9.87654321);
        let b = Complex64::new(0.314159, 2.71829);
        let c = zdiv(a, b);
        // Multiply back: c * b
        let recovered = Complex64::new(c.re * b.re - c.im * b.im, c.re * b.im + c.im * b.re);
        assert!((recovered.re - a.re).abs() < 1e-13);
        assert!((recovered.im - a.im).abs() < 1e-13);
    }

    #[test]
    fn zdiv_large_denominator_no_overflow() {
        let a = Complex64::new(1.0, 1.0);
        let b = Complex64::new(1.0e200, 1.0e200);
        let c = zdiv(a, b);
        // Expected: (1+i)/(1e200+1e200i) = 1/1e200
        let expected_re = 1.0e-200;
        assert!((c.re - expected_re).abs() / expected_re < 1e-14);
        assert!(c.im.abs() < 1e-214);
    }

    // ── f32 tests ──

    #[test]
    fn zabs_f32() {
        use num_complex::Complex32;
        let z = Complex32::new(3.0, 4.0);
        assert!((zabs(z) - 5.0).abs() < 1e-6);
    }

    #[test]
    fn zdiv_f32() {
        use num_complex::Complex32;
        let a = Complex32::new(3.0, 4.0);
        let b = Complex32::new(1.0, 2.0);
        let c = zdiv(a, b);
        assert!((c.re - 2.2).abs() < 1e-5);
        assert!((c.im - (-0.4)).abs() < 1e-5);
    }

    // ── sinpi tests ──

    #[test]
    fn sinpi_integers_are_zero() {
        for n in -5..=5 {
            let x = n as f64;
            assert_eq!(sinpi(x), 0.0, "sinpi({x}) should be exactly 0");
        }
    }

    #[test]
    fn sinpi_half_integers() {
        // sinpi(0.5) = 1, sinpi(1.5) = -1, sinpi(2.5) = 1, ...
        assert_eq!(sinpi(0.5_f64), 1.0);
        assert_eq!(sinpi(1.5_f64), -1.0);
        assert_eq!(sinpi(2.5_f64), 1.0);
        assert_eq!(sinpi(-0.5_f64), -1.0);
        assert_eq!(sinpi(-1.5_f64), 1.0);
    }

    #[test]
    fn sinpi_quarter() {
        let val = sinpi(0.25_f64);
        let expected = core::f64::consts::FRAC_1_SQRT_2;
        assert!((val - expected).abs() < 1e-15);
    }

    #[test]
    fn sinpi_general_values() {
        // sinpi(1/6) = sin(π/6) = 0.5
        let val = sinpi(1.0_f64 / 6.0);
        assert!((val - 0.5).abs() < 1e-15);

        // sinpi(1/3) = sin(π/3) = sqrt(3)/2
        let val = sinpi(1.0_f64 / 3.0);
        assert!((val - 3.0_f64.sqrt() / 2.0).abs() < 1e-15);
    }

    #[test]
    fn sinpi_large_argument() {
        // Large integer: sinpi(1e15) = 0
        assert_eq!(sinpi(1e15_f64), 0.0);
        // Large half-integer: sinpi(1e15 + 0.5) = ±1
        assert!(sinpi(1e15_f64 + 0.5).abs() == 1.0);
    }

    #[test]
    fn sinpi_f32() {
        assert_eq!(sinpi(0.0_f32), 0.0);
        assert_eq!(sinpi(0.5_f32), 1.0);
        assert_eq!(sinpi(1.0_f32), 0.0);
        assert_eq!(sinpi(1.5_f32), -1.0);
    }

    // ── cospi tests ──

    #[test]
    fn cospi_integers() {
        // cospi(0) = 1, cospi(1) = -1, cospi(2) = 1, ...
        assert_eq!(cospi(0.0_f64), 1.0);
        assert_eq!(cospi(1.0_f64), -1.0);
        assert_eq!(cospi(2.0_f64), 1.0);
        assert_eq!(cospi(-1.0_f64), -1.0);
        assert_eq!(cospi(-2.0_f64), 1.0);
    }

    #[test]
    fn cospi_half_integers_are_zero() {
        for n in -5..=5 {
            let x = n as f64 + 0.5;
            assert_eq!(cospi(x), 0.0, "cospi({x}) should be exactly 0");
        }
    }

    #[test]
    fn cospi_quarter() {
        let val = cospi(0.25_f64);
        let expected = core::f64::consts::FRAC_1_SQRT_2;
        assert!((val - expected).abs() < 1e-15);
    }

    #[test]
    fn cospi_general_values() {
        // cospi(1/3) = cos(π/3) = 0.5
        let val = cospi(1.0_f64 / 3.0);
        assert!((val - 0.5).abs() < 1e-15);

        // cospi(1/6) = cos(π/6) = sqrt(3)/2
        let val = cospi(1.0_f64 / 6.0);
        assert!((val - 3.0_f64.sqrt() / 2.0).abs() < 1e-15);
    }

    #[test]
    fn cospi_large_argument() {
        // Large even integer: cospi(1e15) = ±1 (1e15 is even)
        assert!(cospi(1e15_f64).abs() == 1.0);
        // Large half-integer: cospi(1e15 + 0.5) = 0
        assert_eq!(cospi(1e15_f64 + 0.5), 0.0);
    }

    #[test]
    fn cospi_f32() {
        assert_eq!(cospi(0.0_f32), 1.0);
        assert_eq!(cospi(0.5_f32), 0.0);
        assert_eq!(cospi(1.0_f32), -1.0);
        assert_eq!(cospi(1.5_f32), 0.0);
    }
}
