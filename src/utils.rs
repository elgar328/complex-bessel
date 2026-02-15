//! Overflow-safe complex arithmetic utilities.
//!
//! Translations of Amos ZABS and ZDIV from TOMS 644.

use num_complex::Complex;

use crate::machine::BesselFloat;

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
        let b = Complex64::new(0.314159, 2.71828);
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
}
