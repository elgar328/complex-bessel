//! Complex hyperbolic sine and cosine.
//!
//! Translation of Fortran ZSHCH from TOMS 644 (zbsubs.f lines 3082-3103).

use num_complex::Complex;

use crate::machine::BesselFloat;

/// Compute sinh(z) and cosh(z) for complex z = x + iy.
///
/// Returns `(sinh(z), cosh(z))` using the identities:
///   sinh(x+iy) = sinh(x)cos(y) + i·cosh(x)sin(y)
///   cosh(x+iy) = cosh(x)cos(y) + i·sinh(x)sin(y)
///
/// Equivalent to Fortran ZSHCH in TOMS 644.
#[inline]
pub(crate) fn zshch<T: BesselFloat>(z: Complex<T>) -> (Complex<T>, Complex<T>) {
    let sh = z.re.sinh();
    let ch = z.re.cosh();
    let sn = z.im.sin();
    let cn = z.im.cos();
    let csh = Complex::new(sh * cn, ch * sn);
    let cch = Complex::new(ch * cn, sh * sn);
    (csh, cch)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    #[test]
    fn shch_real() {
        // For real z: sinh(z) and cosh(z) should be real
        let z = Complex64::new(1.0, 0.0);
        let (csh, cch) = zshch(z);
        assert!((csh.re - 1.0_f64.sinh()).abs() < 1e-15);
        assert!(csh.im.abs() < 1e-15);
        assert!((cch.re - 1.0_f64.cosh()).abs() < 1e-15);
        assert!(cch.im.abs() < 1e-15);
    }

    #[test]
    fn shch_imaginary() {
        // sinh(iy) = i·sin(y), cosh(iy) = cos(y)
        let y = 1.5;
        let z = Complex64::new(0.0, y);
        let (csh, cch) = zshch(z);
        assert!(csh.re.abs() < 1e-15);
        assert!((csh.im - y.sin()).abs() < 1e-15);
        assert!((cch.re - y.cos()).abs() < 1e-15);
        assert!(cch.im.abs() < 1e-15);
    }

    #[test]
    fn shch_identity() {
        // cosh²(z) - sinh²(z) = 1
        let z = Complex64::new(1.5, 2.3);
        let (csh, cch) = zshch(z);
        let lhs = cch * cch - csh * csh;
        assert!((lhs.re - 1.0).abs() < 1e-13);
        assert!(lhs.im.abs() < 1e-13);
    }

    #[test]
    fn shch_zero() {
        let z = Complex64::new(0.0, 0.0);
        let (csh, cch) = zshch(z);
        assert!(csh.re.abs() < 1e-15);
        assert!(csh.im.abs() < 1e-15);
        assert!((cch.re - 1.0).abs() < 1e-15);
        assert!(cch.im.abs() < 1e-15);
    }
}
