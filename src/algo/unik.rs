//! Region 1 uniform asymptotic expansion parameter computation.
//!
//! Translation of Fortran ZUNIK from TOMS 644 / SLATEC (zbsubs.f lines 4806-5016).
//! Computes parameters PHI, ZETA1, ZETA2, SUM for the uniform asymptotic
//! expansions of the I and K functions.

// Exact Fortran constants — preserve verbatim.
#![allow(clippy::excessive_precision)]

use num_complex::Complex;

use crate::machine::BesselFloat;
use crate::types::{IkFlag, SumOption};
use crate::utils::zdiv;

// CON(1) = 1/sqrt(2*pi), CON(2) = sqrt(pi/2)
// Fortran zbsubs.f lines 4837-4838
const CON: [f64; 2] = [3.98942280401432678e-01, 1.25331413731550025e+00];

// Asymptotic expansion coefficients C(1)..C(120).
// Fortran zbsubs.f lines 4839-4914
#[rustfmt::skip]
const C_COEFFS: [f64; 120] = [
    1.00000000000000000e+00,   -2.08333333333333333e-01,
    1.25000000000000000e-01,    3.34201388888888889e-01,
   -4.01041666666666667e-01,    7.03125000000000000e-02,
   -1.02581259645061728e+00,    1.84646267361111111e+00,
   -8.91210937500000000e-01,    7.32421875000000000e-02,
    4.66958442342624743e+00,   -1.12070026162229938e+01,
    8.78912353515625000e+00,   -2.36408691406250000e+00,
    1.12152099609375000e-01,   -2.82120725582002449e+01,
    8.46362176746007346e+01,   -9.18182415432400174e+01,
    4.25349987453884549e+01,   -7.36879435947963170e+00,
    2.27108001708984375e-01,    2.12570130039217123e+02,
   -7.65252468141181642e+02,    1.05999045252799988e+03,
   -6.99579627376132541e+02,    2.18190511744211590e+02,
   -2.64914304869515555e+01,    5.72501420974731445e-01,
   -1.91945766231840700e+03,    8.06172218173730938e+03,
   -1.35865500064341374e+04,    1.16553933368645332e+04,
   -5.30564697861340311e+03,    1.20090291321635246e+03,
   -1.08090919788394656e+02,    1.72772750258445740e+00,
    2.02042913309661486e+04,   -9.69805983886375135e+04,
    1.92547001232531532e+05,   -2.03400177280415534e+05,
    1.22200464983017460e+05,   -4.11926549688975513e+04,
    7.10951430248936372e+03,   -4.93915304773088012e+02,
    6.07404200127348304e+00,   -2.42919187900551333e+05,
    1.31176361466297720e+06,   -2.99801591853810675e+06,
    3.76327129765640400e+06,   -2.81356322658653411e+06,
    1.26836527332162478e+06,   -3.31645172484563578e+05,
    4.52187689813627263e+04,   -2.49983048181120962e+03,
    2.43805296995560639e+01,    3.28446985307203782e+06,
   -1.97068191184322269e+07,    5.09526024926646422e+07,
   -7.41051482115326577e+07,    6.63445122747290267e+07,
   -3.75671766607633513e+07,    1.32887671664218183e+07,
   -2.78561812808645469e+06,    3.08186404612662398e+05,
   -1.38860897537170405e+04,    1.10017140269246738e+02,
   -4.93292536645099620e+07,    3.25573074185765749e+08,
   -9.39462359681578403e+08,    1.55359689957058006e+09,
   -1.62108055210833708e+09,    1.10684281682301447e+09,
   -4.95889784275030309e+08,    1.42062907797533095e+08,
   -2.44740627257387285e+07,    2.24376817792244943e+06,
   -8.40054336030240853e+04,    5.51335896122020586e+02,
    8.14789096118312115e+08,   -5.86648149205184723e+09,
    1.86882075092958249e+10,   -3.46320433881587779e+10,
    4.12801855797539740e+10,   -3.30265997498007231e+10,
    1.79542137311556001e+10,   -6.56329379261928433e+09,
    1.55927986487925751e+09,   -2.25105661889415278e+08,
    1.73951075539781645e+07,   -5.49842327572288687e+05,
    3.03809051092238427e+03,   -1.46792612476956167e+10,
    1.14498237732025810e+11,   -3.99096175224466498e+11,
    8.19218669548577329e+11,   -1.09837515608122331e+12,
    1.00815810686538209e+12,   -6.45364869245376503e+11,
    2.87900649906150589e+11,   -8.78670721780232657e+10,
    1.76347306068349694e+10,   -2.16716498322379509e+09,
    1.43157876718888981e+08,   -3.87183344257261262e+06,
    1.82577554742931747e+04,    2.86464035717679043e+11,
   -2.40629790002850396e+12,    9.10934118523989896e+12,
   -2.05168994109344374e+13,    3.05651255199353206e+13,
   -3.16670885847851584e+13,    2.33483640445818409e+13,
   -1.23204913055982872e+13,    4.61272578084913197e+12,
   -1.19655288019618160e+12,    2.05914503232410016e+11,
   -2.18229277575292237e+10,    1.24700929351271032e+09,
   -2.91883881222208134e+07,    1.18838426256783253e+05,
];

/// Cached workspace for ZUNIK.
///
/// When a previous call to `zunik` computed the coefficient workspace for
/// a given (z, fnu) pair, the cache can be reused on subsequent calls with
/// a different `ikflg` to avoid redundant computation.
#[derive(Clone, Copy, Debug)]
pub(crate) struct UnikCache<T: BesselFloat> {
    /// Number of terms to sum (Fortran INIT, range 2..=15, or 0 if not valid).
    pub init: usize,
    /// Coefficient workspace (16 elements, 0-based).
    /// [0..14] = asymptotic coefficients, [15] = sqrt(1/(fnu*sr)).
    pub cwrk: [Complex<T>; 16],
    /// Cached zeta1 value (invariant across ikflg changes).
    pub zeta1: Complex<T>,
    /// Cached zeta2 value (invariant across ikflg changes).
    pub zeta2: Complex<T>,
}

/// Output of the ZUNIK parameter computation.
#[derive(Debug, Clone, Copy)]
pub(crate) struct UnikOutput<T: BesselFloat> {
    pub phi: Complex<T>,
    pub zeta1: Complex<T>,
    pub zeta2: Complex<T>,
    pub sum: Complex<T>,
    pub cache: UnikCache<T>,
}

/// Compute the sum from cached coefficients (Fortran label 40).
///
/// For I function: straight sum of cwrk[0..init].
/// For K function: alternating sum of cwrk[0..init].
#[inline]
fn compute_sum<T: BesselFloat>(ikflg: IkFlag, cache: &UnikCache<T>) -> UnikOutput<T> {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);

    let (sum, con_idx) = match ikflg {
        IkFlag::I => {
            // I function: straight sum (Fortran lines 4988-4993)
            let mut sr = czero;
            for i in 0..cache.init {
                sr = sr + cache.cwrk[i];
            }
            (sr, 0usize)
        }
        IkFlag::K => {
            // K function: alternating sum (Fortran lines 5003-5010)
            let mut sr = czero;
            let mut tr = one;
            for i in 0..cache.init {
                sr = sr + cache.cwrk[i] * tr;
                tr = -tr;
            }
            (sr, 1usize)
        }
    };

    let con_val = T::from_f64(CON[con_idx]);
    let phi = cache.cwrk[15] * con_val;

    UnikOutput {
        phi,
        zeta1: cache.zeta1,
        zeta2: cache.zeta2,
        sum,
        cache: *cache,
    }
}

/// Compute parameters for the uniform asymptotic expansion of I and K functions.
///
/// Equivalent to Fortran ZUNIK in TOMS 644 (zbsubs.f lines 4806-5016).
///
/// # Parameters
/// - `zr`: complex argument (with Re(z) >= 0 for standard use)
/// - `fnu`: order ν > 0
/// - `ikflg`: `IkFlag::I` for I function, `IkFlag::K` for K function
/// - `ipmtr`: `SumOption::Full` = compute all, `SumOption::SkipSum` = phi/zeta only
/// - `tol`: tolerance for convergence
/// - `cache`: optional cache from a previous call with the same (z, fnu)
pub(crate) fn zunik<T: BesselFloat>(
    zr: Complex<T>,
    fnu: T,
    ikflg: IkFlag,
    ipmtr: SumOption,
    tol: T,
    cache: Option<UnikCache<T>>,
) -> UnikOutput<T> {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);
    let cone = Complex::from(one);

    // If cache is valid, skip initialization (Fortran: IF(INIT.NE.0) GO TO 40)
    if let Some(c) = cache {
        if c.init > 0 {
            return compute_sum(ikflg, &c);
        }
    }

    // ── Initialize all variables (Fortran lines 4920-4982) ──

    let rfn = one / fnu;

    // Overflow test: z/fnu too small (Fortran lines 4924-4933)
    let test = T::MACH_TINY * T::from_f64(1.0e3);
    let ac = fnu * test;
    if zr.re.abs() <= ac && zr.im.abs() <= ac {
        // Early return with special values
        let zeta1 = Complex::new(T::from_f64(2.0) * test.ln().abs() + fnu, zero);
        let zeta2 = Complex::new(fnu, zero);
        return UnikOutput {
            phi: cone,
            zeta1,
            zeta2,
            sum: czero,
            cache: UnikCache {
                init: 0,
                cwrk: [czero; 16],
                zeta1,
                zeta2,
            },
        };
    }

    // t = z / fnu (Fortran lines 4935-4936)
    let t = zr * rfn;

    // s = 1 + t^2 (Fortran lines 4937-4938)
    let s = cone + t * t;

    // sr = sqrt(s) = sqrt(1 + (z/fnu)^2) (Fortran line 4939)
    let sr = s.sqrt();

    // zn = (1 + sr) / t (Fortran lines 4940-4942)
    let zn = zdiv(cone + sr, t);

    // zeta1 = fnu * log(zn) (Fortran lines 4943-4945)
    let ln_zn = zn.ln();
    let zeta1 = ln_zn * fnu;

    // zeta2 = fnu * sr (Fortran lines 4946-4947)
    let zeta2 = sr * fnu;

    // sr_new = 1/(fnu*sr) (Fortran lines 4948-4950)
    let t_inv = zdiv(cone, sr);
    let sr_new = t_inv * rfn;

    // cwrk[15] = sqrt(1/(fnu*sr)) (Fortran line 4951)
    let mut cwrk = [czero; 16];
    cwrk[15] = sr_new.sqrt();

    // phi = cwrk[15] * CON[ikflg-1] (Fortran lines 4952-4953)
    let con_val = T::from_f64(match ikflg {
        IkFlag::I => CON[0],
        IkFlag::K => CON[1],
    });
    let phi = cwrk[15] * con_val;

    // If ipmtr == SkipSum, return phi/zeta1/zeta2 only (Fortran line 4954)
    if ipmtr == SumOption::SkipSum {
        return UnikOutput {
            phi,
            zeta1,
            zeta2,
            sum: czero,
            cache: UnikCache {
                init: 0,
                cwrk,
                zeta1,
                zeta2,
            },
        };
    }

    // ── Coefficient computation (Fortran lines 4955-4982) ──

    // t2 = 1/s = 1/(1 + (z/fnu)^2) (Fortran line 4955)
    let t2 = zdiv(cone, s);

    // cwrk[0] = 1.0 (Fortran lines 4956-4957)
    cwrk[0] = cone;

    // crfn accumulates powers of sr_new = 1/(fnu*sr)
    let mut crfn = cone;
    let mut ac_val = one;
    let mut l: usize = 0; // 0-based C_COEFFS index (Fortran L starts at 1)
    let mut init: usize = 15; // default: all 15 terms (Fortran K=15)

    // Fortran DO 20 K=2,15 → Rust k=1..=14 (cwrk[1]..cwrk[14])
    #[allow(clippy::needless_range_loop)]
    for k in 1..=14 {
        // Inner polynomial evaluation: s = sum_{j=0}^{k} C[l+j+1] * t2^(k-j)
        // Fortran DO 10 J=1,K (K=k+1 iterations)
        let mut s_sum = czero;
        for _j in 0..=k {
            l += 1;
            let c_val = T::from_f64(C_COEFFS[l]);
            // s = s * t2 + c[l] (Fortran lines 4967-4969)
            s_sum = s_sum * t2 + Complex::new(c_val, zero);
        }

        // crfn = crfn * sr_new (Fortran lines 4971-4973)
        crfn = crfn * sr_new;

        // cwrk[k] = crfn * s_sum (Fortran lines 4974-4975)
        cwrk[k] = crfn * s_sum;

        // Convergence test (Fortran lines 4976-4978)
        ac_val = ac_val * rfn;
        let test_val = cwrk[k].re.abs() + cwrk[k].im.abs();
        if ac_val < tol && test_val < tol {
            init = k + 1; // Fortran INIT = K, where K = k+1
            break;
        }
    }

    // Compute sum based on ikflg (Fortran labels 40-70)
    let cache = UnikCache {
        init,
        cwrk,
        zeta1,
        zeta2,
    };
    compute_sum(ikflg, &cache)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{IkFlag, SumOption};
    use num_complex::Complex64;

    const TOL: f64 = 2.220446049250313e-16;

    #[test]
    fn zunik_cache_reuse() {
        // First call with IkFlag::K, then reuse cache with IkFlag::I
        let zr = Complex64::new(3.0, 1.0);
        let result_k = zunik(zr, 15.0, IkFlag::K, SumOption::Full, TOL, None);
        let result_i = zunik(
            zr,
            15.0,
            IkFlag::I,
            SumOption::Full,
            TOL,
            Some(result_k.cache),
        );

        // zeta1, zeta2 should be identical (from cache)
        assert_eq!(result_i.zeta1.re, result_k.zeta1.re);
        assert_eq!(result_i.zeta1.im, result_k.zeta1.im);
        assert_eq!(result_i.zeta2.re, result_k.zeta2.re);
        assert_eq!(result_i.zeta2.im, result_k.zeta2.im);

        // phi should use CON(1) for I function
        let con1 = CON[0];
        let con2 = CON[1];
        // When z is complex, phi has imaginary part too, so check via norms
        let phi_i_mag =
            (result_i.phi.re * result_i.phi.re + result_i.phi.im * result_i.phi.im).sqrt();
        let phi_k_mag =
            (result_k.phi.re * result_k.phi.re + result_k.phi.im * result_k.phi.im).sqrt();
        assert!(
            (phi_i_mag / phi_k_mag - con1 / con2).abs() < 1e-12,
            "phi magnitude ratio: {}",
            phi_i_mag / phi_k_mag
        );

        // Sum should differ (I: straight sum, K: alternating sum)
        let result_i_fresh = zunik(zr, 15.0, IkFlag::I, SumOption::Full, TOL, None);
        assert!((result_i.sum.re - result_i_fresh.sum.re).abs() < 1e-14);
    }

    #[test]
    fn zunik_ipmtr_1() {
        // SumOption::SkipSum: compute only phi, zeta1, zeta2 (no sum)
        let zr = Complex64::new(4.0, 2.0);
        let result = zunik(zr, 12.0, IkFlag::I, SumOption::SkipSum, TOL, None);

        // phi should be nonzero
        assert!(result.phi.re.abs() > 1e-10 || result.phi.im.abs() > 1e-10);

        // cache.init should be 0 (no coefficient computation)
        assert_eq!(result.cache.init, 0);
    }

    #[test]
    fn zunik_overflow_precheck() {
        // Very small z relative to fnu → overflow pre-check path
        let zr = Complex64::new(1e-310, 1e-310);
        let result = zunik(zr, 1.0, IkFlag::I, SumOption::Full, TOL, None);

        // Should return special values
        assert_eq!(result.phi.re, 1.0);
        assert_eq!(result.phi.im, 0.0);
        assert_eq!(result.zeta2.re, 1.0);
        assert!(result.zeta1.re > result.zeta2.re);
        assert_eq!(result.cache.init, 0);
    }
}
