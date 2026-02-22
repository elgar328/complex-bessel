//! Natural logarithm of the Gamma function.
//!
//! Translation of Fortran DGAMLN from TOMS 644 (zbsubs.f lines 4489-4677).
//! Uses table lookup for integer arguments 1..100 and Stirling's
//! asymptotic expansion for general z > 0.

// Constants and tables are exact Fortran values — preserve verbatim.
#![allow(clippy::excessive_precision)]
#![allow(clippy::approx_constant)]

use crate::algo::constants::R1M5;
use crate::machine::BesselFloat;
use crate::types::Error;

/// ln(2π), used in Stirling's formula.
const CON: f64 = 1.83787706640934548;

/// ln(Γ(n)) for n = 1, 2, ..., 100.
///
/// Used for fast table lookup when the argument is a positive integer ≤ 100.
/// Γ(1) = Γ(2) = 1, so ln(Γ(1)) = ln(Γ(2)) = 0.
#[rustfmt::skip]
const GLN_TABLE: [f64; 100] = [
    0.00000000000000000e+00,  0.00000000000000000e+00,  // Γ(1)=1, Γ(2)=1
    6.93147180559945309e-01,  1.79175946922805500e+00,  // Γ(3)=2, Γ(4)=6
    3.17805383034794562e+00,  4.78749174278204599e+00,
    6.57925121201010100e+00,  8.52516136106541430e+00,
    1.06046029027452502e+01,  1.28018274800814696e+01,
    1.51044125730755153e+01,  1.75023078458738858e+01,
    1.99872144956618861e+01,  2.25521638531234229e+01,
    2.51912211827386815e+01,  2.78992713838408916e+01,
    3.06718601060806728e+01,  3.35050734501368889e+01,
    3.63954452080330536e+01,  3.93398841871994940e+01,
    4.23356164607534850e+01,  4.53801388984769080e+01,
    4.84711813518352239e+01,  5.16066755677643736e+01,
    5.47847293981123192e+01,  5.80036052229805199e+01,
    6.12617017610020020e+01,  6.45575386270063311e+01,
    6.78897431371815350e+01,  7.12570389671680090e+01,
    7.46582363488301644e+01,  7.80922235533153106e+01,
    8.15579594561150372e+01,  8.50544670175815174e+01,
    8.85808275421976788e+01,  9.21361756036870925e+01,
    9.57196945421432025e+01,  9.93306124547874269e+01,
    1.02968198614513813e+02,  1.06631760260643459e+02,
    1.10320639714757395e+02,  1.14034211781461703e+02,
    1.17771881399745072e+02,  1.21533081515438634e+02,
    1.25317271149356895e+02,  1.29123933639127215e+02,
    1.32952575035616310e+02,  1.36802722637326368e+02,
    1.40673923648234259e+02,  1.44565743946344886e+02,
    1.48477766951773032e+02,  1.52409592584497358e+02,
    1.56360836303078785e+02,  1.60331128216630907e+02,
    1.64320112263195181e+02,  1.68327445448427652e+02,
    1.72352797139162802e+02,  1.76395848406997352e+02,
    1.80456291417543771e+02,  1.84533828861449491e+02,
    1.88628173423671591e+02,  1.92739047287844902e+02,
    1.96866181672889994e+02,  2.01009316399281527e+02,
    2.05168199482641199e+02,  2.09342586752536836e+02,
    2.13532241494563261e+02,  2.17736934113954227e+02,
    2.21956441819130334e+02,  2.26190548323727593e+02,
    2.30439043565776952e+02,  2.34701723442818268e+02,
    2.38978389561834323e+02,  2.43268849002982714e+02,
    2.47572914096186884e+02,  2.51890402209723194e+02,
    2.56221135550009525e+02,  2.60564940971863209e+02,
    2.64921649798552801e+02,  2.69291097651019823e+02,
    2.73673124285693704e+02,  2.78067573440366143e+02,
    2.82474292687630396e+02,  2.86893133295426994e+02,
    2.91323950094270308e+02,  2.95766601350760624e+02,
    3.00220948647014132e+02,  3.04686856765668715e+02,
    3.09164193580146922e+02,  3.13652829949879062e+02,
    3.18152639620209327e+02,  3.22663499126726177e+02,
    3.27185287703775217e+02,  3.31717887196928473e+02,
    3.36261181979198477e+02,  3.40815058870799018e+02,
    3.45379407062266854e+02,  3.49954118040770237e+02,
    3.54539085519440809e+02,  3.59134205369575399e+02,
];

/// Coefficients of the asymptotic expansion for ln(Γ(z)).
///
/// These are related to the Bernoulli numbers B_{2k}:
///   CF(k) = B_{2k} / (2k * (2k-1))
/// for k = 1, 2, ..., 22.
#[rustfmt::skip]
const CF_TABLE: [f64; 22] = [
     8.33333333333333333e-02,   // B2/(1*2)    = 1/12
    -2.77777777777777778e-03,   // B4/(3*4)    = -1/360
     7.93650793650793651e-04,   // B6/(5*6)
    -5.95238095238095238e-04,   // B8/(7*8)
     8.41750841750841751e-04,   // B10/(9*10)
    -1.91752691752691753e-03,   // B12/(11*12)
     6.41025641025641026e-03,   // B14/(13*14)
    -2.95506535947712418e-02,   // B16/(15*16)
     1.79644372368830573e-01,   // B18/(17*18)
    -1.39243221690590112e+00,   // B20/(19*20)
     1.34028640441683920e+01,   // B22/(21*22)
    -1.56848284626002017e+02,
     2.19310333333333333e+03,
    -3.61087712537249894e+04,
     6.91472268851313067e+05,
    -1.52382215394074162e+07,
     3.82900751391414141e+08,
    -1.08822660357843911e+10,
     3.47320283765002252e+11,
    -1.23696021422692745e+13,
     4.88788064793079335e+14,
    -2.13203339609193739e+16,
];

/// Compute ln(Γ(z)) for z > 0.
///
/// Algorithm:
/// 1. For positive integers 1..100: table lookup (exact).
/// 2. For general z: Stirling's asymptotic expansion with recursion
///    reduction when z is below the convergence threshold ZMIN.
///
/// Equivalent to Fortran DGAMLN in SLATEC / TOMS 644.
pub(crate) fn gamln<T: BesselFloat>(z: T) -> Result<T, Error> {
    let zero = T::zero();
    let one = T::one();

    if z <= zero {
        return Err(Error::InvalidInput);
    }

    // Table lookup for positive integers 1..100
    // Fortran: NZ = Z; FZ = Z - NZ; IF (FZ.GT.0) GO TO 10; IF (NZ.GT.100) GO TO 10
    if z <= T::from_f64(101.0) {
        // Safety: z is checked > 0 and <= 101 here
        let nz = z.floor().to_i32().unwrap();
        let fz = z - T::from_f64(nz as f64);
        if fz == zero && (1..=100).contains(&nz) {
            return Ok(T::from_f64(GLN_TABLE[(nz - 1) as usize]));
        }
    }

    // D1MACH(4) = B^(2-T) = 2 * MACH_EPSILON for binary IEEE 754
    let two = T::from_f64(2.0);
    let wdtol = (two * T::MACH_EPSILON).max(T::from_f64(0.5e-18));

    // Compute ZMIN: minimum argument for asymptotic convergence
    // R1M5 = D1MACH(5) = log10(2)
    // RLN = R1M5 * I1MACH(14) = log10(2) * DIGITS
    let r1m5 = T::from_f64(R1M5);
    let rln = r1m5 * T::from_f64(T::MACH_DIGITS as f64);
    let fln = rln.min(T::from_f64(20.0)).max(T::from_f64(3.0)) - T::from_f64(3.0);
    let zm = T::from_f64(1.8) + T::from_f64(0.3875) * fln;
    // Safety: zm is a small positive value derived from log-based formula
    let mz = zm.to_i32().unwrap() + 1;
    let zmin = T::from_f64(mz as f64);

    // If z < ZMIN, push z upward using the recursion Γ(z+n) = z(z+1)...(z+n-1)·Γ(z)
    // Safety: z is checked > 0 at function entry and is finite
    let nz_trunc = z.floor().to_i32().unwrap(); // integer part of original z
    let (zdmy, zinc) = if z < zmin {
        let zinc_val = T::from_f64((mz - nz_trunc) as f64);
        (z + zinc_val, zinc_val)
    } else {
        (z, zero)
    };

    // Evaluate Stirling's asymptotic series:
    //   S(z) = Σ_{k=1}^{22} CF(k) * z^{-(2k-1)}
    let zp_init = one / zdmy;
    let t1 = T::from_f64(CF_TABLE[0]) * zp_init;
    let mut s = t1;

    if zp_init >= wdtol {
        let zsq = zp_init * zp_init;
        let tst = t1 * wdtol;
        let mut zp = zp_init;
        for cf in &CF_TABLE[1..] {
            zp = zp * zsq;
            let trm = T::from_f64(*cf) * zp;
            if trm.abs() < tst {
                break;
            }
            s = s + trm;
        }
    }

    // Stirling's formula: ln Γ(z) = z(ln z - 1) + 0.5(ln(2π) - ln z) + S(z)
    let half = T::from_f64(0.5);
    let con = T::from_f64(CON);

    if zinc == zero {
        // No recursion needed
        let tlg = z.ln();
        Ok(z * (tlg - one) + half * (con - tlg) + s)
    } else {
        // Undo recursion: ln Γ(z) = ln Γ(zdmy) - ln(z·(z+1)·...·(z+zinc-1))
        // Safety: zinc = mz - nz_trunc, a small non-negative integer
        let nz_int = zinc.to_i32().unwrap();
        let mut product = one;
        for i in 0..nz_int {
            product = product * (z + T::from_f64(i.into()));
        }
        let tlg = zdmy.ln();
        Ok(zdmy * (tlg - one) - product.ln() + half * (con - tlg) + s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gamln_integers() {
        // ln(Γ(1)) = 0
        assert_eq!(gamln(1.0_f64).unwrap(), 0.0);
        // ln(Γ(2)) = ln(1) = 0
        assert_eq!(gamln(2.0_f64).unwrap(), 0.0);
        // ln(Γ(3)) = ln(2) = 0.6931471805599453
        assert!((gamln(3.0_f64).unwrap() - 2.0_f64.ln()).abs() < 1e-15);
        // ln(Γ(4)) = ln(6) = 1.7917594692280550
        assert!((gamln(4.0_f64).unwrap() - 6.0_f64.ln()).abs() < 1e-14);
        // ln(Γ(7)) = ln(720) = 6.5792512120101010
        assert!((gamln(7.0_f64).unwrap() - 720.0_f64.ln()).abs() < 1e-13);
    }

    #[test]
    fn gamln_table_boundary() {
        // n=100: last table entry
        let val = gamln(100.0_f64).unwrap();
        assert!((val - 3.59134205369575399e+02).abs() < 1e-10);
        // n=101: should NOT use table (101 > 100), goes to asymptotic
        let val101 = gamln(101.0_f64).unwrap();
        // ln(Γ(101)) = ln(100!) ≈ 363.739...
        assert!((val101 - 363.73937555556347).abs() < 1e-9);
    }

    #[test]
    fn gamln_half_integer() {
        // Γ(0.5) = √π, so ln(Γ(0.5)) = 0.5*ln(π) ≈ 0.5723649429247
        let val = gamln(0.5_f64).unwrap();
        let expected = 0.5 * core::f64::consts::PI.ln();
        assert!((val - expected).abs() < 1e-14);
    }

    #[test]
    fn gamln_small_positive() {
        // Γ(0.1) ≈ 9.51351, ln(Γ(0.1)) ≈ 2.25271...
        let val = gamln(0.1_f64).unwrap();
        assert!((val - 2.2527126517342055).abs() < 1e-13);
    }

    #[test]
    fn gamln_large() {
        // ln(Γ(150)) = ln(149!) — tests asymptotic expansion (no recursion needed)
        // Reference: sum of ln(k) for k=1..149
        let val = gamln(150.0_f64).unwrap();
        let mut expected = 0.0_f64;
        for k in 1..150 {
            expected += (k as f64).ln();
        }
        assert!((val - expected).abs() / expected < 1e-14);
    }

    #[test]
    fn gamln_fractional_medium() {
        // Γ(5.5) = 4.5 * 3.5 * 2.5 * 1.5 * 0.5 * Γ(0.5)
        // ln(Γ(5.5)) = ln(4.5*3.5*2.5*1.5*0.5) + 0.5*ln(π)
        let val = gamln(5.5_f64).unwrap();
        let product: f64 = 4.5 * 3.5 * 2.5 * 1.5 * 0.5;
        let expected = product.ln() + 0.5 * core::f64::consts::PI.ln();
        assert!((val - expected).abs() < 1e-13);
    }

    #[test]
    fn gamln_negative_returns_error() {
        assert!(gamln(-1.0_f64).is_err());
        assert!(gamln(0.0_f64).is_err());
    }

    #[test]
    fn gamln_f32() {
        // Γ(3) = 2, ln(Γ(3)) = ln(2)
        let val = gamln(3.0_f32).unwrap();
        assert!((val - 2.0_f32.ln()).abs() < 1e-6);
    }
}
