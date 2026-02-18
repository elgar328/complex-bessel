//! Region 2 uniform asymptotic expansion for the I function.
//!
//! Translation of Fortran ZUNI2 from TOMS 644 / SLATEC (zbsubs.f lines 7045-7312).
//! Computes I(fnu,z) in the region |arg(z)| > pi/3 using the uniform
//! asymptotic expansion for J(fnu, z*exp(m*hpi*i)).

#![allow(clippy::too_many_arguments)]

use num_complex::Complex;

use crate::airy::zairy;
use crate::algo::uchk::zuchk;
use crate::algo::unhj::zunhj;
use crate::algo::uoik::zuoik;
use crate::machine::BesselFloat;
use crate::types::{AiryDerivative, IkFlag, Scaling, SumOption};
use crate::utils::{reciprocal_z, zabs};

use crate::algo::constants::{AIC, HPI};

/// CIP rotation table (Fortran lines 7077-7078)
/// CIPR + i*CIPI: {1, i, -1, -i}
const CIPR: [f64; 4] = [1.0, 0.0, -1.0, 0.0];
const CIPI: [f64; 4] = [0.0, 1.0, 0.0, -1.0];

/// Output of ZUNI2.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Uni2Output {
    /// Underflow count (number of zeroed trailing members).
    /// -1 indicates overflow.
    pub nz: i32,
    /// If nonzero, remaining items need a different method.
    pub nlast: i32,
}

/// Compute I(fnu,z) via Region 2 uniform asymptotic expansion.
///
/// Equivalent to Fortran ZUNI2 in TOMS 644 (zbsubs.f lines 7045-7312).
pub(crate) fn zuni2<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    y: &mut [Complex<T>],
    fnul: T,
    tol: T,
    elim: T,
    alim: T,
) -> Uni2Output {
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);

    let n = y.len();
    let mut nz: i32 = 0;
    let mut nd = n;
    let mut nlast: i32 = 0;

    for v in y.iter_mut() {
        *v = czero;
    }

    // ── 3-level scaling (Fortran lines 7090-7097) ──
    let cscl = one / tol;
    let crsc = tol;
    let cssr = [cscl, one, crsc];
    let csrr = [crsc, one, cscl];
    let bry0 = T::from_f64(1.0e3) * T::MACH_TINY / tol;

    // ── Rotate z to right half plane (Fortran lines 7102-7106) ──
    let mut znr = z.im; // ZNR = ZI
    let zni = -z.re; // ZNI = -ZR
    let zbr = z.re;
    let mut zbi = z.im;
    let mut cidi = -one; // CIDI = -1

    let inu = fnu.to_i32().unwrap() as usize;
    let ang = T::from_f64(HPI) * (fnu - T::from_f64(inu as f64));
    let c2_base = Complex::new(ang.cos(), ang.sin());

    // IN = (INU + N - 1) mod 4 (0-based)
    let in_idx = (inu + n - 1) % 4;
    let mut c2 = c2_base * Complex::new(T::from_f64(CIPR[in_idx]), T::from_f64(CIPI[in_idx]));

    if z.im <= zero {
        // Fortran lines 7119-7122
        znr = -znr;
        zbi = -zbi;
        cidi = -cidi;
        c2 = c2.conj();
    }

    // ── Check first member for overflow/underflow (Fortran lines 7127-7144) ──
    let fn_val = fnu.max(one);
    let result0 = zunhj(Complex::new(znr, zni), fn_val, SumOption::SkipSum, tol);

    let rs1 = if kode == Scaling::Exponential {
        let str = zbr + result0.zeta2.re;
        let sti = zbi + result0.zeta2.im;
        let rast = fn_val / zabs(Complex::new(str, sti));
        let str2 = str * rast * rast;
        -result0.zeta1.re + str2
    } else {
        -result0.zeta1.re + result0.zeta2.re
    };

    if rs1.abs() > elim {
        if rs1 > zero {
            return Uni2Output { nz: -1, nlast: 0 };
        }
        // All underflow
        nz = n as i32;
        for item in y.iter_mut() {
            *item = czero;
        }
        return Uni2Output { nz, nlast: 0 };
    }

    // ── Label 40: main computation loop (Fortran lines 7145-7263) ──
    'label40: loop {
        let nn = nd.min(2);
        let mut iflag: usize = 2;
        let mut cy = [czero; 2];

        let mut computed_ok = true;

        #[allow(clippy::needless_range_loop)]
        for i in 0..nn {
            let fn_val = fnu + T::from_f64((nd - 1 - i) as f64);
            let result = zunhj(Complex::new(znr, zni), fn_val, SumOption::Full, tol);

            let (s1r, s1i) = if kode == Scaling::Exponential {
                // KODE=2 path (Fortran lines 7151-7158)
                // NOTE: s1i += |zi| (absolute value!), NOT zi like ZUNI1
                let str = zbr + result.zeta2.re;
                let sti = zbi + result.zeta2.im;
                let rast = fn_val / zabs(Complex::new(str, sti));
                let str2 = str * rast * rast;
                let sti2 = -sti * rast * rast;
                (
                    -result.zeta1.re + str2,
                    -result.zeta1.im + sti2 + z.im.abs(),
                )
            } else {
                (
                    -result.zeta1.re + result.zeta2.re,
                    -result.zeta1.im + result.zeta2.im,
                )
            };

            // ── Overflow/underflow test (Fortran lines 7167-7203) ──
            let mut rs1 = s1r;
            if rs1.abs() > elim {
                if rs1 > zero {
                    return Uni2Output { nz: -1, nlast: 0 };
                }
                // Underflow (Fortran label 120)
                y[nd - 1] = czero;
                nz += 1;
                nd -= 1;
                if nd == 0 {
                    return Uni2Output { nz, nlast: 0 };
                }
                let nuf = zuoik(z, fnu, kode, IkFlag::I, &mut y[..nd], tol, elim, alim);
                if nuf < 0 {
                    return Uni2Output { nz: -1, nlast: 0 };
                }
                nd -= nuf as usize;
                nz += nuf;
                if nd == 0 {
                    return Uni2Output { nz, nlast: 0 };
                }
                let fn_check = fnu + T::from_f64((nd - 1) as f64);
                if fn_check >= fnul {
                    // Recalculate C2 (Fortran lines 7292-7296)
                    let in_idx2 = (inu + nd - 1) % 4;
                    c2 = c2_base
                        * Complex::new(T::from_f64(CIPR[in_idx2]), T::from_f64(CIPI[in_idx2]));
                    if z.im <= zero {
                        c2 = c2.conj();
                    }
                    computed_ok = false;
                    break;
                }
                nlast = nd as i32;
                return Uni2Output { nz, nlast };
            }

            if i == 0 {
                iflag = 2;
            }
            if rs1.abs() >= alim {
                // Refine test (Fortran lines 7175-7181)
                let aphi = zabs(result.phi);
                let aarg = zabs(result.arg);
                rs1 = rs1 + aphi.ln() - T::from_f64(0.25) * aarg.ln() - T::from_f64(AIC);
                if rs1.abs() > elim {
                    if rs1 > zero {
                        return Uni2Output { nz: -1, nlast: 0 };
                    }
                    y[nd - 1] = czero;
                    nz += 1;
                    nd -= 1;
                    if nd == 0 {
                        return Uni2Output { nz, nlast: 0 };
                    }
                    let nuf = zuoik(z, fnu, kode, IkFlag::I, &mut y[..nd], tol, elim, alim);
                    if nuf < 0 {
                        return Uni2Output { nz: -1, nlast: 0 };
                    }
                    nd -= nuf as usize;
                    nz += nuf;
                    if nd == 0 {
                        return Uni2Output { nz, nlast: 0 };
                    }
                    let fn_check = fnu + T::from_f64((nd - 1) as f64);
                    if fn_check >= fnul {
                        let in_idx2 = (inu + nd - 1) % 4;
                        c2 = c2_base
                            * Complex::new(T::from_f64(CIPR[in_idx2]), T::from_f64(CIPI[in_idx2]));
                        if z.im <= zero {
                            c2 = c2.conj();
                        }
                        computed_ok = false;
                        break;
                    }
                    nlast = nd as i32;
                    return Uni2Output { nz, nlast };
                }
                if i == 0 {
                    iflag = 1;
                }
                if rs1 < zero {
                    // iflag stays
                } else if i == 0 {
                    iflag = 3;
                }
            }

            // ── Scale S1 and compute S2 (Fortran lines 7187-7216) ──
            // S2 = PHI * (Ai*ASUM + Ai'*BSUM)
            let (ai, _nai) = zairy(result.arg, AiryDerivative::Value, Scaling::Exponential)
                .unwrap_or((czero, 0));
            let (dai, _ndai) = zairy(result.arg, AiryDerivative::Derivative, Scaling::Exponential)
                .unwrap_or((czero, 0));

            let s2_airy = result.phi * (ai * result.asum + dai * result.bsum);

            let str_exp = s1r.exp() * cssr[iflag - 1];
            let s1_scaled = Complex::new(str_exp * s1i.cos(), str_exp * s1i.sin());
            let mut s2_val = s2_airy * s1_scaled;

            if iflag == 1 && zuchk(s2_val, bry0, tol) {
                // Underflow (label 120)
                y[nd - 1] = czero;
                nz += 1;
                nd -= 1;
                if nd == 0 {
                    return Uni2Output { nz, nlast: 0 };
                }
                let nuf = zuoik(z, fnu, kode, IkFlag::I, &mut y[..nd], tol, elim, alim);
                if nuf < 0 {
                    return Uni2Output { nz: -1, nlast: 0 };
                }
                nd -= nuf as usize;
                nz += nuf;
                if nd == 0 {
                    return Uni2Output { nz, nlast: 0 };
                }
                let fn_check = fnu + T::from_f64((nd - 1) as f64);
                if fn_check >= fnul {
                    let in_idx2 = (inu + nd - 1) % 4;
                    c2 = c2_base
                        * Complex::new(T::from_f64(CIPR[in_idx2]), T::from_f64(CIPI[in_idx2]));
                    if z.im <= zero {
                        c2 = c2.conj();
                    }
                    computed_ok = false;
                    break;
                }
                nlast = nd as i32;
                return Uni2Output { nz, nlast };
            }

            // Conjugation for ZI <= 0 (Fortran line 7205)
            if z.im <= zero {
                s2_val = s2_val.conj();
            }

            // Phase rotation by C2 (Fortran lines 7206-7208)
            s2_val = s2_val * c2;

            // Store (Fortran lines 7209-7213)
            cy[i] = s2_val;
            let j_idx = nd - 1 - i;
            y[j_idx] = s2_val * csrr[iflag - 1];

            // C2 *= i*CIDI (Fortran lines 7214-7216)
            c2 = c2 * Complex::new(zero, cidi);
        }

        if !computed_ok {
            continue 'label40;
        }

        // ── Forward recurrence for remaining terms (Fortran lines 7218-7263) ──
        if nd <= 2 {
            return Uni2Output { nz, nlast };
        }

        let rz = reciprocal_z(z);

        let bry1 = one / bry0;
        let bry2 = T::MACH_HUGE;

        let mut s1 = cy[0];
        let mut s2 = cy[1];
        let mut c1r = csrr[iflag - 1];
        let mut ascle = if iflag == 1 {
            bry0
        } else if iflag == 2 {
            bry1
        } else {
            bry2
        };

        let mut k = nd - 2;
        let mut fn_rec = T::from_f64(k as f64);

        for _i in 2..nd {
            let prev = s2;
            let cfn = fnu + fn_rec;
            s2 = s1 + rz * prev * cfn;
            s1 = prev;
            let c2_scaled = s2 * c1r;
            y[k] = c2_scaled;
            k -= 1;
            fn_rec = fn_rec - one;

            if iflag >= 3 {
                continue;
            }
            let c2m = c2_scaled.re.abs().max(c2_scaled.im.abs());
            if c2m <= ascle {
                continue;
            }
            iflag += 1;
            ascle = if iflag == 2 { bry1 } else { bry2 };
            s1 = s1 * c1r;
            s2 = c2_scaled;
            s1 = s1 * cssr[iflag - 1];
            s2 = s2 * cssr[iflag - 1];
            c1r = csrr[iflag - 1];
        }

        return Uni2Output { nz, nlast };
    }
}
