//! Region 1 uniform asymptotic expansion for the K function + analytic continuation.
//!
//! Translation of Fortran ZUNK1 from TOMS 644 / SLATEC (zbsubs.f lines 5732-6158).
//! Computes K(fnu,z) and its analytic continuation from the right half plane
//! to the left half plane using the uniform asymptotic expansion.

#![allow(clippy::too_many_arguments)]

use num_complex::Complex;

use crate::algo::constants::PI;
use crate::algo::s1s2::zs1s2;
use crate::algo::uchk::zuchk;
use crate::algo::unik::{UnikCache, zunik};
use crate::machine::BesselFloat;
use crate::types::{IkFlag, Scaling, SumOption};
use crate::utils::{mul_add, reciprocal_z, zabs};

/// Compute K(fnu,z) via Region 1 uniform asymptotic expansion + analytic continuation.
///
/// Equivalent to Fortran ZUNK1 in TOMS 644 (zbsubs.f lines 5732-6158).
///
/// # Parameters
/// - `z`: complex argument
/// - `fnu`: starting order v >= 0
/// - `kode`: scaling mode
/// - `mr`: analytic continuation direction (0 = none, +/-1 = left half plane)
/// - `y`: output slice for sequence members
/// - `tol`, `elim`, `alim`: machine-derived thresholds
///
/// # Returns
/// `nz` where nz = -1 indicates overflow.
pub(crate) fn zunk1<T: BesselFloat>(
    z: Complex<T>,
    fnu: T,
    kode: Scaling,
    mr: i32,
    y: &mut [Complex<T>],
    tol: T,
    elim: T,
    alim: T,
) -> i32 {
    let n = y.len();
    let zero = T::zero();
    let one = T::one();
    let czero = Complex::new(zero, zero);
    let pi_t = T::from_f64(PI);

    y.fill(czero);
    let mut kdflg: usize = 1;
    let mut nz: i32 = 0;

    // -- 3-level scaling (Fortran lines 5769-5779) --
    let cscl = one / tol;
    let crsc = tol;
    let cssr = [cscl, one, crsc];
    let csrr = [crsc, one, cscl];
    let bry = [
        T::from_f64(1.0e3) * T::MACH_TINY / tol,
        one / (T::from_f64(1.0e3) * T::MACH_TINY / tol),
        T::MACH_HUGE,
    ];

    // -- Reflect to right half plane (Fortran lines 5780-5784) --
    let zr = if z.re >= zero { z } else { -z };

    // -- Phase 1: K function computation (Fortran lines 5786-5865) --
    let mut j: usize = 1; // Fortran J starts at 2, first flip gives 1 -> Rust: start 1, first flip gives 0
    let mut phi_arr = [czero; 2];
    let mut zeta1_arr = [czero; 2];
    let mut zeta2_arr = [czero; 2];
    let mut sum_arr = [czero; 2];
    let mut init_arr = [0usize; 2];
    let mut caches: [Option<UnikCache<T>>; 3] = [None, None, None];
    let mut cy = [czero; 2];
    let mut kflag: usize = 2;

    let mut i_exit: usize = n; // 1-based Fortran I at loop exit (set to N by default)
    let mut fn_at_exit = fnu;

    for i in 0..n {
        // J flip-flop (Fortran: J = 3 - J) -> Rust: j = 1 - j
        j = 1 - j;
        let fn_val = fnu + T::from_f64(i as f64);

        // Fresh ZUNIK call with IKFLG=2 (K function), IPMTR=0
        let result = zunik(zr, fn_val, IkFlag::K, SumOption::Full, tol, None);
        phi_arr[j] = result.phi;
        zeta1_arr[j] = result.zeta1;
        zeta2_arr[j] = result.zeta2;
        sum_arr[j] = result.sum;
        init_arr[j] = result.cache.init;
        caches[j] = Some(result.cache);

        // Compute S1 exponent (Fortran lines 5797-5808)
        let s1_exp = if kode == Scaling::Exponential {
            let st = zr + result.zeta2;
            let rast = fn_val / zabs(st);
            result.zeta1 - st.conj() * (rast * rast)
        } else {
            result.zeta1 - result.zeta2
        };

        let mut rs1 = s1_exp.re;

        // -- Overflow/underflow test (Fortran lines 5814-5864) --
        if rs1.abs() > elim {
            // label 60: underflow or overflow
            if rs1 > zero {
                return -1; // label 300: overflow
            }
            if z.re < zero {
                return -1; // label 300
            }
            // Underflow
            kdflg = 1;
            y[i] = czero;
            nz += 1;
            if i > 0 && y[i - 1] == czero {
                // already zeroed
            } else if i > 0 {
                y[i - 1] = czero;
                nz += 1;
            }
            continue;
        }

        if kdflg == 1 {
            kflag = 2;
        }
        if rs1.abs() >= alim {
            // Refine test (Fortran lines 5820-5825)
            let aphi = zabs(result.phi);
            rs1 = rs1 + aphi.ln();
            if rs1.abs() > elim {
                // label 60
                if rs1 > zero {
                    return -1;
                }
                if z.re < zero {
                    return -1;
                }
                kdflg = 1;
                y[i] = czero;
                nz += 1;
                if i > 0 && y[i - 1] != czero {
                    y[i - 1] = czero;
                    nz += 1;
                }
                continue;
            }
            if kdflg == 1 {
                kflag = 1;
            }
            if rs1 >= zero && kdflg == 1 {
                kflag = 3;
            }
        }

        // -- Scale and store (Fortran lines 5831-5848) --
        let s2_raw = result.phi * result.sum;
        let s1_scaled = s1_exp.exp() * cssr[kflag - 1];
        let s2 = s2_raw * s1_scaled;

        if kflag == 1 && zuchk(s2, bry[0], tol) {
            // label 60: underflow
            if rs1 > zero {
                return -1;
            }
            if z.re < zero {
                return -1;
            }
            kdflg = 1;
            y[i] = czero;
            nz += 1;
            if i > 0 && y[i - 1] != czero {
                y[i - 1] = czero;
                nz += 1;
            }
            continue;
        }

        // label 50
        cy[kdflg - 1] = s2;
        y[i] = s2 * csrr[kflag - 1];

        if kdflg == 2 {
            i_exit = i + 1; // Fortran 1-based I
            fn_at_exit = fn_val;
            break;
        }
        kdflg = 2;
    }
    // If loop completed without early exit: i_exit = N (Fortran), fn_at_exit = fnu + N - 1

    if i_exit == n {
        fn_at_exit = fnu + T::from_f64((n - 1) as f64);
    }

    // -- Compute RZ and CK (Fortran lines 5868-5874) --
    let rz = reciprocal_z(zr);
    let mut ck = rz * fn_at_exit;

    let ib = i_exit + 1; // Fortran IB = I + 1 (1-based index of next item)

    // -- Test last member and forward recurrence (Fortran lines 5876-5961) --
    if ib <= n {
        // Test last member for underflow (Fortran lines 5881-5921)
        let fn_last = fnu + T::from_f64((n - 1) as f64); // Fortran line 5881
        let ipard = if mr != 0 {
            SumOption::Full
        } else {
            SumOption::SkipSum
        };
        let result_last = zunik(zr, fn_last, IkFlag::K, ipard, tol, None);
        caches[2] = Some(result_last.cache);

        let s1_last = if kode == Scaling::Exponential {
            let st = zr + result_last.zeta2;
            let rast = fn_last / zabs(st);
            result_last.zeta1 - st.conj() * (rast * rast)
        } else {
            result_last.zeta1 - result_last.zeta2
        };

        let mut rs1 = s1_last.re;
        if rs1.abs() > elim {
            // label 95
            if rs1.abs() > zero && rs1 > zero {
                return -1;
            }
            if z.re < zero {
                return -1;
            }
            nz = n as i32;
            y.fill(czero);
            return nz;
        }
        if rs1.abs() >= alim {
            let aphi = zabs(result_last.phi);
            rs1 = rs1 + aphi.ln();
            if rs1.abs() >= elim {
                // label 95
                if rs1 > zero {
                    return -1;
                }
                if z.re < zero {
                    return -1;
                }
                nz = n as i32;
                y.fill(czero);
                return nz;
            }
        }

        // label 100: forward recurrence (Fortran lines 5925-5961)
        let mut s1 = cy[0];
        let mut s2 = cy[1];
        let mut c1r = csrr[kflag - 1];
        let mut ascle = bry[kflag - 1];

        for y_item in y.iter_mut().take(n).skip(ib - 1) {
            let prev = s2;
            s2 = mul_add(ck, prev, s1);
            s1 = prev;
            ck = ck + rz;
            let c2_scaled = s2 * c1r;
            *y_item = c2_scaled;

            if kflag >= 3 {
                continue;
            }
            let c2m = c2_scaled.re.abs().max(c2_scaled.im.abs());
            if c2m <= ascle {
                continue;
            }
            kflag += 1;
            ascle = bry[kflag - 1];
            s1 = s1 * c1r;
            s2 = c2_scaled;
            s1 = s1 * cssr[kflag - 1];
            s2 = s2 * cssr[kflag - 1];
            c1r = csrr[kflag - 1];
        }
    }

    // -- Phase 2: Analytic continuation (Fortran lines 5963-6153) --
    if mr == 0 {
        return nz;
    }

    nz = 0;
    let fmr = T::from_f64(mr as f64);
    let sgn = if fmr < zero { pi_t } else { -pi_t };
    let csgni = sgn;

    // Safety: fnu is finite and < ~1e15 per upper-interface checks
    let inu = fnu.to_i32().unwrap() as usize;
    let fnf = fnu - T::from_f64(inu as f64);
    let ifn = inu + n - 1;
    let ang = fnf * sgn;
    let mut cspn = Complex::new(ang.cos(), ang.sin());
    if ifn % 2 != 0 {
        cspn = -cspn;
    }

    let asc = bry[0];
    let mut iuf: i32 = 0;
    let mut kk = n; // 1-based index counting down
    kdflg = 1;
    let ib_ac = ib.saturating_sub(1); // Fortran: IB = IB - 1 (1-based)
    let ic = ib_ac.saturating_sub(1); // Fortran: IC = IB - 1

    let mut iflag: usize = 2;

    // Phase 2 loop: K=1..N (Fortran labels 172-270)
    let mut k_exit = n; // Fortran K at loop exit
    for k in 0..n {
        let fn_val = fnu + T::from_f64((kk - 1) as f64);

        // -- Cache reuse logic (Fortran labels 172/175/180) --
        let mut m: usize = 2; // workspace slot (0-based), default = slot 2 (Fortran M=3)
        let mut use_cache: Option<UnikCache<T>> = None;

        if n <= 2 {
            // label 172: always try to reuse
            // Copy from J slot, then flip J
            use_cache = caches[j];
            m = j;
            j = 1 - j;
        } else {
            // label 175
            if kk == n && ib_ac < n {
                // First iteration with IB < N: fresh computation (M=3)
                // (INITD=0 implicit because use_cache=None)
            } else if kk == ib_ac || kk == ic {
                // label 172: reuse cache from phase 1
                use_cache = caches[j];
                m = j;
                j = 1 - j;
            }
            // else: INITD=0, fresh computation with M=2 (slot 3)
        }

        // label 180: call ZUNIK with IKFLG=1 (I function)
        let result_i = zunik(zr, fn_val, IkFlag::I, SumOption::Full, tol, use_cache);
        let phid = result_i.phi;
        let zet1d = result_i.zeta1;
        let zet2d = result_i.zeta2;
        let sumd = result_i.sum;
        caches[m] = Some(result_i.cache);

        // Compute S1 exponent for I function (Fortran lines 6019-6030)
        let s1_exp = if kode == Scaling::Exponential {
            let st = zr + zet2d;
            let rast = fn_val / zabs(st);
            st.conj() * (rast * rast) - zet1d
        } else {
            zet2d - zet1d
        };

        let mut rs1 = s1_exp.re;

        // -- Overflow/underflow test (Fortran lines 6035-6096) --
        if rs1.abs() > elim {
            // label 260
            if rs1 > zero {
                return -1; // label 300
            }
            // Underflow: S2 = 0 (goto 230)
            let s2_i = czero;
            let c2_save = s2_i;

            // Fall through to label 230
            cy[kdflg - 1] = s2_i;
            let s2_scaled = s2_i; // already zero

            // Add I and K functions (Fortran lines 6074-6081)
            let mut s1_k = y[kk - 1]; // Y(KK), 1-based
            let mut s2_k = s2_scaled;
            if kode == Scaling::Exponential {
                let s1s2_result = zs1s2(zr, s1_k, s2_k, asc, alim, iuf);
                s1_k = s1s2_result.s1;
                s2_k = s1s2_result.s2;
                nz += s1s2_result.nz;
                iuf = s1s2_result.iuf;
            }
            y[kk - 1] = mul_add(cspn, s1_k, s2_k);
            kk -= 1;
            cspn = -cspn;

            if c2_save == czero {
                kdflg = 1;
            } else if kdflg == 2 {
                k_exit = k + 1;
                break;
            } else {
                kdflg = 2;
            }
            continue;
        }

        if kdflg == 1 {
            iflag = 2;
        }
        if rs1.abs() >= alim {
            // Refine (Fortran lines 6042-6047)
            let aphi = zabs(phid);
            rs1 = rs1 + aphi.ln();
            if rs1.abs() > elim {
                // label 260
                if rs1 > zero {
                    return -1;
                }
                // S2 = 0 path
                let c2_save = czero;
                cy[kdflg - 1] = czero;
                let s2_scaled = czero;

                let mut s1_k = y[kk - 1];
                let mut s2_k = s2_scaled;
                if kode == Scaling::Exponential {
                    let s1s2_result = zs1s2(zr, s1_k, s2_k, asc, alim, iuf);
                    s1_k = s1s2_result.s1;
                    s2_k = s1s2_result.s2;
                    nz += s1s2_result.nz;
                    iuf = s1s2_result.iuf;
                }
                y[kk - 1] = mul_add(cspn, s1_k, s2_k);
                kk -= 1;
                cspn = -cspn;

                if c2_save == czero {
                    kdflg = 1;
                } else if kdflg == 2 {
                    k_exit = k + 1;
                    break;
                } else {
                    kdflg = 2;
                }
                continue;
            }
            if kdflg == 1 {
                iflag = 1;
            }
            if rs1 >= zero && kdflg == 1 {
                iflag = 3;
            }
        }

        // -- Compute I function contribution (Fortran lines 6048-6070) --
        // S2 = -CSGNI * Im(PHI*SUM) + i * CSGNI * Re(PHI*SUM)
        let ps = phid * sumd;
        let mut s2_i = Complex::new(zero, csgni) * ps;

        let s1_sc = s1_exp.exp() * cssr[iflag - 1];
        s2_i = s2_i * s1_sc;

        if iflag == 1 && zuchk(s2_i, bry[0], tol) {
            // NW != 0: set S2 to zero (Fortran lines 6062-6063)
            s2_i = czero;
        }

        // label 230
        cy[kdflg - 1] = s2_i;
        let c2_save = s2_i;
        let s2_scaled = s2_i * csrr[iflag - 1];

        // -- Add I and K functions (Fortran lines 6074-6091) --
        let mut s1_k = y[kk - 1]; // Y(KK)
        let mut s2_k = s2_scaled;

        if kode == Scaling::Exponential {
            let s1s2_result = zs1s2(zr, s1_k, s2_k, asc, alim, iuf);
            s1_k = s1s2_result.s1;
            s2_k = s1s2_result.s2;
            nz += s1s2_result.nz;
            iuf = s1s2_result.iuf;
        }

        // Y(KK) = CSPN * K_part + I_part (Fortran lines 6080-6081)
        y[kk - 1] = mul_add(cspn, s1_k, s2_k);
        kk -= 1;
        cspn = -cspn;

        // KDFLG state machine (Fortran lines 6085-6091)
        if c2_save == czero {
            kdflg = 1;
        } else if kdflg == 2 {
            k_exit = k + 1;
            break;
        } else {
            kdflg = 2;
        }
    }

    // -- Backward recurrence for remaining I terms (Fortran lines 6097-6153) --
    // Fortran: K = N after loop, or K = k_exit if early exit
    // IL = N - K (number of remaining items)
    let il = n - k_exit;
    if il == 0 {
        return nz;
    }

    // Recurrence with scaling (Fortran lines 6107-6153)
    let mut s1 = cy[0];
    let mut s2 = cy[1];
    let mut csr = csrr[iflag - 1];
    let mut ascle = bry[iflag - 1];
    let mut fn_rec = T::from_f64((inu + il) as f64);

    for _i in 0..il {
        let prev = s2;
        let cfn = fn_rec + fnf;
        s2 = s1 + rz * prev * cfn;
        s1 = prev;
        fn_rec = fn_rec - one;
        let ck_val = s2 * csr;

        // Add K function (Fortran lines 6126-6133)
        let mut c1_k = y[kk - 1]; // Y(KK)
        let mut c2_k = ck_val;

        if kode == Scaling::Exponential {
            let s1s2_result = zs1s2(zr, c1_k, c2_k, asc, alim, iuf);
            c1_k = s1s2_result.s1;
            c2_k = s1s2_result.s2;
            nz += s1s2_result.nz;
            iuf = s1s2_result.iuf;
        }

        y[kk - 1] = mul_add(cspn, c1_k, c2_k);
        kk -= 1;
        cspn = -cspn;

        if iflag >= 3 {
            continue;
        }
        let c2m = ck_val.re.abs().max(ck_val.im.abs());
        if c2m <= ascle {
            continue;
        }
        iflag += 1;
        ascle = bry[iflag - 1];
        s1 = s1 * csr;
        s2 = ck_val;
        s1 = s1 * cssr[iflag - 1];
        s2 = s2 * cssr[iflag - 1];
        csr = csrr[iflag - 1];
    }

    nz
}
