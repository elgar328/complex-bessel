//! Componentwise accuracy of values on the real and imaginary axes.
//!
//! On the half-axes z = x, ix, −x, −ix (x > 0), every function of real order
//! ν is a combination α·A + β·B of two real functions of order m = |ν| at
//! the real argument x, either (J_m, Y_m) or (I_m, K_m), with coefficients
//! that depend on ν only (DLMF 10.4.6, 10.11.2/5, 10.27.6/8/11, 10.34.2).
//! When A and B differ greatly in size, one component of the result is
//! tiny. Amos computes the value as one complex number (e.g. H from
//! K(−iz)), so that component carries the rounding error of the other,
//! ε·|large|: accurate normwise, but possibly without a single correct digit
//! (e.g. `hankel1(10, 1).re` = 8.4e-9 instead of J₁₀(1) = 2.6e-10).
//!
//! This module recomputes such a small component from A and B. The
//! coefficient table already contains the reflection formulas for ν < 0:
//! rotating a corrected value of order |ν| afterwards would cancel again
//! (e.g. H1_{−3.3}(−x) = −H2_{3.3}(x)).
//!
//! Only the values below have a small component that Amos gets wrong (apart
//! from orders next to an integer or a half-integer, where a nearly vanishing
//! cos(νπ) or sin(νπ) times one real function keeps a residue relative to
//! |f(z)|; see the crate docs); the others on the axes are a single real
//! function times a phase, or Amos already evaluates the same two-term formula (Y on the imaginary axis, K on
//! the negative real axis). Scaled values whose scaling factor is a phase
//! (H on the real axis, K on the imaginary axis) mix the components, so they
//! are left alone.

use num_complex::Complex;

use crate::machine::BesselFloat;
use crate::types::{HankelKind, Scaling};
use crate::utils::{self, Workspace, cospi, sinpi};
use crate::{besh, besi, besj, besk};

/// The function a value on an axis belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AxisFn {
    Hankel(HankelKind),
    K,
    Y,
}

/// A (function, half-axis) pair whose small component is recomputed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Row {
    /// H^(1)(x), H^(2)(x)
    H1Pos,
    H2Pos,
    /// H^(1)(−x), H^(2)(−x), Y(−x)
    H1Neg,
    H2Neg,
    YNeg,
    /// K(ix), K(−ix)
    KPosI,
    KNegI,
    /// H^(2)(ix), H^(1)(−ix)
    H2PosI,
    H1NegI,
}

impl Row {
    /// Whether the basis is (I, K) rather than (J, Y).
    fn ik(self) -> bool {
        matches!(self, Row::H2PosI | Row::H1NegI)
    }
}

// Recompute the small component when |small| < RATIO·|large|. Measured on
// ~1e5 small components (orders 0.05–300 of both signs, including integers,
// half-integers and orders 1e-12…1e-4 away from them; mpmath references):
// - (J, Y) for x < m: with 0.1 the largest remaining error is 2.3e-11
//   (20 components above 1e-12, worst increase 5.5e-13 → 3.2e-12); with
//   1e-3 it is 2.8e-11 with 174 above 1e-12. Above m, next to the zeros of
//   J_m and Y_m (which all lie beyond m), J and Y computed separately are
//   in general no more accurate than Amos's components.
// - On the negative real axis, Y(−x) and H^(2)(−x) are the exceptions for
//   x ≥ m, in two cases:
//   - Next to a half-integer order the real part of Y(−x) has small
//     coefficients (e.g. Re Y_ν(−x) = cos(νπ)·Y_ν(x)), and Amos's has a
//     relative error of ε/|cos(νπ)| or more. Recomputed when |α| + |β| <
//     COEF_Y_NEG. In 1.2e5 values beyond m (half of the orders next to
//     half-integers, the others next to integers or generic) this fixes
//     errors up to 54 (f32 150); the worst increase is 8.2e-14 → 1.9e-13
//     (f32 1.4e-3 → 8.4e-3), and nothing gets worse normwise.
//   - A component that is Y_m alone, |α| < RATIO_J_NEG·|β|, for
//     m < ORDER_NEG: Re Y_ν(−x) = cos(νπ)·Y_ν(x) for ν ≥ 0, and at or
//     within ~1e-4 of an integer or a half-integer order the others whose
//     J coefficient vanishes there (e.g. Im H^(2)_n(−x) = −(−1)^n·Y_n(x)).
//     Measured in units of ε·(|α| + |β|)·|H_m| (the relative error is
//     meaningless next to a zero). For Re Y(−x), ν ≥ 0, on 1.5e4 values
//     next to zeros, Amos's errors reach 617 (p99) and 1937 units, the
//     recomputed ones 5 and 20, and 5 of 3620 get worse, by at most 15. For
//     the others, on 5.0e4 components next to zeros with x up to 2010 and
//     dense where Amos's separate J errs most (f32 x ≈ 33–46, f64
//     x ≈ 84–102), Amos's reach 126 (p99) and 269 units (f32 1.4e4), the
//     recomputed ones 7 and 21, and 61 get worse by more than 5 units, by
//     at most 13. Above m ≈ 18, Y_m computed separately is less accurate
//     (its error grows with m): ~9% get worse, by up to 80 units. A larger
//     ratio lets J's error in (up to 5.8e4 units for f32 at x ≈ 35–45, as
//     in Fortran's CBESJ; +184 units with ratios up to 0.1). Of the
//     components whose J coefficient is not small, Im Y(−x) gets worse more
//     often than better, and Re Y(−x) for ν < 0 away from integers, though
//     mostly better, gets worse by up to 132 units (f32 1.3e4). H^(1)(−x)
//     is left alone: Amos's components are already at that level, and so
//     are all those on the positive real axis, where J computed separately
//     is in general less accurate than Amos's component.
// - (I, K): re-measured including f32 and a separate sample of orders next
//   to integers (where these small components occur). Raising the ratio
//   from 1e-3 to 1e-2 lowers the largest error left unfixed from 1.6e-4 to
//   2.2e-5 for f32 and from 2.7e-13 to 1.4e-13 for f64. The worst case where
//   the recomputed component is less accurate than Amos's stays the same.
//   Ratios above 1e-2 add recomputations for no measurable gain in f64.
// Next to x = m the recomputed component can be a few times less accurate
// than Amos's (e.g. 6.7e-14 → 7.4e-13 for Im H^(2)_{−10.7}(−10)); the
// ratios above keep such cases rare and small.
const RATIO_JY: f64 = 0.1;
const RATIO_IK: f64 = 1e-2;
const COEF_Y_NEG: f64 = 0.1;
const RATIO_J_NEG: f64 = 1e-3;
const ORDER_NEG: f64 = 18.0;

/// The row of `(f, z, scaling)` and x = |z|, if z is on a half-axis and the
/// value has a small component to recompute. Signed zeros are treated like
/// +0, as everywhere in the crate (the negative real axis is the upper side
/// of the branch cut).
fn row<T: BesselFloat>(f: AxisFn, z: Complex<T>, scaling: Scaling) -> Option<(Row, T)> {
    let zero = T::zero();
    let unscaled = scaling == Scaling::Unscaled;
    let first = AxisFn::Hankel(HankelKind::First);
    let second = AxisFn::Hankel(HankelKind::Second);
    let row = if z.im == zero && z.re > zero {
        match f {
            _ if !unscaled => return None,
            _ if f == first => Row::H1Pos,
            _ if f == second => Row::H2Pos,
            _ => return None,
        }
    } else if z.im == zero && z.re < zero {
        match f {
            // Y's scaling factor is exp(−|Im z|) = 1 here
            AxisFn::Y => Row::YNeg,
            _ if !unscaled => return None,
            _ if f == first => Row::H1Neg,
            _ if f == second => Row::H2Neg,
            _ => return None,
        }
    } else if z.re == zero && z.im > zero {
        match f {
            // H^(2)'s scaling factor is exp(iz) = exp(−x), real
            _ if f == second => Row::H2PosI,
            AxisFn::K if unscaled => Row::KPosI,
            _ => return None,
        }
    } else if z.re == zero && z.im < zero {
        match f {
            _ if f == first => Row::H1NegI,
            AxisFn::K if unscaled => Row::KNegI,
            _ => return None,
        }
    } else {
        return None;
    };
    let x = if z.im == zero { z.re.abs() } else { z.im.abs() };
    Some((row, x))
}

/// Index of the small component of `v` (0 = real, 1 = imaginary) if it is
/// to be recomputed for order ν = ±m (`neg` for ν < 0).
fn small_component<T: BesselFloat>(
    row: Row,
    neg: bool,
    x: T,
    m: T,
    v: Complex<T>,
) -> Option<usize> {
    let (a, b) = (v.re.abs(), v.im.abs());
    if !(a.is_finite() && b.is_finite()) {
        return None;
    }
    let (k, small, large) = if a < b { (0, a, b) } else { (1, b, a) };
    if large == T::zero() {
        return None;
    }
    let ratio = if row.ik() { RATIO_IK } else { RATIO_JY };
    if small >= T::from_f64(ratio) * large {
        return None;
    }
    // (J, Y) beyond m only for Y(−x) and H^(2)(−x): components that are Y_m
    // alone below order ORDER_NEG, or the real part of Y(−x) with small
    // coefficients; see RATIO_JY.
    if !row.ik() && x >= m {
        if !matches!(row, Row::YNeg | Row::H2Neg) {
            return None;
        }
        let (ca, cb) = coefficients(row, neg, m, k);
        let y_alone = m < T::from_f64(ORDER_NEG) && ca.abs() < T::from_f64(RATIO_J_NEG) * cb.abs();
        let small_coef =
            row == Row::YNeg && k == 0 && ca.abs() + cb.abs() < T::from_f64(COEF_Y_NEG);
        if !(y_alone || small_coef) {
            return None;
        }
    }
    Some(k)
}

/// sin(1.5πm) and cos(1.5πm) without the rounding of the product 1.5·m.
fn sincospi_three_halves<T: BesselFloat>(m: T) -> (T, T) {
    // 1.5·4k = 6k is even, so m can be reduced mod 4 (exactly)
    let r = m % T::from_f64(4.0);
    let h = r * T::from_f64(0.5);
    // TwoSum: s + e = r + h exactly
    let s = r + h;
    let bb = s - r;
    let e = (r - (s - bb)) + (h - bb);
    let pe = T::from_f64(core::f64::consts::PI) * e;
    let (sp, cp) = (sinpi(s), cospi(s));
    (sp + pe * cp, cp - pe * sp)
}

/// Coefficients of A and B for component k (0 = real, 1 = imaginary) of the
/// value of order ν = ±m (`neg` for ν < 0). The table was derived with
/// sympy from the DLMF formulas and checked against them at 300 digits.
fn coefficients<T: BesselFloat>(row: Row, neg: bool, m: T, k: usize) -> (T, T) {
    let one = T::one();
    let two = T::from_f64(2.0);
    let pi = T::from_f64(core::f64::consts::PI);
    let hpi = pi / two;
    let tpi = two / pi;
    // (A re, A im, B re, B im)
    let (ar, ai, br, bi) = match row {
        Row::H1Pos | Row::H2Pos => {
            let sign = if row == Row::H1Pos { one } else { -one };
            if neg {
                let (c, s) = (cospi(m), sinpi(m));
                (c, sign * s, -s, sign * c)
            } else {
                (one, T::zero(), T::zero(), sign)
            }
        }
        Row::H1Neg => {
            if neg {
                (-one, T::zero(), T::zero(), one)
            } else {
                let (c, s) = (cospi(m), sinpi(m));
                (-c, s, s, c)
            }
        }
        Row::H2Neg => {
            if neg {
                let (c2, s2) = (cospi(two * m), sinpi(two * m));
                (c2 + two, -s2, -s2, -c2)
            } else {
                let (c, s) = (cospi(m), sinpi(m));
                (T::from_f64(3.0) * c, s, -s, -c)
            }
        }
        Row::YNeg => {
            let c = cospi(m);
            if neg {
                let half_s2 = sinpi(two * m) / two;
                (half_s2, c * c + one, c * c, -half_s2)
            } else {
                (T::zero(), two * c, c, -sinpi(m))
            }
        }
        Row::KPosI | Row::KNegI => {
            // K_{−ν} = K_ν, so the sign of ν does not matter
            let (ch, sh) = (cospi(m / two), sinpi(m / two));
            let conj = if row == Row::KPosI { one } else { -one };
            (-hpi * sh, -conj * hpi * ch, -hpi * ch, conj * hpi * sh)
        }
        Row::H2PosI | Row::H1NegI => {
            // H^(1)(−ix) = conj(H^(2)(ix)) for real ν
            let conj = if row == Row::H2PosI { one } else { -one };
            let (ch, sh) = (cospi(m / two), sinpi(m / two));
            if neg {
                let (s3, c3) = sincospi_three_halves(m);
                (two * ch, -conj * two * sh, tpi * s3, conj * tpi * c3)
            } else {
                (two * ch, conj * two * sh, tpi * sh, conj * tpi * ch)
            }
        }
    };
    if k == 0 { (ar, br) } else { (ai, bi) }
}

/// The basis (A, B) at order m and the real argument x: (J_m, Y_m) or
/// (I_m, K_m), scaled by exp(−x) for the scaled Hankel rows. `None` if an
/// Amos call fails (the original value is kept then).
fn basis<T: BesselFloat>(row: Row, x: T, m: T, scaling: Scaling) -> Option<(T, T)> {
    let zero = T::zero();
    let czero = Complex::new(zero, zero);
    let z = Complex::new(x, zero);
    let mut a = [czero];
    let mut b = [czero];
    if row.ik() {
        besi::zbesi(z, m, scaling, &mut a).ok()?;
        besk::zbesk(z, m, scaling, &mut b).ok()?;
        Some((a[0].re, scale_k(b[0].re, x, scaling)))
    } else {
        besj::zbesj(z, m, Scaling::Unscaled, &mut a).ok()?;
        besh::zbesh(z, m, HankelKind::First, Scaling::Unscaled, &mut b).ok()?;
        Some((a[0].re, b[0].im))
    }
}

/// exp(−x)·K from Amos's scaled exp(x)·K, in two steps (exp(−2x) alone is
/// subnormal above x ≈ 354 for f64).
#[inline]
fn scale_k<T: BesselFloat>(ks: T, x: T, scaling: Scaling) -> T {
    if scaling == Scaling::Exponential {
        let e = (-x).exp();
        ks * e * e
    } else {
        ks
    }
}

/// Set component k of `v` from the basis.
#[inline]
fn apply<T: BesselFloat>(row: Row, neg: bool, m: T, k: usize, a: T, b: T, v: &mut Complex<T>) {
    let (ca, cb) = coefficients(row, neg, m, k);
    let s = ca * a + cb * b;
    if k == 0 {
        v.re = s;
    } else {
        v.im = s;
    }
}

/// Recompute the small component of a single value `v` = f_ν(z) if needed.
pub(crate) fn fix_single<T: BesselFloat>(
    f: AxisFn,
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    v: Complex<T>,
) -> Complex<T> {
    let Some((row, x)) = row(f, z, scaling) else {
        return v;
    };
    let m = nu.abs();
    let neg = nu < T::zero();
    let Some(k) = small_component(row, neg, x, m, v) else {
        return v;
    };
    // For H on the positive real axis with ν ≥ 0, the small component is
    // Re = J: Y is not needed
    let ab = if !neg && matches!(row, Row::H1Pos | Row::H2Pos) {
        let mut a = [Complex::new(T::zero(), T::zero())];
        besj::zbesj(Complex::new(x, T::zero()), m, Scaling::Unscaled, &mut a)
            .ok()
            .map(|_| (a[0].re, T::zero()))
    } else {
        basis(row, x, m, scaling)
    };
    let mut out = v;
    if let Some((a, b)) = ab {
        apply(row, neg, m, k, a, b, &mut out);
    }
    out
}

/// Recompute the small components of a sequence `out[j]` = f_{ν+j}(z).
///
/// The orders below zero and the others each form a lattice of |orders|;
/// for each, the basis is computed as sequences over the range of elements
/// that need it, in chunks that fit a stack buffer (so `_seq` and
/// `_seq_into` agree). A chunk whose Amos call fails falls back to single
/// values (e.g. a chunk of K can overflow where the whole call does not).
pub(crate) fn fix_seq<T: BesselFloat>(
    f: AxisFn,
    nu: T,
    z: Complex<T>,
    scaling: Scaling,
    out: &mut [Complex<T>],
) {
    let Some((row, x)) = row(f, z, scaling) else {
        return;
    };
    let n = out.len();
    let zero = T::zero();
    // Number of negative orders; |ν + j| = |ν| − j for those (exact)
    let nc = if nu < zero {
        nu.abs().ceil().to_usize().unwrap_or(n).min(n)
    } else {
        0
    };
    let abs_nu = nu.abs();
    let order = |j: usize| -> T {
        if j < nc {
            abs_nu - T::from_f64(j as f64)
        } else {
            nu + T::from_f64(j as f64)
        }
    };
    let needs = |j: usize, v: Complex<T>| small_component(row, j < nc, x, order(j), v);
    for (lo, hi, neg) in [(0, nc, true), (nc, n, false)] {
        let Some(first) = (lo..hi).find(|&j| needs(j, out[j]).is_some()) else {
            continue;
        };
        let last = (first..hi)
            .rev()
            .find(|&j| needs(j, out[j]).is_some())
            .unwrap_or(first);
        let len = last - first + 1;
        // Lowest |order| of the range: ascending orders run with j for the
        // non-negative part and against j for the negative one
        let (m_lo, j_of) = if neg {
            (order(last), last)
        } else {
            (order(first), first)
        };
        let j_at = |i: usize| if neg { j_of - i } else { j_of + i };
        // H on the positive real axis with ν ≥ 0: Re = J, Y is not needed
        let j_only = !neg && matches!(row, Row::H1Pos | Row::H2Pos);
        utils::with_workspace(Workspace::Stack, 0, |buf: &mut [Complex<T>]| {
            let (abuf, bbuf) = buf.split_at_mut(buf.len() / 2);
            let chunk = abuf.len();
            let mut c0 = 0;
            while c0 < len {
                let c1 = (c0 + chunk).min(len);
                let cl = c1 - c0;
                // Same ulp grid as the orders, so exact
                let start = m_lo + T::from_f64(c0 as f64);
                let zx = Complex::new(x, zero);
                let ok = if row.ik() {
                    besi::zbesi(zx, start, scaling, &mut abuf[..cl]).is_ok()
                        && besk::zbesk(zx, start, scaling, &mut bbuf[..cl]).is_ok()
                } else {
                    besj::zbesj(zx, start, Scaling::Unscaled, &mut abuf[..cl]).is_ok()
                        && (j_only
                            || besh::zbesh(
                                zx,
                                start,
                                HankelKind::First,
                                Scaling::Unscaled,
                                &mut bbuf[..cl],
                            )
                            .is_ok())
                };
                for i in 0..cl {
                    let j = j_at(c0 + i);
                    let Some(k) = needs(j, out[j]) else {
                        continue;
                    };
                    let m = order(j);
                    if ok {
                        let a = abuf[i].re;
                        let b = if row.ik() {
                            scale_k(bbuf[i].re, x, scaling)
                        } else if j_only {
                            zero
                        } else {
                            bbuf[i].im
                        };
                        apply(row, neg, m, k, a, b, &mut out[j]);
                    } else {
                        let nu_j = if neg { -m } else { m };
                        out[j] = fix_single(f, nu_j, z, scaling, out[j]);
                    }
                }
                c0 = c1;
            }
        });
    }
}

/// J and Y of real order are real on the positive real axis; drop the
/// rounding noise Amos leaves in the imaginary part.
#[inline]
pub(crate) fn real_on_positive_axis<T: BesselFloat>(z: Complex<T>, v: Complex<T>) -> Complex<T> {
    if z.im == T::zero() && z.re > T::zero() {
        Complex::new(v.re, T::zero())
    } else {
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// sin and cos of 1.5πm where 1.5·m is not representable and lies next
    /// to an integer: the naive product rounds onto the integer (sin = 0).
    #[test]
    fn sincospi_three_halves_keeps_the_rounding() {
        // (m, sin(1.5πm), cos(1.5πm)) from mpmath with the exact product
        let cases = [
            (2.6666666666666665f64, -6.975736996017264e-16, 1.0),
            (3.3333333333333335, -6.975736996017264e-16, -1.0),
            (67.33333333333333, 2.2322358387255244e-14, -1.0),
        ];
        for (m, s, c) in cases {
            let (gs, gc) = sincospi_three_halves(m);
            assert!(
                ((gs - s) / s).abs() < 1e-12,
                "sin(1.5π·{m}) = {gs:e}, expected {s:e}"
            );
            assert!(
                (gc - c).abs() < 1e-15,
                "cos(1.5π·{m}) = {gc:e}, expected {c:e}"
            );
        }
        // Representable products are unchanged
        for m in [0.3f64, 2.0, 5.5, 1.999999] {
            let (gs, gc) = sincospi_three_halves(m);
            assert!((gs - sinpi(1.5 * m)).abs() <= 2.0 * f64::EPSILON, "m={m}");
            assert!((gc - cospi(1.5 * m)).abs() <= 2.0 * f64::EPSILON, "m={m}");
        }
    }
}
