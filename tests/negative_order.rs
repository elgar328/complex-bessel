//! Tests for negative order support and scaled convenience functions.

use complex_bessel::*;
use num_complex::Complex;

type C64 = Complex<f64>;

const TOL: f64 = 5e-14;

fn rel_err(a: C64, b: C64) -> f64 {
    let diff = (a - b).norm();
    let denom = b.norm();
    if denom == 0.0 { diff } else { diff / denom }
}

// ── K_{-ν} = K_ν (K is even in ν) ──

#[test]
fn k_negative_integer_order() {
    let z = C64::new(2.0, 1.0);
    let k_neg = besselk(-3.0, z).unwrap();
    let k_pos = besselk(3.0, z).unwrap();
    assert!(
        rel_err(k_neg, k_pos) < TOL,
        "K_-3 != K_3: {k_neg} vs {k_pos}"
    );
}

#[test]
fn k_negative_fractional_order() {
    let z = C64::new(1.5, 0.5);
    let k_neg = besselk(-0.25, z).unwrap();
    let k_pos = besselk(0.25, z).unwrap();
    assert!(
        rel_err(k_neg, k_pos) < TOL,
        "K_-0.25 != K_0.25: {k_neg} vs {k_pos}"
    );
}

// ── J_{-n} = (-1)^n J_n (integer order) ──

#[test]
fn j_negative_integer_even() {
    let z = C64::new(1.0, 1.0);
    let j_neg = besselj(-2.0, z).unwrap();
    let j_pos = besselj(2.0, z).unwrap();
    // (-1)^2 = 1
    assert!(
        rel_err(j_neg, j_pos) < TOL,
        "J_-2 != J_2: {j_neg} vs {j_pos}"
    );
}

#[test]
fn j_negative_integer_odd() {
    let z = C64::new(1.0, 1.0);
    let j_neg = besselj(-3.0, z).unwrap();
    let j_pos = besselj(3.0, z).unwrap();
    // (-1)^3 = -1
    assert!(
        rel_err(j_neg, -j_pos) < TOL,
        "J_-3 != -J_3: {j_neg} vs {j_pos}"
    );
}

// ── Y_{-n} = (-1)^n Y_n (integer order) ──

#[test]
fn y_negative_integer_even() {
    let z = C64::new(2.0, 1.0);
    let y_neg = bessely(-2.0, z).unwrap();
    let y_pos = bessely(2.0, z).unwrap();
    assert!(
        rel_err(y_neg, y_pos) < TOL,
        "Y_-2 != Y_2: {y_neg} vs {y_pos}"
    );
}

#[test]
fn y_negative_integer_odd() {
    let z = C64::new(2.0, 1.0);
    let y_neg = bessely(-3.0, z).unwrap();
    let y_pos = bessely(3.0, z).unwrap();
    assert!(
        rel_err(y_neg, -y_pos) < TOL,
        "Y_-3 != -Y_3: {y_neg} vs {y_pos}"
    );
}

// ── I_{-n} = I_n (integer order) ──

#[test]
fn i_negative_integer_order() {
    let z = C64::new(1.0, 1.0);
    let i_neg = besseli(-2.0, z).unwrap();
    let i_pos = besseli(2.0, z).unwrap();
    assert!(
        rel_err(i_neg, i_pos) < TOL,
        "I_-2 != I_2: {i_neg} vs {i_pos}"
    );
}

// ── Fractional order reflection formulas ──

#[test]
fn j_negative_fractional_reflection() {
    // J_{-ν}(z) = cos(νπ)*J_ν(z) - sin(νπ)*Y_ν(z)
    let z = C64::new(1.0, 2.0);
    let nu = 0.75;
    let j_neg = besselj(-nu, z).unwrap();

    let j_pos = besselj(nu, z).unwrap();
    let y_pos = bessely(nu, z).unwrap();
    let nu_pi = nu * std::f64::consts::PI;
    let expected = j_pos * nu_pi.cos() - y_pos * nu_pi.sin();

    assert!(
        rel_err(j_neg, expected) < TOL,
        "J_-0.75 reflection: {j_neg} vs {expected}"
    );
}

#[test]
fn y_negative_fractional_reflection() {
    // Y_{-ν}(z) = sin(νπ)*J_ν(z) + cos(νπ)*Y_ν(z)
    let z = C64::new(1.0, 2.0);
    let nu = 0.75;
    let y_neg = bessely(-nu, z).unwrap();

    let j_pos = besselj(nu, z).unwrap();
    let y_pos = bessely(nu, z).unwrap();
    let nu_pi = nu * std::f64::consts::PI;
    let expected = j_pos * nu_pi.sin() + y_pos * nu_pi.cos();

    assert!(
        rel_err(y_neg, expected) < TOL,
        "Y_-0.75 reflection: {y_neg} vs {expected}"
    );
}

#[test]
fn i_negative_fractional_reflection() {
    // I_{-ν}(z) = I_ν(z) + (2/π)*sin(νπ)*K_ν(z)
    let z = C64::new(1.0, 2.0);
    let nu = 1.5;
    let i_neg = besseli(-nu, z).unwrap();

    let i_pos = besseli(nu, z).unwrap();
    let k_pos = besselk(nu, z).unwrap();
    let nu_pi = nu * std::f64::consts::PI;
    let expected = i_pos + k_pos * (2.0 / std::f64::consts::PI * nu_pi.sin());

    assert!(
        rel_err(i_neg, expected) < TOL,
        "I_-1.5 reflection: {i_neg} vs {expected}"
    );
}

// ── Hankel negative order ──

#[test]
fn h1_negative_order() {
    // H^(1)_{-ν}(z) = exp(νπi) * H^(1)_ν(z)
    let z = C64::new(2.0, 1.0);
    let nu = 1.25;
    let h_neg = hankel(HankelKind::First, -nu, z).unwrap();

    let h_pos = hankel(HankelKind::First, nu, z).unwrap();
    let nu_pi = nu * std::f64::consts::PI;
    let rotation = C64::new(nu_pi.cos(), nu_pi.sin());
    let expected = h_pos * rotation;

    assert!(
        rel_err(h_neg, expected) < TOL,
        "H1_-1.25 reflection: {h_neg} vs {expected}"
    );
}

#[test]
fn h2_negative_order() {
    // H^(2)_{-ν}(z) = exp(-νπi) * H^(2)_ν(z)
    let z = C64::new(2.0, 1.0);
    let nu = 1.25;
    let h_neg = hankel(HankelKind::Second, -nu, z).unwrap();

    let h_pos = hankel(HankelKind::Second, nu, z).unwrap();
    let nu_pi = nu * std::f64::consts::PI;
    let rotation = C64::new(nu_pi.cos(), -nu_pi.sin());
    let expected = h_pos * rotation;

    assert!(
        rel_err(h_neg, expected) < TOL,
        "H2_-1.25 reflection: {h_neg} vs {expected}"
    );
}

// ── Scaled functions ──

#[test]
fn scaled_besselk_consistency() {
    // exp(z)*K_ν(z) should equal besselk_scaled(ν, z)
    let z = C64::new(1.0, 2.0);
    let k_unscaled = besselk(1.0, z).unwrap();
    let k_scaled = besselk_scaled(1.0, z).unwrap();
    let expected = k_unscaled * z.exp();
    assert!(
        rel_err(k_scaled, expected) < TOL,
        "scaled K: {k_scaled} vs {expected}"
    );
}

#[test]
fn scaled_besseli_consistency() {
    // exp(-|Re(z)|)*I_ν(z) should equal besseli_scaled(ν, z)
    let z = C64::new(1.0, 2.0);
    let i_unscaled = besseli(0.5, z).unwrap();
    let i_scaled = besseli_scaled(0.5, z).unwrap();
    let expected = i_unscaled * (-z.re.abs()).exp();
    assert!(
        rel_err(i_scaled, expected) < TOL,
        "scaled I: {i_scaled} vs {expected}"
    );
}

#[test]
fn scaled_besselj_consistency() {
    // exp(-|Im(z)|)*J_ν(z) should equal besselj_scaled(ν, z)
    let z = C64::new(1.0, 2.0);
    let j_unscaled = besselj(0.5, z).unwrap();
    let j_scaled = besselj_scaled(0.5, z).unwrap();
    let expected = j_unscaled * (-z.im.abs()).exp();
    assert!(
        rel_err(j_scaled, expected) < TOL,
        "scaled J: {j_scaled} vs {expected}"
    );
}

#[test]
fn scaled_bessely_consistency() {
    // exp(-|Im(z)|)*Y_ν(z) should equal bessely_scaled(ν, z)
    let z = C64::new(1.0, 2.0);
    let y_unscaled = bessely(1.0, z).unwrap();
    let y_scaled = bessely_scaled(1.0, z).unwrap();
    let expected = y_unscaled * (-z.im.abs()).exp();
    assert!(
        rel_err(y_scaled, expected) < TOL,
        "scaled Y: {y_scaled} vs {expected}"
    );
}

#[test]
fn scaled_hankel_h1_consistency() {
    // exp(-iz)*H^(1)_ν(z) should equal hankel_scaled(First, ν, z)
    let z = C64::new(2.0, 1.0);
    let h_unscaled = hankel(HankelKind::First, 0.0, z).unwrap();
    let h_scaled = hankel_scaled(HankelKind::First, 0.0, z).unwrap();
    let iz = C64::new(-z.im, z.re); // i*z
    let expected = h_unscaled * (-iz).exp();
    assert!(
        rel_err(h_scaled, expected) < TOL,
        "scaled H1: {h_scaled} vs {expected}"
    );
}

// ── Scaled + negative order ──

#[test]
fn scaled_k_negative_order() {
    let z = C64::new(2.0, 1.0);
    let k_neg = besselk_scaled(-1.5, z).unwrap();
    let k_pos = besselk_scaled(1.5, z).unwrap();
    assert!(
        rel_err(k_neg, k_pos) < TOL,
        "scaled K_-1.5 != K_1.5: {k_neg} vs {k_pos}"
    );
}

#[test]
fn scaled_j_negative_integer() {
    let z = C64::new(1.0, 2.0);
    let j_neg = besselj_scaled(-3.0, z).unwrap();
    let j_pos = besselj_scaled(3.0, z).unwrap();
    assert!(
        rel_err(j_neg, -j_pos) < TOL,
        "scaled J_-3 != -J_3: {j_neg} vs {j_pos}"
    );
}

// ── Wronskian identity with negative order ──

#[test]
fn wronskian_jy_negative_order() {
    // J_{ν+1}(z)*Y_ν(z) - J_ν(z)*Y_{ν+1}(z) = 2/(πz) (DLMF 10.5.3)
    // Use ν = -0.5
    let z = C64::new(2.0, 1.0);
    let nu = -0.5;
    let j0 = besselj(nu, z).unwrap();
    let j1 = besselj(nu + 1.0, z).unwrap();
    let y0 = bessely(nu, z).unwrap();
    let y1 = bessely(nu + 1.0, z).unwrap();

    let wronskian = j1 * y0 - j0 * y1;
    let expected = C64::new(2.0 / std::f64::consts::PI, 0.0) / z;

    assert!(
        rel_err(wronskian, expected) < 1e-13,
        "Wronskian J/Y at nu={nu}: {wronskian} vs {expected}"
    );
}

#[test]
fn wronskian_ik_negative_order() {
    // I_ν(z)*K_{ν+1}(z) + I_{ν+1}(z)*K_ν(z) = 1/z
    // Use ν = -1.5
    let z = C64::new(2.0, 1.0);
    let nu = -1.5;
    let i0 = besseli(nu, z).unwrap();
    let i1 = besseli(nu + 1.0, z).unwrap();
    let k0 = besselk(nu, z).unwrap();
    let k1 = besselk(nu + 1.0, z).unwrap();

    let wronskian = i0 * k1 + i1 * k0;
    let expected = C64::new(1.0, 0.0) / z;

    assert!(
        rel_err(wronskian, expected) < 1e-13,
        "Wronskian I/K at nu={nu}: {wronskian} vs {expected}"
    );
}

// ── Sequence functions reject negative order ──

#[test]
fn seq_rejects_negative_order() {
    let z = C64::new(1.0, 1.0);
    assert!(besselj_seq(-1.0, z, 1, Scaling::Unscaled).is_err());
    assert!(bessely_seq(-1.0, z, 1, Scaling::Unscaled).is_err());
    assert!(besseli_seq(-1.0, z, 1, Scaling::Unscaled).is_err());
    assert!(besselk_seq(-1.0, z, 1, Scaling::Unscaled).is_err());
    assert!(hankel_seq(HankelKind::First, -1.0, z, 1, Scaling::Unscaled).is_err());
}
