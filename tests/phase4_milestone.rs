// Integration tests for Phase 4 (I, J, Y Bessel functions + ZACON).
//
// Tests zbesi, zbesj, zbesy against Fortran TOMS 644 reference values
// stored in tests/reference_values/{besi,besj,besy}_f64.json.
//
// Also tests K and H function analytic continuation (ZACON) via
// existing besk/besh reference files with left half-plane arguments.
//
// Tolerance: 2e-14 (same algorithm, different operation ordering
// between Fortran and Rust).

use complex_bessel::types::{HankelKind, Scaling};
use complex_bessel::{besseli_seq, besselj_seq, besselk_seq, bessely_seq, hankel_seq};
use num_complex::Complex;
use serde::Deserialize;

// ── JSON schema ──

#[derive(Deserialize)]
#[allow(dead_code)]
struct RefFile {
    generator: String,
    machine: String,
    tests: Vec<RefTest>,
}

#[derive(Deserialize)]
struct RefTest {
    function: String,
    label: String,
    inputs: RefInputs,
    outputs: RefOutputs,
}

#[derive(Deserialize)]
struct RefInputs {
    z_re: f64,
    z_im: f64,
    fnu: f64,
    kode: i32,
    #[serde(default)]
    m: Option<i32>,
    n: i32,
}

#[derive(Deserialize)]
#[allow(dead_code)]
struct RefOutputs {
    nz: Option<i32>,
    ierr: Option<i32>,
    cy_re: Option<Vec<f64>>,
    cy_im: Option<Vec<f64>>,
    y_re: Option<Vec<f64>>,
    y_im: Option<Vec<f64>>,
}

const FORTRAN_TOL: f64 = 2e-14;

fn assert_complex_rel(label: &str, got: Complex<f64>, re_ref: f64, im_ref: f64, tol: f64) {
    let expected = Complex::new(re_ref, im_ref);
    let scale = expected.norm();
    if scale == 0.0 {
        assert!(
            got.re.abs() < tol && got.im.abs() < tol,
            "{label}: expected ~0, got ({}, {})",
            got.re,
            got.im
        );
        return;
    }
    let err = (got - expected).norm() / scale;
    assert!(
        err < tol,
        "{label}: relative error {err:.2e} exceeds tolerance {tol:.0e}\n  \
         expected: ({re_ref:.17e}, {im_ref:.17e})\n  \
         got:      ({:.17e}, {:.17e})",
        got.re,
        got.im,
    );
}

// ── Test: zbesi ──

#[test]
fn besi_vs_fortran() {
    let json = include_str!("reference_values/besi_f64.json");
    let refs: RefFile = serde_json::from_str(json).expect("failed to parse besi_f64.json");

    let tests: Vec<&RefTest> = refs
        .tests
        .iter()
        .filter(|t| t.function == "zbesi" && t.outputs.ierr == Some(0))
        .collect();

    assert!(!tests.is_empty(), "no zbesi tests found");

    for t in &tests {
        let z = Complex::new(t.inputs.z_re, t.inputs.z_im);
        let fnu = t.inputs.fnu;
        let scaling = match t.inputs.kode {
            1 => Scaling::Unscaled,
            2 => Scaling::Exponential,
            _ => panic!("{}: unknown kode", t.label),
        };
        let n = t.inputs.n as usize;

        let result = besseli_seq(fnu, z, n, scaling)
            .unwrap_or_else(|e| panic!("{}: unexpected error {:?}", t.label, e));

        let cy_re = t.outputs.cy_re.as_ref().expect("missing cy_re");
        let cy_im = t.outputs.cy_im.as_ref().expect("missing cy_im");

        for j in 0..n {
            let sub = if n > 1 {
                format!("{}[{}]", t.label, j)
            } else {
                t.label.clone()
            };
            assert_complex_rel(&sub, result.values[j], cy_re[j], cy_im[j], FORTRAN_TOL);
        }
    }

    eprintln!(
        "  besi_vs_fortran: {} test cases passed (tol={:.0e})",
        tests.len(),
        FORTRAN_TOL
    );
}

// ── Test: zbesj ──

#[test]
fn besj_vs_fortran() {
    let json = include_str!("reference_values/besj_f64.json");
    let refs: RefFile = serde_json::from_str(json).expect("failed to parse besj_f64.json");

    let tests: Vec<&RefTest> = refs
        .tests
        .iter()
        .filter(|t| t.function == "zbesj" && t.outputs.ierr == Some(0))
        .collect();

    assert!(!tests.is_empty(), "no zbesj tests found");

    for t in &tests {
        let z = Complex::new(t.inputs.z_re, t.inputs.z_im);
        let fnu = t.inputs.fnu;
        let scaling = match t.inputs.kode {
            1 => Scaling::Unscaled,
            2 => Scaling::Exponential,
            _ => panic!("{}: unknown kode", t.label),
        };
        let n = t.inputs.n as usize;

        let result = besselj_seq(fnu, z, n, scaling)
            .unwrap_or_else(|e| panic!("{}: unexpected error {:?}", t.label, e));

        let cy_re = t.outputs.cy_re.as_ref().expect("missing cy_re");
        let cy_im = t.outputs.cy_im.as_ref().expect("missing cy_im");

        for j in 0..n {
            let sub = if n > 1 {
                format!("{}[{}]", t.label, j)
            } else {
                t.label.clone()
            };
            assert_complex_rel(&sub, result.values[j], cy_re[j], cy_im[j], FORTRAN_TOL);
        }
    }

    eprintln!(
        "  besj_vs_fortran: {} test cases passed (tol={:.0e})",
        tests.len(),
        FORTRAN_TOL
    );
}

// ── Test: zbesy ──

#[test]
fn besy_vs_fortran() {
    let json = include_str!("reference_values/besy_f64.json");
    let refs: RefFile = serde_json::from_str(json).expect("failed to parse besy_f64.json");

    let tests: Vec<&RefTest> = refs
        .tests
        .iter()
        .filter(|t| t.function == "zbesy" && t.outputs.ierr == Some(0))
        .collect();

    assert!(!tests.is_empty(), "no zbesy tests found");

    for t in &tests {
        let z = Complex::new(t.inputs.z_re, t.inputs.z_im);
        let fnu = t.inputs.fnu;
        let scaling = match t.inputs.kode {
            1 => Scaling::Unscaled,
            2 => Scaling::Exponential,
            _ => panic!("{}: unknown kode", t.label),
        };
        let n = t.inputs.n as usize;

        let result = bessely_seq(fnu, z, n, scaling)
            .unwrap_or_else(|e| panic!("{}: unexpected error {:?}", t.label, e));

        let cy_re = t.outputs.cy_re.as_ref().expect("missing cy_re");
        let cy_im = t.outputs.cy_im.as_ref().expect("missing cy_im");

        for j in 0..n {
            let sub = if n > 1 {
                format!("{}[{}]", t.label, j)
            } else {
                t.label.clone()
            };
            assert_complex_rel(&sub, result.values[j], cy_re[j], cy_im[j], FORTRAN_TOL);
        }
    }

    eprintln!(
        "  besy_vs_fortran: {} test cases passed (tol={:.0e})",
        tests.len(),
        FORTRAN_TOL
    );
}

// ── Test: K function left half-plane (ZACON) ──

#[test]
fn besk_left_half_plane_via_fortran() {
    let json = include_str!("reference_values/besk_f64.json");
    let refs: RefFile = serde_json::from_str(json).expect("failed to parse besk_f64.json");

    // Filter for left half-plane tests (z_re < 0)
    let lhp_tests: Vec<&RefTest> = refs
        .tests
        .iter()
        .filter(|t| t.function == "zbesk" && t.outputs.ierr == Some(0) && t.inputs.z_re < 0.0)
        .collect();

    if lhp_tests.is_empty() {
        eprintln!("  besk_left_half_plane: no left-half-plane tests in reference file (skipped)");
        return;
    }

    for t in &lhp_tests {
        let z = Complex::new(t.inputs.z_re, t.inputs.z_im);
        let fnu = t.inputs.fnu;
        let scaling = match t.inputs.kode {
            1 => Scaling::Unscaled,
            2 => Scaling::Exponential,
            _ => panic!("{}: unknown kode", t.label),
        };
        let n = t.inputs.n as usize;

        let result = besselk_seq(fnu, z, n, scaling)
            .unwrap_or_else(|e| panic!("{}: unexpected error {:?}", t.label, e));

        let cy_re = t.outputs.cy_re.as_ref().expect("missing cy_re");
        let cy_im = t.outputs.cy_im.as_ref().expect("missing cy_im");

        for j in 0..n {
            let sub = if n > 1 {
                format!("{}[{}]", t.label, j)
            } else {
                t.label.clone()
            };
            assert_complex_rel(&sub, result.values[j], cy_re[j], cy_im[j], FORTRAN_TOL);
        }
    }

    eprintln!(
        "  besk_left_half_plane: {} test cases passed (tol={:.0e})",
        lhp_tests.len(),
        FORTRAN_TOL
    );
}

// ── Test: H function full domain (ZACON) ──

#[test]
fn besh_full_domain_via_fortran() {
    let json = include_str!("reference_values/besh_f64.json");
    let refs: RefFile = serde_json::from_str(json).expect("failed to parse besh_f64.json");

    let tests: Vec<&RefTest> = refs
        .tests
        .iter()
        .filter(|t| t.function == "zbesh" && t.outputs.ierr == Some(0))
        .collect();

    assert!(!tests.is_empty(), "no zbesh tests found");

    for t in &tests {
        let z = Complex::new(t.inputs.z_re, t.inputs.z_im);
        let fnu = t.inputs.fnu;
        let m = t.inputs.m.unwrap_or(1);
        let kind = if m == 1 {
            HankelKind::First
        } else {
            HankelKind::Second
        };
        let scaling = match t.inputs.kode {
            1 => Scaling::Unscaled,
            2 => Scaling::Exponential,
            _ => panic!("{}: unknown kode", t.label),
        };
        let n = t.inputs.n as usize;

        let result = hankel_seq(kind, fnu, z, n, scaling)
            .unwrap_or_else(|e| panic!("{}: unexpected error {:?}", t.label, e));

        let cy_re = t.outputs.cy_re.as_ref().expect("missing cy_re");
        let cy_im = t.outputs.cy_im.as_ref().expect("missing cy_im");

        for j in 0..n {
            let sub = if n > 1 {
                format!("{}[{}]", t.label, j)
            } else {
                t.label.clone()
            };
            assert_complex_rel(&sub, result.values[j], cy_re[j], cy_im[j], FORTRAN_TOL);
        }
    }

    eprintln!(
        "  besh_full_domain: {} test cases passed (tol={:.0e})",
        tests.len(),
        FORTRAN_TOL
    );
}

// ── Test: Wronskian identities ──

#[test]
fn wronskian_ik_identity() {
    // I(v,z)*K(v+1,z) + I(v+1,z)*K(v,z) = 1/z
    let test_points = [
        (0.0, Complex::new(1.0, 0.0)),
        (0.5, Complex::new(2.0, 1.0)),
        (1.0, Complex::new(0.5, 0.5)),
        (2.5, Complex::new(3.0, 2.0)),
    ];

    for (fnu, z) in &test_points {
        let i_result = besseli_seq(*fnu, *z, 2, Scaling::Unscaled).unwrap();
        let k_result = besselk_seq(*fnu, *z, 2, Scaling::Unscaled).unwrap();

        let wronskian =
            i_result.values[0] * k_result.values[1] + i_result.values[1] * k_result.values[0];
        let expected = Complex::new(1.0, 0.0) / z;
        let err = (wronskian - expected).norm() / expected.norm();
        assert!(
            err < 1e-13,
            "I*K Wronskian error at fnu={fnu}, z={z}: {err:.2e}"
        );
    }

    eprintln!("  wronskian_ik: {} test points passed", test_points.len());
}

#[test]
fn wronskian_jy_identity() {
    // J(v,z)*Y(v+1,z) - J(v+1,z)*Y(v,z) = -2/(pi*z)
    let pi = std::f64::consts::PI;

    let test_points = [
        (0.0, Complex::new(1.0, 0.0)),
        (0.5, Complex::new(2.0, 0.0)),
        (1.0, Complex::new(3.0, 1.0)),
        (2.5, Complex::new(5.0, 2.0)),
    ];

    for (fnu, z) in &test_points {
        let j_result = besselj_seq(*fnu, *z, 2, Scaling::Unscaled).unwrap();
        let y_result = bessely_seq(*fnu, *z, 2, Scaling::Unscaled).unwrap();

        let wronskian =
            j_result.values[0] * y_result.values[1] - j_result.values[1] * y_result.values[0];
        let expected = Complex::new(-2.0 / pi, 0.0) / z;
        let err = (wronskian - expected).norm() / expected.norm();
        assert!(
            err < 1e-13,
            "J*Y Wronskian error at fnu={fnu}, z={z}: {err:.2e}"
        );
    }

    eprintln!("  wronskian_jy: {} test points passed", test_points.len());
}
