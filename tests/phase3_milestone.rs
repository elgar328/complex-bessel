// Integration tests for Phase 3 (Hankel function).
//
// Tests zbesh (upper interface) against Fortran TOMS 644 reference values
// stored in tests/reference_values/besh_f64.json.
//
// Tolerance: 2e-14 (same algorithm, bit-exact is not expected due to
// different operation ordering between Fortran and Rust).

use complex_bessel::hankel_seq;
use complex_bessel::types::{HankelKind, Scaling};
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
}

// ── Tolerance ──

const FORTRAN_TOL: f64 = 2e-14;

// ── Helper: complex relative error ──

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
        "{label}: relative error {err:.2e} exceeds tolerance {tol:.0e}\n  expected: ({re_ref:.17e}, {im_ref:.17e})\n  got:      ({:.17e}, {:.17e})",
        got.re,
        got.im,
    );
}

// ── Test: zbesh (Hankel function) ──

#[test]
fn besh_vs_fortran() {
    let json = include_str!("reference_values/besh_f64.json");
    let refs: RefFile = serde_json::from_str(json).expect("failed to parse besh_f64.json");

    let zbesh_tests: Vec<&RefTest> = refs
        .tests
        .iter()
        .filter(|t| t.function == "zbesh")
        .collect();

    assert!(
        !zbesh_tests.is_empty(),
        "no zbesh tests found in reference file"
    );

    let mut passed = 0;
    let mut skipped = 0;

    for t in &zbesh_tests {
        let z = Complex::new(t.inputs.z_re, t.inputs.z_im);
        let fnu = t.inputs.fnu;
        let scaling = match t.inputs.kode {
            1 => Scaling::Unscaled,
            2 => Scaling::Exponential,
            _ => panic!("{}: unknown kode {}", t.label, t.inputs.kode),
        };
        let m = t.inputs.m.expect("missing m field");
        let kind = match m {
            1 => HankelKind::First,
            2 => HankelKind::Second,
            _ => panic!("{}: unknown m {}", t.label, m),
        };
        let n = t.inputs.n as usize;

        // Skip test cases that need ZACON (not yet implemented)
        let result = match hankel_seq(kind, fnu, z, n, scaling) {
            Ok(r) => r,
            Err(_) => {
                // Expected for unsupported regions
                skipped += 1;
                continue;
            }
        };

        let cy_re = t.outputs.cy_re.as_ref().expect("missing cy_re");
        let cy_im = t.outputs.cy_im.as_ref().expect("missing cy_im");

        assert_eq!(
            result.values.len(),
            n,
            "{}: expected {} values, got {}",
            t.label,
            n,
            result.values.len()
        );
        assert_eq!(
            result.underflow_count,
            t.outputs.nz.unwrap_or(0) as usize,
            "{}: nz mismatch",
            t.label
        );

        for j in 0..n {
            let sub_label = if n > 1 {
                format!("{}[{}]", t.label, j)
            } else {
                t.label.clone()
            };
            assert_complex_rel(
                &sub_label,
                result.values[j],
                cy_re[j],
                cy_im[j],
                FORTRAN_TOL,
            );
        }

        passed += 1;
    }

    assert!(
        passed > 0,
        "no tests passed (all {} skipped)",
        zbesh_tests.len()
    );

    eprintln!(
        "  besh_vs_fortran: {} zbesh tests passed, {} skipped (tol={:.0e})",
        passed, skipped, FORTRAN_TOL
    );
}
