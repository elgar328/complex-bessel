// Integration tests for Phase 6 (Airy functions Ai, Bi and their derivatives).
//
// Tests zairy and zbiry against Fortran TOMS 644 reference values
// stored in tests/reference_values/airy_f64.json.
//
// Also tests the Wronskian identity: Ai(z)*Bi'(z) - Ai'(z)*Bi(z) = 1/pi.
//
// Tolerance: 2e-14 (same algorithm, different operation ordering
// between Fortran and Rust).

use complex_bessel::types::AiryDerivative;
use complex_bessel::{airy, airy_scaled, biry, biry_scaled};
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
    #[serde(default)]
    id: Option<i32>,
    kode: i32,
    #[allow(dead_code)]
    #[serde(default)]
    fnu: Option<f64>,
    #[allow(dead_code)]
    #[serde(default)]
    n: Option<i32>,
    #[allow(dead_code)]
    #[serde(default)]
    m: Option<i32>,
}

#[derive(Deserialize)]
#[allow(dead_code)]
struct RefOutputs {
    nz: Option<i32>,
    ierr: Option<i32>,
    cy_re: Option<Vec<f64>>,
    cy_im: Option<Vec<f64>>,
}

const FORTRAN_TOL: f64 = 3e-14;

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

// ── Test: zairy vs Fortran reference values ──

#[test]
fn airy_vs_fortran() {
    let json = include_str!("reference_values/airy_f64.json");
    let refs: RefFile = serde_json::from_str(json).expect("failed to parse airy_f64.json");

    let zairy_tests: Vec<&RefTest> = refs
        .tests
        .iter()
        .filter(|t| t.function == "zairy" && t.outputs.ierr == Some(0))
        .collect();

    assert!(
        !zairy_tests.is_empty(),
        "no zairy tests found in airy_f64.json"
    );

    for t in &zairy_tests {
        let z = Complex::new(t.inputs.z_re, t.inputs.z_im);
        let id_val = t.inputs.id.unwrap_or(0);
        let deriv = match id_val {
            0 => AiryDerivative::Value,
            1 => AiryDerivative::Derivative,
            _ => panic!("{}: unknown id {}", t.label, id_val),
        };
        let kode = t.inputs.kode;

        let result = match kode {
            1 => airy(z, deriv).unwrap_or_else(|e| panic!("{}: unexpected error {:?}", t.label, e)),
            2 => airy_scaled(z, deriv)
                .unwrap_or_else(|e| panic!("{}: unexpected error {:?}", t.label, e)),
            _ => panic!("{}: unknown kode {}", t.label, kode),
        };

        let re_ref = t.outputs.cy_re.as_ref().unwrap()[0];
        let im_ref = t.outputs.cy_im.as_ref().unwrap()[0];

        assert_complex_rel(&t.label, result, re_ref, im_ref, FORTRAN_TOL);
    }

    println!(
        "zairy: all {} test cases passed (tol={:.0e})",
        zairy_tests.len(),
        FORTRAN_TOL
    );
}

// ── Test: zbiry vs Fortran reference values ──

#[test]
fn biry_vs_fortran() {
    let json = include_str!("reference_values/airy_f64.json");
    let refs: RefFile = serde_json::from_str(json).expect("failed to parse airy_f64.json");

    let zbiry_tests: Vec<&RefTest> = refs
        .tests
        .iter()
        .filter(|t| t.function == "zbiry" && t.outputs.ierr == Some(0))
        .collect();

    assert!(
        !zbiry_tests.is_empty(),
        "no zbiry tests found in airy_f64.json"
    );

    for t in &zbiry_tests {
        let z = Complex::new(t.inputs.z_re, t.inputs.z_im);
        let id_val = t.inputs.id.unwrap_or(0);
        let deriv = match id_val {
            0 => AiryDerivative::Value,
            1 => AiryDerivative::Derivative,
            _ => panic!("{}: unknown id {}", t.label, id_val),
        };
        let kode = t.inputs.kode;

        let result = match kode {
            1 => biry(z, deriv).unwrap_or_else(|e| panic!("{}: unexpected error {:?}", t.label, e)),
            2 => biry_scaled(z, deriv)
                .unwrap_or_else(|e| panic!("{}: unexpected error {:?}", t.label, e)),
            _ => panic!("{}: unknown kode {}", t.label, kode),
        };

        let re_ref = t.outputs.cy_re.as_ref().unwrap()[0];
        let im_ref = t.outputs.cy_im.as_ref().unwrap()[0];

        assert_complex_rel(&t.label, result, re_ref, im_ref, FORTRAN_TOL);
    }

    println!(
        "zbiry: all {} test cases passed (tol={:.0e})",
        zbiry_tests.len(),
        FORTRAN_TOL
    );
}

// ── Test: Wronskian identity Ai(z)*Bi'(z) - Ai'(z)*Bi(z) = 1/pi ──

#[test]
fn airy_wronskian_identity() {
    let inv_pi = 1.0 / core::f64::consts::PI;
    let tol = 1e-13;

    let test_points = [
        Complex::new(0.0, 0.0),
        Complex::new(0.5, 0.0),
        Complex::new(2.0, 0.0),
        Complex::new(-1.0, 0.0),
        Complex::new(-3.0, 0.0),
        Complex::new(1.0, 1.0),
        Complex::new(-1.0, 2.0),
        Complex::new(0.5, -0.5),
        Complex::new(3.0, 2.0),
        Complex::new(-2.0, -1.0),
    ];

    for z in &test_points {
        let ai = airy(*z, AiryDerivative::Value).unwrap();
        let ai_p = airy(*z, AiryDerivative::Derivative).unwrap();
        let bi = biry(*z, AiryDerivative::Value).unwrap();
        let bi_p = biry(*z, AiryDerivative::Derivative).unwrap();

        // W = Ai(z)*Bi'(z) - Ai'(z)*Bi(z) should equal 1/pi
        let w = ai * bi_p - ai_p * bi;

        assert!(
            (w.re - inv_pi).abs() < tol,
            "z=({:.1},{:.1}i): Wronskian real part = {:.17e}, expected {:.17e}",
            z.re,
            z.im,
            w.re,
            inv_pi,
        );
        assert!(
            w.im.abs() < tol,
            "z=({:.1},{:.1}i): Wronskian imag part = {:.2e}, expected 0",
            z.re,
            z.im,
            w.im,
        );
    }

    println!(
        "Wronskian identity verified at {} points (tol={:.0e})",
        test_points.len(),
        tol
    );
}
