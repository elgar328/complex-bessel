use complex_bessel::*;
use num_complex::Complex;

fn main() {
    let z = Complex::new(1.0_f64, 2.0);

    // -- Single-value Bessel functions --
    println!("=== Single-value functions (f64) ===");
    let j = besselj(0.5, z).unwrap();
    println!("J_0.5({z}) = {j}");

    let y = bessely(1.0, z).unwrap();
    println!("Y_1({z}) = {y}");

    let i = besseli(0.0, z).unwrap();
    println!("I_0({z}) = {i}");

    let k = besselk(1.0, z).unwrap();
    println!("K_1({z}) = {k}");

    let h1 = hankel1(0.0, z).unwrap();
    println!("H^(1)_0({z}) = {h1}");

    let h2 = hankel2(0.0, z).unwrap();
    println!("H^(2)_0({z}) = {h2}");

    // -- Negative order --
    println!("\n=== Negative order ===");
    let j_neg = besselj(-0.5, z).unwrap();
    println!("J_-0.5({z}) = {j_neg}");

    let k_neg = besselk(-3.0, z).unwrap();
    let k_pos = besselk(3.0, z).unwrap();
    println!("K_-3({z}) = {k_neg}");
    println!("K_3({z})  = {k_pos}  (should be equal)");

    // -- Scaled computation --
    println!("\n=== Scaled functions ===");
    let k_sc = besselk_scaled(1.0, z).unwrap();
    println!("exp(z)*K_1({z}) = {k_sc}");

    let j_sc = besselj_scaled(0.5, z).unwrap();
    println!("exp(-|Im(z)|)*J_0.5({z}) = {j_sc}");

    // -- Sequence computation --
    println!("\n=== Sequence: K_0, K_1, K_2 ===");
    let seq = besselk_seq(0.0, z, 3, Scaling::Unscaled).unwrap();
    for (j, val) in seq.values.iter().enumerate() {
        println!("  K_{j}({z}) = {val}");
    }

    println!("  underflow_count: {}", seq.underflow_count);
    println!("  status: {:?}", seq.status);

    // -- Airy functions --
    println!("\n=== Airy functions ===");
    let ai = airy(z).unwrap();
    println!("Ai({z}) = {ai}");

    let ai_prime = airyprime(z).unwrap();
    println!("Ai'({z}) = {ai_prime}");

    let bi = biry(z).unwrap();
    println!("Bi({z}) = {bi}");

    // -- f32 support --
    println!("\n=== f32 support ===");
    let z32 = Complex::new(1.0_f32, 2.0);
    let j32 = besselj(0.5, z32).unwrap();
    println!("J_0.5({z32}) = {j32} (f32)");
}
