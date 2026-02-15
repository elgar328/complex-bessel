#!/usr/bin/env python3
"""Generate high-precision reference values using mpmath.

These values serve as cross-verification against the Fortran TOMS 644 references.
mpmath uses a completely different algorithm with arbitrary precision,
so agreement confirms correctness rather than just matching implementations.

Usage:
    python3 scripts/generate_references.py

Requirements:
    pip install mpmath
"""

import json
import sys

try:
    import mpmath
except ImportError:
    print("Error: mpmath is required. Install with: pip install mpmath", file=sys.stderr)
    sys.exit(1)


def compute_besselk(nu, z_re, z_im, digits=50):
    """Compute K_nu(z) using mpmath with high precision."""
    mpmath.mp.dps = digits
    z = mpmath.mpc(z_re, z_im)
    nu_mp = mpmath.mpf(nu)
    result = mpmath.besselk(nu_mp, z)
    return float(result.real), float(result.imag)


def generate_besk_references():
    """Generate K Bessel function reference values."""
    tests = []

    # Same test points as the Fortran generator for cross-verification
    test_cases = [
        # (label, z_re, z_im, fnu, kode, n)
        # Real axis
        ("K_0(1.0)", 1.0, 0.0, 0.0, 1, 1),
        ("K_1(1.0)", 1.0, 0.0, 1.0, 1, 1),
        ("K_0(0.01)", 0.01, 0.0, 0.0, 1, 1),
        ("K_0(0.1)", 0.1, 0.0, 0.0, 1, 1),
        ("K_0(2.0)", 2.0, 0.0, 0.0, 1, 1),
        ("K_1(2.0)", 2.0, 0.0, 1.0, 1, 1),
        ("K_0(5.0)", 5.0, 0.0, 0.0, 1, 1),
        ("K_0(10.0)", 10.0, 0.0, 0.0, 1, 1),
        ("K_0(50.0)", 50.0, 0.0, 0.0, 1, 1),
        ("K_0.5(1.0)", 1.0, 0.0, 0.5, 1, 1),
        ("K_0.25(1.0)", 1.0, 0.0, 0.25, 1, 1),
        ("K_2.5(3.0)", 3.0, 0.0, 2.5, 1, 1),
        ("K_5(2.0)", 2.0, 0.0, 5.0, 1, 1),
        # Complex argument
        ("K_0(1+i)", 1.0, 1.0, 0.0, 1, 1),
        ("K_1(1+i)", 1.0, 1.0, 1.0, 1, 1),
        ("K_0(1.5+0.5i)", 1.5, 0.5, 0.0, 1, 1),
        ("K_0.25(1.5+0.5i)", 1.5, 0.5, 0.25, 1, 1),
        ("K_0(5+3i)", 5.0, 3.0, 0.0, 1, 1),
        ("K_0(0.5+2i)", 0.5, 2.0, 0.0, 1, 1),
        ("K_2(3+4i)", 3.0, 4.0, 2.0, 1, 1),
    ]

    for label, z_re, z_im, fnu, kode, n in test_cases:
        if n == 1:
            re, im = compute_besselk(fnu, z_re, z_im)
            cy_re = [re]
            cy_im = [im]
        else:
            cy_re = []
            cy_im = []
            for j in range(n):
                re, im = compute_besselk(fnu + j, z_re, z_im)
                cy_re.append(re)
                cy_im.append(im)

        test = {
            "function": "besselk",
            "label": label,
            "inputs": {
                "z_re": z_re,
                "z_im": z_im,
                "fnu": fnu,
                "kode": kode,
                "n": n,
            },
            "outputs": {
                "cy_re": cy_re,
                "cy_im": cy_im,
            },
        }
        tests.append(test)

    # Sequence tests
    for label, z_re, z_im, fnu, kode, n in [
        ("K_{0,1,2}(2.0)", 2.0, 0.0, 0.0, 1, 3),
        ("K_{1,2,3}(1+i)", 1.0, 1.0, 1.0, 1, 3),
        ("K_{0.5,1.5,2.5}(3.0)", 3.0, 0.0, 0.5, 1, 3),
    ]:
        cy_re = []
        cy_im = []
        for j in range(n):
            re, im = compute_besselk(fnu + j, z_re, z_im)
            cy_re.append(re)
            cy_im.append(im)

        test = {
            "function": "besselk",
            "label": label,
            "inputs": {
                "z_re": z_re,
                "z_im": z_im,
                "fnu": fnu,
                "kode": kode,
                "n": n,
            },
            "outputs": {
                "cy_re": cy_re,
                "cy_im": cy_im,
            },
        }
        tests.append(test)

    return {
        "generator": f"mpmath {mpmath.__version__} (dps=50)",
        "machine": "Python float64 output from arbitrary precision",
        "tests": tests,
    }


def main():
    data = generate_besk_references()
    json.dump(data, sys.stdout, indent=2)
    sys.stdout.write("\n")
    print(f"Generated {len(data['tests'])} test cases", file=sys.stderr)


if __name__ == "__main__":
    main()
