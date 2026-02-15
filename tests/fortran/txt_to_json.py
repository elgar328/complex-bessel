#!/usr/bin/env python3
"""Convert structured text output from Fortran reference generators to JSON.

Usage: python3 txt_to_json.py < gen_besk_refs.txt > besk_f64.json
"""

import json
import re
import sys


def parse_fortran_output(text):
    """Parse the structured text output from Fortran generators."""
    tests = []
    current = None
    generator = "Fortran TOMS 644 (zbsubs.f, revision 930101)"
    machine = "IEEE 754 double precision"

    for line in text.splitlines():
        line = line.strip()

        # Skip comments and empty lines
        if not line or line.startswith("#"):
            continue

        # Begin a new test case
        if line.startswith("BEGIN "):
            parts = line[6:].split(None, 1)
            func = parts[0]
            label = parts[1] if len(parts) > 1 else ""
            current = {"function": func, "label": label, "inputs": {}, "outputs": {}}
            continue

        # End current test case
        if line == "END":
            if current:
                tests.append(current)
                current = None
            continue

        if current is None:
            continue

        # Parse key-value pairs
        # Input fields
        if line.startswith("z_re  "):
            current["inputs"]["z_re"] = _parse_float(line[6:])
        elif line.startswith("z_im  "):
            current["inputs"]["z_im"] = _parse_float(line[6:])
        elif line.startswith("fnu   "):
            current["inputs"]["fnu"] = _parse_float(line[6:])
        elif line.startswith("kode  "):
            current["inputs"]["kode"] = int(line[6:].strip())
        elif line.startswith("id    "):
            current["inputs"]["id"] = int(line[6:].strip())
        elif line.startswith("m     "):
            current["inputs"]["m"] = int(line[6:].strip())
        elif line.startswith("n     "):
            current["inputs"]["n"] = int(line[6:].strip())

        # Output fields
        elif line.startswith("nz    "):
            current["outputs"]["nz"] = int(line[6:].strip())
        elif line.startswith("ierr  "):
            current["outputs"]["ierr"] = int(line[6:].strip())

        # Array outputs: cy_re(J), cy_im(J), y_re(J), y_im(J)
        else:
            m = re.match(r"(cy_re|cy_im|y_re|y_im)\(\s*(\d+)\)\s+(.*)", line)
            if m:
                key = m.group(1)
                # idx = int(m.group(2))  # 1-based, but we just append in order
                val = _parse_float(m.group(3))
                current["outputs"].setdefault(key, [])
                current["outputs"][key].append(val)

    return {"generator": generator, "machine": machine, "tests": tests}


def _parse_float(s):
    """Parse Fortran double precision format (D-exponent) to Python float.

    Handles cases where gfortran drops the 'E' for 3-digit exponents,
    e.g., '0.12345678901234567+110' -> '0.12345678901234567E+110'
    """
    s = s.strip().replace("D", "E").replace("d", "e")
    # Fix missing 'E' before +/- exponent (3-digit exponents in gfortran)
    s = re.sub(r'(\d)([+-]\d{3})$', r'\1E\2', s)
    return float(s)


def main():
    text = sys.stdin.read()
    data = parse_fortran_output(text)
    json.dump(data, sys.stdout, indent=2)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
