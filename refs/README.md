# Reference Materials Index

This directory contains reference materials for the complex-bessel project.
Downloaded files are **not tracked by git** — run `bash refs/download.sh` to fetch them.

---

## Sources Overview

### 1. TOMS 644 Fortran — `fortran/toms644/`

- **Origin**: [Netlib TOMS 644](https://netlib.org/toms-2014-06-10/644) (revision 930101)
- **License**: ACM copyright — algorithm reimplementation OK, direct code copying not permitted
- **Contents**:

| File | Description |
|------|-------------|
| `zbsubs.f` | Double precision subroutines (all routines in one file) |
| `cbsubs.f` | Single precision subroutines (reference only) |
| `zqcbes.f` | Double precision quick-check test programs |
| `cqcbes.f` | Single precision quick-check test programs |
| `machcon.f` | Machine constants (D1MACH, I1MACH, R1MACH) |

- **Key**: This is the **final version** with both the 1990 bug fix (K function, ln(Gamma) precision)
  and the 1995 Y function 25% speed optimization.

### 2. SLATEC Fortran — `fortran/slatec/`

- **Origin**: [Netlib SLATEC](https://www.netlib.org/slatec/) (revision 920811)
- **License**: **Public domain**
- **Contents**: 33 individual `.f` files (see mapping table below)
- **Key**: Contains 1990 bug fix but **not** the 1995 Y optimization.
  Each subroutine is a separate file, making it ideal for 1:1 mapping to Rust modules.

### 3. jpcima/zbessel C++ — `cpp/zbessel/`

- **Origin**: [github.com/jpcima/zbessel](https://github.com/jpcima/zbessel) (MIT license)
- **License**: MIT
- **Contents**:

| File | Description |
|------|-------------|
| `zbessel.cc` | Full implementation (includes `.x` fragments via `#include`) |
| `zbessel.h` | C API header |
| `zbessel.hh` | Header-only C++ API |
| `zbessel/*.x` | Per-subroutine source fragments (30 files, e.g. `zbknu.x`, `zacon.x`) |
| `zbessel/zbsubr.h` | Internal subroutine declarations |
| `zbessel/zops.h` | Complex arithmetic operations |

- **Key**: The `.x` files in `zbessel/` are per-subroutine source fragments that map 1:1
  to SLATEC Fortran files, making them ideal for side-by-side comparison when porting.
  Useful for understanding Fortran GOTO → C structured control flow conversion patterns.
  Thread-safe implementation (no static/global state).

### 4. Papers — `papers/`

| File | Description |
|------|-------------|
| `zaghloul-johnson-2025.pdf` | Zaghloul & Johnson (2025), arXiv:2505.09770 — I_nu optimization (54-80% speedup vs Amos) |
| `SAND83-0086.pdf` | Amos: "Computation of Bessel Functions of Complex Argument" (often cited as SAND83-0083) |
| `SAND83-0643.pdf` | Amos: "Computation of Bessel Functions of Complex Argument and Large Order" |

**Not available online**:
- SAND85-1018: "A Subroutine Package for Bessel Functions of a Complex Argument and Nonnegative Order"
  — Published as ACM TOMS Algorithm 644 ([ACM DL](https://dl.acm.org/doi/10.1145/7921.214331))

---

## SLATEC Fortran ↔ Rust Module Mapping

### Top-level functions (7)

| SLATEC | Fortran routine | Rust module | Description |
|--------|----------------|-------------|-------------|
| `zbesj.f` | `ZBESJ` | `besj.rs` | J_nu(z), Bessel 1st kind |
| `zbesy.f` | `ZBESY` | `besy.rs` | Y_nu(z), Bessel 2nd kind |
| `zbesi.f` | `ZBESI` | `besi.rs` | I_nu(z), modified Bessel 1st kind |
| `zbesk.f` | `ZBESK` | `besk.rs` | K_nu(z), modified Bessel 2nd kind |
| `zbesh.f` | `ZBESH` | `besh.rs` | H_nu^(m)(z), Hankel functions |
| `zairy.f` | `ZAIRY` | `airy.rs` | Ai(z), Ai'(z) |
| `zbiry.f` | `ZBIRY` | `airy.rs` | Bi(z), Bi'(z) |

### Core algorithms (21)

| SLATEC | Fortran routine | Rust module | Description |
|--------|----------------|-------------|-------------|
| `zbknu.f` | `ZBKNU` | `algo/bknu.rs` | K function core computation |
| `zbinu.f` | `ZBINU` | `algo/binu.rs` | I function core (series/Miller/asymptotic dispatch) |
| `zacon.f` | `ZACON` | `algo/acon.rs` | Analytic continuation (Re(z) < 0) |
| `zacai.f` | `ZACAI` | `algo/acon.rs` | Analytic continuation for Airy (simplified ZACON) |
| `zbunk.f` | `ZBUNK` | `algo/bunk.rs` | K uniform asymptotic dispatch |
| `zbuni.f` | `ZBUNI` | `algo/buni.rs` | I/J uniform asymptotic dispatch |
| `zseri.f` | `ZSERI` | `algo/seri.rs` | Power series (small \|z\|) |
| `zmlri.f` | `ZMLRI` | `algo/mlri.rs` | Miller algorithm |
| `zasyi.f` | `ZASYI` | `algo/asyi.rs` | Asymptotic expansion (large \|z\|) |
| `zuoik.f` | `ZUOIK` | `algo/uoik.rs` | Overflow/underflow test |
| `zuchk.f` | `ZUCHK` | `algo/uchk.rs` | Underflow check |
| `zkscl.f` | `ZKSCL` | `algo/kscl.rs` | K function scaling |
| `zshch.f` | `ZSHCH` | `algo/shch.rs` | sinh/cosh |
| `zs1s2.f` | `ZS1S2` | `algo/s1s2.rs` | Stokes multiplier |
| `zwrsk.f` | `ZWRSK` | `algo/wrsk.rs` | Wronskian relation |
| `zrati.f` | `ZRATI` | `algo/rati.rs` | Continued fraction ratios |
| `zuni1.f` | `ZUNI1` | `algo/uni1.rs` | I/J uniform asymptotic, region 1 |
| `zuni2.f` | `ZUNI2` | `algo/uni2.rs` | I/J uniform asymptotic, region 2 |
| `zunik.f` | `ZUNIK` | `algo/unik.rs` | K Airy parameters |
| `zunhj.f` | `ZUNHJ` | `algo/unhj.rs` | H/J Airy parameters |
| `zunk1.f` | `ZUNK1` | `algo/unk1.rs` | K uniform asymptotic, region 1 |
| `zunk2.f` | `ZUNK2` | `algo/unk2.rs` | K uniform asymptotic, region 2 |

### Utilities (9)

| SLATEC | Fortran routine | Rust module | Description |
|--------|----------------|-------------|-------------|
| `zabs.f` | `ZABS` | `utils.rs` | \|z\| (overflow-safe) |
| `zdiv.f` | `ZDIV` | `utils.rs` | z1/z2 (overflow-safe) |
| `zexp.f` | `ZEXP` | *(intrinsic)* | Complex exp — `Complex::exp()` in Rust |
| `zlog.f` | `ZLOG` | *(intrinsic)* | Complex log — `Complex::ln()` in Rust |
| `zmlt.f` | `ZMLT` | *(intrinsic)* | Complex multiply — `*` operator in Rust |
| `zsqrt.f` | `ZSQRT` | *(intrinsic)* | Complex sqrt — `Complex::sqrt()` in Rust |
| `dgamln.f` | `DGAMLN` | `algo/gamln.rs` | ln(Gamma(x)) |
| `d1mach.f` | `D1MACH` | `machine.rs` | Double precision machine constants |
| `i1mach.f` | `I1MACH` | `machine.rs` | Integer machine constants |

---

## TOMS 644 vs SLATEC Differences

Both derive from Amos's original code, but differ in revision date:

| Aspect | TOMS 644 (930101) | SLATEC (920811) |
|--------|-------------------|-----------------|
| 1990 bug fix | Included | Included |
| 1995 Y optimization | **Included** | **Not included** |
| File structure | Single combined file per precision | Individual files per subroutine |
| License | ACM copyright | Public domain |

### 1995 Y Optimization Details

The optimization affects `ZBESY` (Y function), providing ~25% speed improvement.
To identify the exact changes:

```bash
# Extract ZBESY from TOMS 644 combined file, then diff
grep -A 9999 'SUBROUTINE ZBESY' refs/fortran/toms644/zbsubs.f | \
  sed '/SUBROUTINE ZB[^E]/,$d' > /tmp/toms_zbesy.f
diff refs/fortran/slatec/zbesy.f /tmp/toms_zbesy.f
```

For a broader comparison of all routines:

```bash
# Compare any specific routine (example: ZBKNU)
grep -A 9999 'SUBROUTINE ZBKNU' refs/fortran/toms644/zbsubs.f | \
  sed '/^      SUBROUTINE Z[^B]/,$d; /^      SUBROUTINE ZB[^K]/,$d' > /tmp/toms_zbknu.f
diff refs/fortran/slatec/zbknu.f /tmp/toms_zbknu.f
```

**Recommendation**: Use TOMS 644 `zbsubs.f` as the primary reference for algorithm correctness
(latest version), and SLATEC individual files for understanding module boundaries.

---

## Phase-by-Phase Reference Guide

### Phase 1: Foundation

| Task | Primary reference | Secondary reference |
|------|------------------|-------------------|
| `machine.rs` (BesselFloat) | `slatec/d1mach.f`, `slatec/i1mach.f` | `toms644/machcon.f` |
| `utils.rs` (zabs, zdiv) | `slatec/zabs.f`, `slatec/zdiv.f` | `zbessel.cc` (C++ patterns) |

### Phase 2: K Function Core

| Task | Primary reference | Secondary reference |
|------|------------------|-------------------|
| `algo/gamln.rs` | `slatec/dgamln.f` | `toms644/zbsubs.f` (DGAMLN section) |
| `algo/bknu.rs` | `slatec/zbknu.f` | `toms644/zbsubs.f` (ZBKNU section) |
| `algo/shch.rs` | `slatec/zshch.f` | — |
| `algo/s1s2.rs` | `slatec/zs1s2.f` | — |
| `algo/kscl.rs` | `slatec/zkscl.f` | — |
| `algo/uchk.rs` | `slatec/zuchk.f` | — |
| `besk.rs` | `slatec/zbesk.f` | `toms644/zbsubs.f` (ZBESK section) |

### Phase 3: H, J, Y Functions

| Task | Primary reference | Secondary reference |
|------|------------------|-------------------|
| `algo/wrsk.rs` | `slatec/zwrsk.f` | — |
| `besh.rs` | `slatec/zbesh.f` | `toms644/zbsubs.f` (ZBESH section) |
| `besj.rs` | `slatec/zbesj.f` | `toms644/zbsubs.f` (ZBESJ section) |
| `besy.rs` | **`toms644/zbsubs.f`** (ZBESY) | `slatec/zbesy.f` (diff for 1995 changes) |

### Phase 4: I Function + Analytic Continuation

| Task | Primary reference | Secondary reference |
|------|------------------|-------------------|
| `algo/seri.rs` | `slatec/zseri.f` | — |
| `algo/mlri.rs` | `slatec/zmlri.f` | — |
| `algo/asyi.rs` | `slatec/zasyi.f` | — |
| `algo/rati.rs` | `slatec/zrati.f` | — |
| `algo/uoik.rs` | `slatec/zuoik.f` | — |
| `algo/binu.rs` | `slatec/zbinu.f` | `zbessel.cc` (control flow) |
| `besi.rs` | `slatec/zbesi.f` | — |
| `algo/acon.rs` | `slatec/zacon.f` | `zbessel.cc` (control flow) |

### Phase 5: Uniform Asymptotic Expansion

| Task | Primary reference | Secondary reference |
|------|------------------|-------------------|
| `algo/unik.rs` | `slatec/zunik.f` | SAND83-0643 (theory) |
| `algo/unhj.rs` | `slatec/zunhj.f` | SAND83-0643 (theory) |
| `algo/uni1.rs` | `slatec/zuni1.f` | — |
| `algo/uni2.rs` | `slatec/zuni2.f` | — |
| `algo/buni.rs` | `slatec/zbuni.f` | — |
| `algo/unk1.rs` | `slatec/zunk1.f` | — |
| `algo/unk2.rs` | `slatec/zunk2.f` | — |
| `algo/bunk.rs` | `slatec/zbunk.f` | — |

### Phase 6: Airy Functions

| Task | Primary reference | Secondary reference |
|------|------------------|-------------------|
| `airy.rs` (zairy) | `slatec/zairy.f` | `toms644/zbsubs.f` (ZAIRY section) |
| `airy.rs` (zbiry) | `slatec/zbiry.f` | `toms644/zbsubs.f` (ZBIRY section) |

### Phase 8+: I_nu Optimization (optional)

| Task | Primary reference |
|------|------------------|
| I_nu algorithm replacement | `papers/zaghloul-johnson-2025.pdf` |

---

## Call Relationship Diagram

```
zbesj ──→ zbinu ─┬→ zseri   (power series, small |z|)
  │              ├→ zmlri   (Miller algorithm)
  │              ├→ zasyi   (asymptotic, large |z|)
  │              └→ zuoik   (overflow/underflow test)
  ├───→ zbuni ─┬→ zuni1 ─→ zunik, zunhj
  │            └→ zuni2 ─→ zunik, zunhj
  └───→ zbknu              (K core, Wronskian path)

zbesy ──→ zbesh             (via Hankel functions)

zbesi ──→ zbinu             (series/Miller/asymptotic)
  ├───→ zbuni               (uniform asymptotic)
  └───→ zuoik

zbesk ──→ zbknu             (K core, Re(z) >= 0)
  ├───→ zacon ─→ zbinu      (analytic continuation, Re(z) < 0)
  ├───→ zbunk ─┬→ zunk1 ─→ zunik
  │            └→ zunk2 ─→ zunik
  └───→ zuoik

zbesh ──→ zbknu             (K core, direct)
  ├───→ zacon               (analytic continuation)
  ├───→ zbunk               (uniform asymptotic)
  └───→ zuoik

zairy ──→ zseri / zbknu / zuni1 / zuni2 / zunik / zunhj
zbiry ──→ zbinu / zbknu / zuni1 / zuni2 / zunik / zunhj

Common utilities used throughout:
  zabs    — overflow-safe |z|
  zdiv    — overflow-safe complex division
  dgamln  — ln(Gamma(x))
  zuchk   — underflow check
  zshch   — sinh/cosh
  zs1s2   — Stokes multiplier
  zkscl   — K function scaling
  zrati   — continued fraction ratios
  zwrsk   — Wronskian relation
```

---

## SAND Reports

| Report | Title | Available | URL |
|--------|-------|-----------|-----|
| SAND83-0086 | Computation of Bessel Functions of Complex Argument | PDF downloaded | [OSTI](https://doi.org/10.2172/6148453) |
| SAND83-0643 | Computation of Bessel Functions of Complex Argument and Large Order | PDF downloaded | [OSTI](https://doi.org/10.2172/5903937) |
| SAND85-1018 | A Subroutine Package for Bessel Functions of a Complex Argument and Nonnegative Order | **Not available** | Published as [ACM TOMS Alg. 644](https://dl.acm.org/doi/10.1145/7921.214331) |

> **Note**: SAND83-0086 is often cited as SAND83-0083 in older source code comments.
> The correct number per OSTI records is SAND83-0086.
