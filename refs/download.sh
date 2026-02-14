#!/bin/bash
# download.sh — Download reference materials for complex-bessel project
#
# Usage: bash refs/download.sh
#        (run from project root)
#
# Downloads:
#   1. TOMS 644 Fortran (5 files, split from combined archive)
#   2. SLATEC Fortran (31 individual files)
#   3. jpcima/zbessel C++ (2 files)
#   4. Zaghloul-Johnson 2025 paper (PDF)
#   5. Amos SAND reports (2 available PDFs)

set -euo pipefail

REFS_DIR="$(cd "$(dirname "$0")" && pwd)"
CURL_OPTS="-fsSL --retry 3 --retry-delay 2"

# Color output helpers
green() { printf "\033[32m%s\033[0m\n" "$1"; }
yellow() { printf "\033[33m%s\033[0m\n" "$1"; }
red() { printf "\033[31m%s\033[0m\n" "$1"; }

download() {
    local url="$1"
    local dest="$2"
    if [ -f "$dest" ]; then
        yellow "  SKIP (exists): $(basename "$dest")"
        return 0
    fi
    if curl $CURL_OPTS -o "$dest" "$url" 2>/dev/null; then
        green "  OK: $(basename "$dest")"
        return 0
    else
        red "  FAIL: $(basename "$dest") ($url)"
        return 1
    fi
}

# Track failures
FAILURES=0

# ============================================================
# 1. TOMS 644 Fortran (combined file → split into 5 files)
# ============================================================
echo ""
echo "=== TOMS 644 Fortran (netlib, revision 930101) ==="
TOMS_DIR="$REFS_DIR/fortran/toms644"
mkdir -p "$TOMS_DIR"

TOMS_COMBINED="$TOMS_DIR/.combined_644"
if [ -f "$TOMS_DIR/zbsubs.f" ] && [ -f "$TOMS_DIR/cbsubs.f" ] && \
   [ -f "$TOMS_DIR/zqcbes.f" ] && [ -f "$TOMS_DIR/cqcbes.f" ] && \
   [ -f "$TOMS_DIR/machcon.f" ]; then
    yellow "  SKIP (all files exist): toms644/"
else
    echo "  Downloading combined archive..."
    if curl $CURL_OPTS -o "$TOMS_COMBINED" \
       "https://netlib.org/toms-2014-06-10/644" 2>/dev/null; then
        green "  Downloaded combined file"

        # Split by "C*** filename" markers
        echo "  Splitting into individual files..."
        awk -v outdir="$TOMS_DIR" '
        /^C\*\*\* / {
            fname = $2
            # Skip the Readme
            if (fname == "Readme") { skip = 1; next }
            skip = 0
            next
        }
        !skip && fname { print > (outdir "/" fname) }
        ' "$TOMS_COMBINED"

        # Verify expected files
        for f in zbsubs.f cbsubs.f zqcbes.f cqcbes.f machcon.f; do
            if [ -f "$TOMS_DIR/$f" ]; then
                green "  OK: $f"
            else
                red "  FAIL: $f not created from split"
                FAILURES=$((FAILURES + 1))
            fi
        done

        rm -f "$TOMS_COMBINED"
    else
        red "  FAIL: Could not download combined TOMS 644 archive"
        FAILURES=$((FAILURES + 1))
    fi
fi

# ============================================================
# 2. SLATEC Fortran (31 individual files)
# ============================================================
echo ""
echo "=== SLATEC Fortran (public domain, revision 920811) ==="
SLATEC_DIR="$REFS_DIR/fortran/slatec"
mkdir -p "$SLATEC_DIR"

SLATEC_BASE="https://www.netlib.org/slatec/src"

# Top-level functions (7)
SLATEC_TOP="zbesj zbesy zbesi zbesk zbesh zairy zbiry"

# Core algorithms (21 — includes zacai for zairy's analytic continuation)
SLATEC_CORE="zbknu zbinu zacon zacai zbunk zbuni zseri zmlri zasyi zuoik zuchk
             zkscl zshch zs1s2 zwrsk zrati zuni1 zuni2 zunik zunhj zunk1 zunk2"

# Utilities (9 — includes complex arithmetic wrappers)
SLATEC_UTIL="zabs zdiv zexp zlog zmlt zsqrt dgamln d1mach i1mach"

for name in $SLATEC_TOP $SLATEC_CORE $SLATEC_UTIL; do
    download "$SLATEC_BASE/${name}.f" "$SLATEC_DIR/${name}.f" || FAILURES=$((FAILURES + 1))
done

# ============================================================
# 3. jpcima/zbessel C++ (root files + zbessel/ subroutines)
# ============================================================
echo ""
echo "=== jpcima/zbessel C++ port ==="
CPP_DIR="$REFS_DIR/cpp/zbessel"
CPP_SUB_DIR="$CPP_DIR/zbessel"
mkdir -p "$CPP_SUB_DIR"

ZBESSEL_BASE="https://raw.githubusercontent.com/jpcima/zbessel/master"

# Root-level files
download "$ZBESSEL_BASE/zbessel.cc" "$CPP_DIR/zbessel.cc" || FAILURES=$((FAILURES + 1))
download "$ZBESSEL_BASE/zbessel.h"  "$CPP_DIR/zbessel.h"  || FAILURES=$((FAILURES + 1))
download "$ZBESSEL_BASE/zbessel.hh" "$CPP_DIR/zbessel.hh" || FAILURES=$((FAILURES + 1))

# zbessel/ subroutine files (.x = included source fragments, .h = internal headers)
ZBESSEL_SUBS="zacai zacon zairy zasyi zbesh zbesi zbesj zbesk zbesy zbinu zbiry
              zbknu zbuni zbunk zkscl zmlri zops zrati zs1s2 zseri zshch zuchk
              zunhj zuni1 zuni2 zunik zunk1 zunk2 zuoik zwrsk"

for name in $ZBESSEL_SUBS; do
    download "$ZBESSEL_BASE/zbessel/${name}.x" "$CPP_SUB_DIR/${name}.x" || FAILURES=$((FAILURES + 1))
done
download "$ZBESSEL_BASE/zbessel/zbsubr.h" "$CPP_SUB_DIR/zbsubr.h" || FAILURES=$((FAILURES + 1))
download "$ZBESSEL_BASE/zbessel/zops.h"   "$CPP_SUB_DIR/zops.h"   || FAILURES=$((FAILURES + 1))

# ============================================================
# 4. Zaghloul-Johnson 2025 paper
# ============================================================
echo ""
echo "=== Zaghloul-Johnson 2025 paper (arXiv:2505.09770) ==="
PAPERS_DIR="$REFS_DIR/papers"
mkdir -p "$PAPERS_DIR"
download "https://arxiv.org/pdf/2505.09770" \
         "$PAPERS_DIR/zaghloul-johnson-2025.pdf" || FAILURES=$((FAILURES + 1))

# ============================================================
# 5. Amos SAND reports (try available ones)
# ============================================================
echo ""
echo "=== Amos SAND reports ==="

# SAND83-0086 (often cited as SAND83-0083)
# "Computation of Bessel Functions of Complex Argument"
download "https://www.osti.gov/servlets/purl/6148453" \
         "$PAPERS_DIR/SAND83-0086.pdf" || {
    yellow "  NOTE: SAND83-0086 not available — URL recorded in README"
    FAILURES=$((FAILURES + 1))
}

# SAND83-0643
# "Computation of Bessel Functions of Complex Argument and Large Order"
download "https://www.osti.gov/servlets/purl/5903937" \
         "$PAPERS_DIR/SAND83-0643.pdf" || {
    yellow "  NOTE: SAND83-0643 not available — URL recorded in README"
    FAILURES=$((FAILURES + 1))
}

# SAND85-1018
# "A Subroutine Package for Bessel Functions of a Complex Argument and Nonnegative Order"
# Not available online as PDF
yellow "  NOTE: SAND85-1018 not available online — see README for details"

# ============================================================
# Summary
# ============================================================
echo ""
echo "========================================="
if [ $FAILURES -eq 0 ]; then
    green "All downloads completed successfully."
else
    yellow "$FAILURES download(s) had issues. Check output above."
fi
echo "========================================="
