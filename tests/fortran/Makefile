# Makefile for generating Fortran reference values.
#
# Prerequisites:
#   - gfortran
#   - TOMS 644 source: ../../refs/fortran/toms644/zbsubs.f
#     (download via: bash ../../refs/download.sh)
#
# Usage:
#   make          — build all generators and produce JSON reference files
#   make clean    — remove compiled binaries and intermediate files

FC       = gfortran
FFLAGS   = -O2 -std=legacy
ZBSUBS   = ../../refs/fortran/toms644/zbsubs.f
MACHCON  = machcon_ieee.f
OUTDIR   = ../reference_values

# Fix zbsubs.f separator lines (leading dashes → comment lines)
ZBSUBS_FIXED = zbsubs_fixed.f

.PHONY: all clean

all: $(OUTDIR)/besk_f64.json

# Fix separator lines in zbsubs.f
$(ZBSUBS_FIXED): $(ZBSUBS)
	sed 's/^-\{2,\}/C ---/' $< > $@

# Compile reference generator
gen_besk_refs: gen_besk_refs.f $(ZBSUBS_FIXED) $(MACHCON)
	$(FC) $(FFLAGS) -o $@ gen_besk_refs.f $(ZBSUBS_FIXED) $(MACHCON)

# Generate JSON reference values
$(OUTDIR)/besk_f64.json: gen_besk_refs txt_to_json.py
	@mkdir -p $(OUTDIR)
	./gen_besk_refs | python3 txt_to_json.py > $@
	@echo "Generated $@"

clean:
	rm -f gen_besk_refs $(ZBSUBS_FIXED) gen_*.txt
