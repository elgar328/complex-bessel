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

all: $(OUTDIR)/besk_f64.json $(OUTDIR)/besh_f64.json $(OUTDIR)/besi_f64.json $(OUTDIR)/besj_f64.json $(OUTDIR)/besy_f64.json $(OUTDIR)/airy_f64.json $(OUTDIR)/phase5b_f64.json

# Fix separator lines in zbsubs.f
$(ZBSUBS_FIXED): $(ZBSUBS)
	sed 's/^-\{2,\}/C ---/' $< > $@

# ── Phase 2: K function ──

gen_besk_refs: gen_besk_refs.f $(ZBSUBS_FIXED) $(MACHCON)
	$(FC) $(FFLAGS) -o $@ gen_besk_refs.f $(ZBSUBS_FIXED) $(MACHCON)

$(OUTDIR)/besk_f64.json: gen_besk_refs txt_to_json.py
	@mkdir -p $(OUTDIR)
	./gen_besk_refs | python3 txt_to_json.py > $@
	@echo "Generated $@"

# ── Phase 3: H function ──

gen_besh_refs: gen_besh_refs.f $(ZBSUBS_FIXED) $(MACHCON)
	$(FC) $(FFLAGS) -o $@ gen_besh_refs.f $(ZBSUBS_FIXED) $(MACHCON)

$(OUTDIR)/besh_f64.json: gen_besh_refs txt_to_json.py
	@mkdir -p $(OUTDIR)
	./gen_besh_refs | python3 txt_to_json.py > $@
	@echo "Generated $@"

# ── Phase 4: I function ──

gen_besi_refs: gen_besi_refs.f $(ZBSUBS_FIXED) $(MACHCON)
	$(FC) $(FFLAGS) -o $@ gen_besi_refs.f $(ZBSUBS_FIXED) $(MACHCON)

$(OUTDIR)/besi_f64.json: gen_besi_refs txt_to_json.py
	@mkdir -p $(OUTDIR)
	./gen_besi_refs | python3 txt_to_json.py > $@
	@echo "Generated $@"

# ── Phase 4: J function ──

gen_besj_refs: gen_besj_refs.f $(ZBSUBS_FIXED) $(MACHCON)
	$(FC) $(FFLAGS) -o $@ gen_besj_refs.f $(ZBSUBS_FIXED) $(MACHCON)

$(OUTDIR)/besj_f64.json: gen_besj_refs txt_to_json.py
	@mkdir -p $(OUTDIR)
	./gen_besj_refs | python3 txt_to_json.py > $@
	@echo "Generated $@"

# ── Phase 4: Y function ──

gen_besy_refs: gen_besy_refs.f $(ZBSUBS_FIXED) $(MACHCON)
	$(FC) $(FFLAGS) -o $@ gen_besy_refs.f $(ZBSUBS_FIXED) $(MACHCON)

$(OUTDIR)/besy_f64.json: gen_besy_refs txt_to_json.py
	@mkdir -p $(OUTDIR)
	./gen_besy_refs | python3 txt_to_json.py > $@
	@echo "Generated $@"

# ── Phase 6: Airy functions ──

gen_airy_refs: gen_airy_refs.f $(ZBSUBS_FIXED) $(MACHCON)
	$(FC) $(FFLAGS) -o $@ gen_airy_refs.f $(ZBSUBS_FIXED) $(MACHCON)

$(OUTDIR)/airy_f64.json: gen_airy_refs txt_to_json.py
	@mkdir -p $(OUTDIR)
	./gen_airy_refs | python3 txt_to_json.py > $@
	@echo "Generated $@"

# ── Phase 5b: Large-order tests ──

gen_phase5b_refs: gen_phase5b_refs.f $(ZBSUBS_FIXED) $(MACHCON)
	$(FC) $(FFLAGS) -o $@ gen_phase5b_refs.f $(ZBSUBS_FIXED) $(MACHCON)

$(OUTDIR)/phase5b_f64.json: gen_phase5b_refs txt_to_json.py
	@mkdir -p $(OUTDIR)
	./gen_phase5b_refs | python3 txt_to_json.py > $@
	@echo "Generated $@"

clean:
	rm -f gen_besk_refs gen_besh_refs gen_besi_refs gen_besj_refs gen_besy_refs gen_airy_refs gen_phase5b_refs $(ZBSUBS_FIXED) gen_*.txt
