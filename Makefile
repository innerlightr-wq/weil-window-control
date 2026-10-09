# weil-window-control
# Code: MIT (LICENSE).  Paper and documentation: CC BY 4.0 (LICENSE-CC-BY-4.0).

PY ?= python3

.PHONY: help quick stage1 stage2 gates connes law onset figures paper all clean comb-dips-quick \
        witness-quick witness-cold witness-independent

help:
	@echo "make quick    - self-test, no high-precision runs        (~2 min)"
	@echo "make stage1   - function-field control, genuine spectra  (~2 min)"
	@echo "make stage2   - planted spectra, exact over Q            (~6 min)"
	@echo "make gates    - zeta Gate A and Gate B                   (~45 min)"
	@echo "make connes   - Connes [1,13] window, even/odd, split    (~60 min)"
	@echo "make law      - delta^2-L^3 law: proofs and constants    (~25 min)"
	@echo "make onset    - detection onset vs L_pred and basis      (HOURS)"
	@echo "make paper    - compile the CURRENT paper PDF             (~10 s)"
	@echo "make all      - everything except 'onset'"
	@echo ""
	@echo ""
	@echo "certified four-mode witness (paper S6.1, Thm 11; standard library only):"
	@echo "make witness-quick       - interval + cache audits         (~50 s)"
	@echo "make witness-cold        - COLD certificate, no cache      (~5.5 min)"
	@echo "make witness-independent - separate implementation         (~3.5 min)"
	@echo ""
	@echo "extensions (not part of the technical note):"
	@echo "make comb-dips-quick - comb-dips gate + 2 Stage 2 rows   (~15 s)"

quick:   ; ./scripts/quick_check.sh
stage1:  ; ./scripts/01_stage1_function_field.sh
stage2:  ; ./scripts/02_stage2_planted.sh
gates:   ; ./scripts/03_stage3_zeta_gates.sh
connes:  ; ./scripts/04_stage3_connes_window.sh
law:     ; ./scripts/05_delta2_L3_law.sh
onset:   ; ./scripts/06_detection_onset.sh

comb-dips-quick: ; ./extensions/comb-dips/scripts/quick_check.sh

# --- certified four-mode witness (paper S6.1, Theorem 11) ----------------------
# Standard library only: no mpmath, no numpy. Each target exits non-zero on failure.
WITNESS := extensions/planted-witness-consolidation-2026-10-08

witness-quick:
	cd $(WITNESS) && $(PY) audit_primitives.py && $(PY) audit_cache_invalidation.py

# Deletes the cache first, so the certificate is genuinely recomputed from the
# frozen rational input and cannot reuse a stored value.
witness-cold:
	cd $(WITNESS) && rm -f regenerated_qzeta_cache.json cold_certificate.json \
	    interval_cells.json && $(PY) verify_certificate_consolidated.py

witness-independent:
	cd $(WITNESS) && $(PY) independent_point_check.py

# Builds the CURRENT paper under its own name. paper/weil_window_control.pdf is the
# previous, superseded version and is deliberately not regenerated here.
PAPER_JOB := weil_window_control_rev_certified_witness
paper:
	cd paper && pdflatex -interaction=nonstopmode -jobname=$(PAPER_JOB) weil_window_control.tex >/dev/null && \
	           pdflatex -interaction=nonstopmode -jobname=$(PAPER_JOB) weil_window_control.tex >/dev/null && \
	           echo "built paper/$(PAPER_JOB).pdf"

all: stage1 stage2 gates connes law paper

clean:
	rm -rf src/__pycache__ extensions/comb-dips/src/__pycache__
	cd paper && rm -f *.aux *.log *.out *.toc *.fls *.fdb_latexmk
