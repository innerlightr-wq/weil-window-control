# weil-window-control
# Code: MIT (LICENSE).  Paper and documentation: CC BY 4.0 (LICENSE-CC-BY-4.0).

PY ?= python3

.PHONY: help quick stage1 stage2 gates connes law onset figures paper all clean comb-dips-quick resolution-inheritance-quick

help:
	@echo "make quick    - self-test, no high-precision runs        (~2 min)"
	@echo "make stage1   - function-field control, genuine spectra  (~2 min)"
	@echo "make stage2   - planted spectra, exact over Q            (~6 min)"
	@echo "make gates    - zeta Gate A and Gate B                   (~45 min)"
	@echo "make connes   - Connes [1,13] window, even/odd, split    (~60 min)"
	@echo "make law      - delta^2-L^3 law: proofs and constants    (~25 min)"
	@echo "make onset    - detection onset vs L_pred and basis      (HOURS)"
	@echo "make paper    - compile paper/weil_window_control.pdf    (~10 s)"
	@echo "make all      - everything except 'onset'"
	@echo ""
	@echo "extensions (not part of the technical note):"
	@echo "make comb-dips-quick            - comb-dips gate + 2 Stage 2 rows   (~15 s)"
	@echo "make resolution-inheritance-quick - Weil-form evaluator validation  (~6 s)"

quick:   ; ./scripts/quick_check.sh
stage1:  ; ./scripts/01_stage1_function_field.sh
stage2:  ; ./scripts/02_stage2_planted.sh
gates:   ; ./scripts/03_stage3_zeta_gates.sh
connes:  ; ./scripts/04_stage3_connes_window.sh
law:     ; ./scripts/05_delta2_L3_law.sh
onset:   ; ./scripts/06_detection_onset.sh

comb-dips-quick: ; ./extensions/comb-dips/scripts/quick_check.sh
resolution-inheritance-quick: ; $(PY) extensions/resolution-inheritance/src/validate.py

paper:
	cd paper && pdflatex -interaction=nonstopmode weil_window_control.tex >/dev/null && \
	           pdflatex -interaction=nonstopmode weil_window_control.tex >/dev/null && \
	           echo "built paper/weil_window_control.pdf"

all: stage1 stage2 gates connes law paper

clean:
	rm -rf src/__pycache__ extensions/*/src/__pycache__
	cd paper && rm -f *.aux *.log *.out *.toc
