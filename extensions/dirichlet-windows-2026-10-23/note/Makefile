# Parity-blindness note.  `make` builds the PDF; `make verify` checks every number.
PY := ../../../.venv/bin/python

.PHONY: all verify tables clean

all: siegel_note.pdf

siegel_note.pdf: siegel_note.tex eig_table.tex pole_table.tex onset_table.tex ff_table.tex K_table.tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error siegel_note.tex

# regenerate the machine-generated tables from the extension's data files
tables:
	$(PY) make_tables.py

# re-derive every number printed in the note and fail on any disagreement, any
# undeclared literal, any overstatement of precision, or an unmet prerequisite
verify:
	$(PY) verify_note_numbers.py

clean:
	latexmk -C siegel_note.tex
