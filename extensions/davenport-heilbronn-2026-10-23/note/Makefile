# Davenport-Heilbronn note.  `make` builds the PDF; `make verify` checks every number.
PY := ../../../.venv/bin/python

.PHONY: all verify tables clean

all: dh_note.pdf

dh_note.pdf: dh_note.tex witness_table.tex census_table.tex remaining_table.tex onset_table.tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error dh_note.tex

# regenerate the machine-generated tables from the data files
tables:
	$(PY) make_witness_table.py
	$(PY) make_tables.py

# re-derive every number printed in the note and fail on any disagreement,
# any undeclared literal, or any overstatement of precision
verify:
	$(PY) verify_note_numbers.py

clean:
	latexmk -C dh_note.tex
