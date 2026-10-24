# Lambda_DH bracket note.  `make` builds the PDF; `make verify` checks every number.
PY := ../../../.venv/bin/python

.PHONY: all verify clean

all: lambda_note.pdf

lambda_note.pdf: lambda_note.tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error lambda_note.tex

# re-derive every number printed in the note and fail on any disagreement, any
# undeclared literal, any overstatement of precision, or an unmet prerequisite
verify:
	$(PY) verify_note_numbers.py

clean:
	latexmk -C lambda_note.tex
