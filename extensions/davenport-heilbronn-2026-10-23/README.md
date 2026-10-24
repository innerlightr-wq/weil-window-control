# Davenport–Heilbronn control for the finite-window Weil detector

Extension of [`weil-window-control`](../../README.md). **Carries no DOI, is not part of
the paper, and has not been peer-reviewed.** Non-claims, evidence tiers, full lab record
and retraction log: [`REPORT.md`](REPORT.md). The derivation of everything computed here:
[`DERIVATION.md`](DERIVATION.md). What was and was not read: [`SOURCE_AUDIT.md`](SOURCE_AUDIT.md).

---

## What this is

The repository's finite-window Weil detector, run against the **actual
Davenport–Heilbronn function** instead of a planted spectrum.

`F(s) = [(1−iκ)L(s,χ) + (1+iκ)L(s,χ̄)]/2` with `χ` the odd primitive character mod 5 and
`κ` **derived from the Gauss sum** is the classical example of a Dirichlet series with a
`ζ`-shaped functional equation, **no Euler product, and zeros off the critical line**. It
is the one object on which Weil positivity is *supposed* to fail, so it is the one object
that can tell a working detector from a detector that only ever sees what was planted in
it.

The deliverable is **one explicitly specified compactly supported test function** whose
Davenport–Heilbronn Weil value is negative, certified in ball arithmetic, together with
the decomposition that separates the selected off-line quartet from everything else.

## Reproducing

See [`REPRODUCE.md`](REPRODUCE.md). In short, from the repository root:

```bash
.venv/bin/pip install 'mpmath==1.3.0' python-flint
./extensions/davenport-heilbronn-2026-10-23/run_suite.sh
```

Python **3.14.7**, mpmath **1.3.0**, python-flint **0.9.0** (Arb). No numpy, no scipy, no
float64 anywhere on a path that reaches a reported value.

## The note

A standalone 12-page note, [`note/dh_note.pdf`](note/dh_note.pdf), states the one result
with its non-claims, tiers and retraction log:
**a certified negative finite-window Weil witness for the Davenport-Heilbronn function,
with the remaining-zero form positive.** `make -C note verify` re-derives every number
printed in it from the data files and fails on any disagreement, any undeclared literal,
or any overstatement of precision.

## Layout

```
note/           dh_note.tex -> dh_note.pdf, its generated tables, and the verifier
REPORT.md       the lab record: non-claims first, a tier on every claim, verdict last
DERIVATION.md   every identity that is computed, derived, with the check that closes it
SOURCE_AUDIT.md what was inspected, what was only a bibliographic record, what failed
RETRACTIONS.md  append-only log of what was got wrong and how it was found
REPRODUCE.md    dependency order, commands, which script makes which file
src/            dh.py (the function and the form) · certify.py (the Arb re-implementation)
                census.py · phaseA..G.py · orbitform.py · config.py
tests/          test_orbits.py (orbit conventions) · regression.py (parent repository)
data/           CSV results
certificates/   the frozen witness and the ball certificate
results/        console transcripts and JSON summaries
logs/           the runs that produced the committed results
```

## Conventions, once

`g` real, supported in `[−L, L]`; `G(z) = ∫g(u)e^{zu}du`; `h(v) = ∫g(u+v)g(u)du`;
`H(z) = G(z)G(−z)`. The spectral form runs over the zeros of the **entire completed**
`Λ_F(s) = (5/π)^{(s+1)/2}Γ((s+1)/2)F(s)`, which satisfies `Λ_F(s) = Λ_F(1−s)`:

```
Q_F(g) = Σ_ρ mult(ρ) H(ρ − ½)
       = h(0) log(5/π) + (1/2π)∫|G(it)|² Re ψ(¾+it/2) dt − 2 Σ_{2≤n≤e^{2L}} b_n h(log n)/√n
```

with `−F′/F(s) = Σ b_n n^{−s}`. **No pole term**: `Λ_F` is entire. The window half-length
`L` is the support half-length of **`g`**, so `h` is supported in `[−2L, 2L]`; a reader
comparing with another convention should check which of the two is meant.

## Licence

As the parent repository: **MIT** for code under `src/`, `tests/` and `run_suite.sh`,
**CC BY 4.0** for the Markdown, `data/`, `certificates/` and `results/`.
