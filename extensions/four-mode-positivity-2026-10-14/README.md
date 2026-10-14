# Four-mode Weil positivity, and an extended added-quartet detection interval

`extensions/four-mode-positivity-2026-10-14/` — follow-on work to §6.1 /
Theorem 11 of the note. Preregistered and tiered the same way as the rest of
`extensions/`; **no DOI of its own, no external peer review.**

Three certified theorems and one corrected diagnostic, all at `L = 4/5`,
`gamma = 14`, in the span of four `EVEN_D` modes.

> ### Scope, before the results
>
> The quartet of zeros at `1/2 ± delta ± 14i` is **hypothetical**: it is *added*
> to the true zero multiset of `zeta`, and nothing is removed. **`gamma = 14` is
> a prescribed rational height, not an observed ordinate of `zeta`.** No zero of
> `zeta` is located, constrained or verified anywhere here.
>
> **Four-mode positivity is not universal Weil positivity.** Weil positivity —
> `Q_zeta(f) >= 0` for *all* admissible `f` — is equivalent to RH; a
> four-dimensional subspace at one window length is a measure-zero slice of it,
> and the sign found there is the sign RH predicts. **Nothing here proves,
> weakens, or conditionally establishes the Riemann hypothesis.**
>
> Likewise, **non-detection in a finite basis is not positivity**: a larger
> basis could detect at smaller `delta`, and nothing here excludes that.

---

## What is certified

### Theorem A — detection down to `delta = 9/100`
[`CERTIFICATE_THEOREM.md`](CERTIFICATE_THEOREM.md)

```
   Q_add,delta(f_x) < 0   for EVERY delta in [9/100, 1/4],
```
by **one** frozen rational witness `x in Q^4` over **69 closed rational cells**,
each certified *as a cell* by interval arithmetic in `delta`, not by sampling.
The join at `delta = 1/4` is certified **from both sides**, so with the
inherited published certificate on `[1/4, 49/100]` the certified detection
domain is

```
   delta  in  [ 9/100 , 49/100 ] ,
```
extending the previously certified left endpoint from `1/4` down to `9/100`.

### Theorem O1 — unconditional four-mode positivity
[`FOUR_MODE_POSITIVITY.md`](FOUR_MODE_POSITIVITY.md)

```
   G >= m0 I   on the four-mode subspace,   m0 >= 1.659599e-10 > 0 ,
```
i.e. `Q_zeta(f_x) >= m0 ||x||^2` for every `x in R^4`. **Unconditional**: derived
from the explicit formula alone, with no RH assumption and **no external
positivity input**. The complete matrix `G` is used, with the archimedean tail
included — not the tail-free `G_core`; entries are rigorously enclosed by
polarization, the row-sum loss is bounded in dimension 4, and the eigenvalue
bound is an exact rational `LDL^T` computation. Consequently
`H(delta) = G + 4aa^T >= m0 I` for **every** `delta`.

The spectrum of the enclosure midpoint is
`1.746e-10, 2.711e-06, 3.846e-03, 1.528` — condition number `8.75e+09`. The
form is **nearly singular**, which is why the error budget had to be sharpened
twice.

### Theorem O2 — certified non-detection on `[0, 0.05986]`
[`NON_DETECTION.md`](NON_DETECTION.md)

```
   Q_add,delta(f_x) >= 0   for EVERY x in the subspace and every
   delta in [0, 3831/64000] = [0, 0.05985937…] ,
```
over **6072 closed rational cells**, consecutive and gap-free from `0`. **No
negative witness exists in this subspace there.** It needs Theorem O1 as input,
and an AM–GM reformulation of the cell loss to be computationally feasible at
all.

**It fell short of its target.** The run aimed at `0.065`, just below the
numerical onset, and stopped at `0.05986` on its wall-clock budget. So the gap
`(0.05986, 0.09)` between certified non-detection and certified detection is
**open**, and the numerical onset bracket `(0.0675, 0.0685)` lies **inside**
that gap, uncertified from both sides. That bracket is a midpoint diagnostic,
**not a certified sharp threshold** — no monotonicity and no uniqueness of the
crossing is claimed.

### New analytic input — a sharper archimedean budget
[`ANALYTIC_TAIL.md`](ANALYTIC_TAIL.md)

At the **same** truncation `M = 2^17 - 1`, by two and then three integrations by
parts of an explicit finite trigonometric integral, plus the endpoint identities
`C(U) = 0`, `C'(0) = 0`, `C'(U) = -C'(0)`, `C''(U) = 0`,
`C''(0) = -||f'||^2 = -sum_k w_k^2 x_k^2`:

| | published `O(1/M)` | Prop. T `O(1/M^2)` | Prop. T3 `O(1/M^3)` |
|---|---|---|---|
| `Q_zeta` enclosure width | `9.870e-05` | `3.339e-09` | `2.477e-12` |
| `gamma_E` term `T` width | `3.881e-11` | `3.881e-11` | `1.357e-16` |

**Each enclosure is contained in the previous one** — three valid enclosures of
the same real number, each inside the last. This is **not** a precision
increase: there is no floating-point arithmetic in this code and not one extra
term is summed. The `O(1/M^2)` budget is what makes Theorem A possible; the
`O(1/M^3)` budget is what makes Theorem O1 possible, since the quantity
certified there is `1.75e-10`.

## A correction carried out here

An earlier scan of ours computed the detection diagnostic `Theta(delta)` from
`G_core`, the **tail-free** matrix, rather than from `G`. Since
`lambda_min(G_core) = 5.722182e-06` while `lambda_min(G) = 1.746061e-10` — the
omitted diagonal term `(3/4)S_1 = 5.7221e-06` nearly cancels it — **every
`Theta` value and the onset bracket `(0.085, 0.09)` from that scan are
withdrawn.** The corrected table is in
[`NON_DETECTION.md`](NON_DETECTION.md) §3; the diagnosis is in
[`FOUR_MODE_POSITIVITY.md`](FOUR_MODE_POSITIVITY.md) §4.1.

**The detection certificate never depended on it** — it is a direct sign
evaluation with the tail included and uses no matrix — and the error was in the
safe direction: the true onset is *lower* than was reported, so the certified
detection interval stands further from it, not closer.

## Independent corroboration, and its limit

The repository's separately coded digamma/trigamma implementation, run on the
**new** witness, reaches a negative upper bound of its own at `delta = 1/5` and
`delta = 1/4`, and its `Q_zeta` enclosure **contains** the primary `O(1/M^3)`
one. It is a **second implementation, not external peer review**: the frozen
input, the explicit formula, the canonical `(A_k,B_k)` reduction and the idea of
rational interval arithmetic are shared.

Its budget is ~`1.1e-04` here, roughly `90x` the worst margin of Theorem A, so
**it does not reach the new left endpoint `delta = 9/100`, which rests on the
primary implementation alone.**

## Files

| | |
|---|---|
| Theorem A | [`CERTIFICATE_THEOREM.md`](CERTIFICATE_THEOREM.md) · [`verify_new_certificate.py`](verify_new_certificate.py) · [`results/certificate_new.json`](results/certificate_new.json) |
| Theorem O1 | [`FOUR_MODE_POSITIVITY.md`](FOUR_MODE_POSITIVITY.md) · [`verify_route_o.py`](verify_route_o.py) · [`results/route_o_step1.json`](results/route_o_step1.json) |
| Theorem O2, onset diagnostic | [`NON_DETECTION.md`](NON_DETECTION.md) · [`results/route_o_step2.json`](results/route_o_step2.json) · [`results/theta_table.json`](results/theta_table.json) |
| New analysis, with proofs | [`ANALYTIC_TAIL.md`](ANALYTIC_TAIL.md) |
| Evidence tiering, duplication check, limitations | [`METHODS_AND_LIMITATIONS.md`](METHODS_AND_LIMITATIONS.md) |
| Inputs, hashes, error budget, environment | [`PROVENANCE.md`](PROVENANCE.md) |
| How to run it | [`REPRODUCE.md`](REPRODUCE.md) |
| Second implementation | [`scripts/indep_crosscheck.py`](scripts/indep_crosscheck.py) · [`results/indep_crosscheck.json`](results/indep_crosscheck.json) |

## Quick start (standard library only — no dependencies)

```bash
git clone https://github.com/innerlightr-wq/weil-window-control.git
cd weil-window-control/extensions/four-mode-positivity-2026-10-14

python3 tests/test_core_vs_published.py   # reproduces the PUBLISHED baseline  ~6 min (first run)
python3 tests/test_tail3.py               # Propositions T3, Lemmas E, S1, T    ~1 min
python3 verify_new_certificate.py         # Theorem A                          ~17 min
python3 verify_route_o.py                 # Theorems O1 and O2                 ~45 min
python3 scripts/indep_crosscheck.py       # second implementation              ~8 min
```

Each exits `0` on success and **non-zero** if any gate fails; the gates assert
the inequalities themselves, they do not print a stored label. No `pip install`,
no virtual environment — `mpmath`, `numpy`, `sympy` and `scipy` are **not**
used and were not present on the machine that produced these results. Every
script resolves its inputs relative to its own location and **verifies the
sha256 of the published files it imports**, so a fresh checkout works unchanged
and cannot silently run against a modified copy. See
[`REPRODUCE.md`](REPRODUCE.md).

## Non-claims

1. **No proof, weakening, or conditional version of RH.** See the scope box.
2. **The quartet is fabricated.** `gamma = 14` is prescribed, not observed. No
   statement about any actual zero of `zeta`.
3. **Theorem O1 is finite-dimensional.** It is positivity of `G` on one explicit
   four-dimensional subspace, **not** of the Weil form, and its sign is the one
   RH predicts — a negative value would have been the informative outcome.
4. **Theorem O2 is non-detection in a finite basis, not positivity.**
5. **No threshold, monotonicity, uniqueness or sharpness.** `9/100` and
   `0.05986` are where the constructions stopped. The onset bracket is a
   diagnostic.
6. **Not uniform in `L`, `gamma` or dimension.** One window, one height, four
   modes.
7. **Computer-assisted, not formally verified.** Exact rational interval
   arithmetic run by CPython; no Lean, nothing formalised, no proof assistant.
8. **No novelty** for Proposition 1 (rank-one/Rayleigh threshold),
   Sherman–Morrison, the Schur complement, Sylvester's criterion and law of
   inertia, exact rational `LDL`, Gershgorin/row-sum bounds, polarization,
   AM–GM, Hermite–Hadamard, integration by parts, or interval arithmetic. The
   baseline-weighted matrix criterion in particular is **reused, not renamed** —
   it already existed in the author's own earlier notes and is classical anyway.
   See [`METHODS_AND_LIMITATIONS.md`](METHODS_AND_LIMITATIONS.md) §2.
