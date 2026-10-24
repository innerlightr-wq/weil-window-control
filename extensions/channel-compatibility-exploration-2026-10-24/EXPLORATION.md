# Background-normalized channel interaction on the Davenport–Heilbronn window form

A bounded exploration. Nothing here is committed, published, or offered as a certificate,
and nothing outside this directory was modified.

**Repository state at the start.** `weil-window-control` on `main` at
`e2cf7153ce6903be7c8560ef26c20ea449babbd9`, working tree carrying three pre-existing
untracked paths (`extensions/newtonian-weil-response-2026-10-08/` and two
`paper/*.fdb_latexmk`/`*.fls` build products), all left untouched. `ou-threshold-lean` on
branch `generalize-palindromic-determinant` at
`ff5aebc92ad5bef0e9beb70124621a2b872385db`, clean, read only — not built, and no OU
quantity is imported.

---

## Verdict, in plain language

**Reformulation only.** The channel picture reproduces the existing eigenvalue signs
exactly, and every quantity it produces turns out to be something already available or
already recorded:

* The "compatibility condition" `η ≤ 1` is **algebraically identical** to `det Q ≥ 0`.
  I verified `η = 1 − det Q/det B` to between `10⁻⁶⁷` and `10⁻⁵¹` relative at six
  configurations. Since a rank-one channel can create at most one negative eigenvalue,
  "`Q ⪰ 0 ⟺ det Q ≥ 0`" holds, and the threshold `1` is a normalization of a determinant
  sign test. It is not an independently derived constant.
* The one candidate *mechanism* — that negativity arises because the negative channel
  aligns with weak directions of the background — is **false, and measurably so**. At
  `L = 2.0` the near-null direction of `B` (eigenvalue `2.05×10⁻²⁴`) contributes
  `3.34×10⁻⁶` of `η = 5.602`, i.e. six millionths of one percent. The channel barely
  couples to it.
* What the channel *does* align with is `Q`'s own negative eigenvector, at `0.951`. So
  the reduction says, in different coordinates, that the negative direction of `Q` is the
  negative direction of `Q`.
* The one quantitative structural fact it surfaces — that the channel is orthogonal to
  the remaining-zero form's near-null direction to `9.0×10⁻¹²` — is **numerically the
  same measurement the existing work already reports**: the parent extension's Phase G
  records an overlap of `9.07×10⁻¹²` between that near-null direction and
  `u cos(γ₁u)`, and I measure the channel `s` to be that functional to within
  `1.3×10⁻⁴`.
* There is **no computational advantage**. Deciding positivity by Cholesky on `Q`
  directly took `0.02 s`; the channel route (Cholesky on `B`, then one solve) took
  `0.15 s`; full eigendecomposition took `0.27 s`. The channel route is cheaper than
  eigendecomposition and seven times more expensive than the test it would replace.

The methodological analogy to the OU threshold does **not** survive contact with the
numbers, and the reason is specific. In the OU model the compatibility equation is a
Lyapunov equation whose solution is an independently computable closed form, so the
threshold `3 + 2√2` is an *output*: you can evaluate the criterion without solving the
original existence problem. Here `η = u^T B^{-1} u` requires `B`, which requires `Q`,
which is the whole original object. A criterion that needs the full arithmetic-side
matrix is not an independent explanation of that matrix's sign.

**One thing worth keeping** is the exact nested-enrichment increment
`K_W − K_V = Eᵀ S⁻¹ E`, which localizes precisely what new basis modes add after their
overlap with the old space is projected out. Its nonnegativity is automatic; its value is
that `E` and `S` are separately inspectable. That is a diagnostic, not a lead, unless one
of them can be bounded without solving the problem — see §4.

---

## 1. What was tested and reproduced

### 1.1 The phase removal, derived rather than assumed

The parent extension records the quartet identity
`Q_O(g) = 4(|P|² − |S|²)` with `P = ∫g cosh(δu)e^{iγu}du`,
`S = ∫g sinh(δu)e^{iγu}du`, for the complete orbit
`O = {β ± iγ, 1−β ± iγ}` and `δ = β − ½`. Neither functional is real in general. In a
**single real parity sector** they become real one at a time, and which one is which
swaps with the sector:

* `g` **odd** ⟹ `g·cosh(δu)` is odd and `g·sinh(δu)` is even, so
  `P = i∫g cosh(δu) sin(γu) du` is purely **imaginary** and
  `S = ∫g sinh(δu) cos(γu) du` is purely **real**;
* `g` **even** ⟹ the roles swap: `P` is real, `S` purely imaginary.

Hence in each sector `|P|² = p(g)²` and `|S|² = s(g)²` for two **real** linear
functionals, and `Q_O = 4pp^T − 4ss^T` with real vectors. The phases cancel; they are not
assumed away. `channel_vectors` computes `P = (A+B)/2`, `S = (A−B)/2` from
`A = G(δ+iγ)`, `B = G(−δ+iγ)` and **checks the component that must vanish**. In the odd
sector it vanishes **exactly** — `0.0`, not merely to the working precision — because the
closed form for `G` carries `sin(w L) = 0` identically for the odd basis
(`wL = (k+1)π`). The rank-two matrix was separately checked against the parent's scalar
`quartet_contribution` by polarization: residual `1.6×10⁻⁸³` at 80 digits.

All of this work is in **one sector at a time**; the positive/negative channel split is
*not* the even/odd sector split, and is *not* the prime-power/mixed-prime split of
Phase H.

### 1.2 Basis, nesting, Gram matrix

The declared level-2 and level-3 index sets of the parent extension were read, not
assumed. At each fixed `L` the level-2 set is a **subset** of the level-3 set, with
exactly 20 new modes:

| `L` | `k₀` | level 2 | level 3 | `V ⊂ W`? | new |
|---|---|---|---|---|---|
| 1.6 | 43 | 37 modes | 57 modes | yes | 20 |
| 1.8 | 48 | 37 | 57 | yes | 20 |
| 2.0 | 54 | 37 | 57 | yes | 20 |

Across different `L` the levels are **not** comparable, since the basis frequencies
`w_k = (k+1)π/L` change; no cross-`L` nesting is used or implied anywhere.

The Gram matrix of the normalized basis is computed (through the same cross-correlation
machinery the form itself uses) and found to be the identity to `5.5×10⁻⁸⁰`, so the
reported eigenvalues are already the generalized ones and no direction was discarded or
regularized. (An earlier version of this check compared `norm2(k)/norm2(k)` against `1`,
which is a tautology that cannot fail; it was replaced.)

### 1.3 The frozen witnesses reproduce

From their **stored** coefficients, re-parsed at 80 digits:

| witness | `L` | `N` | `Q_f(g)/‖g‖²` | sign |
|---|---|---|---|---|
| `witness.json` | 2.0 | 37 | `−0.679071346441480590986377586093` | negative |
| `witness_smallest_L.json` | 1.8 | 57 | `−1.68123895523670103326233298941×10⁻¹⁶` | negative |

Both match the parent extension's certified values to all stored digits. These are
benchmark checks; nothing was refitted.

### 1.4 The reduction and its gates

With `Q` the complete arithmetic-side matrix (the full form — the remaining spectrum is
**never** replaced by a finite zero list, and the Phase H support deletion is **not**
repeated), `R = Q − Q_O`, `u = 2s`, and `B = R + 4pp^T = Q + uu^T`:

| `L` | lvl | `N` | `λ_min(Q)` | `#neg Q` | `λ_min(B)` | `B ≻ 0` | `η` | `η ≤ 1` | `Q ⪰ 0` |
|---|---|---|---|---|---|---|---|---|---|
| 1.6 | 2 | 37 | `+2.15385166477×10⁻¹⁶` | 0 | `2.22011119471×10⁻¹⁶` | yes | `0.36560043476094` | yes | yes |
| 1.6 | 3 | 57 | `+1.09593824675×10⁻¹⁸` | 0 | `1.10349747349×10⁻¹⁸` | yes | `0.39016683249573` | yes | yes |
| 1.8 | 2 | 37 | `+1.00754195149×10⁻²⁰` | 0 | `1.25670510595×10⁻²⁰` | yes | `0.85006833019136` | yes | yes |
| 1.8 | 3 | 57 | `−1.68123895524×10⁻¹⁶` | 1 | `1.03630169384×10⁻²⁵` | yes | `1.1236970969958` | no | no |
| 2.0 | 2 | 37 | `−0.679071346441` | 1 | `2.05195003117×10⁻²⁴` | yes | `5.6022332216115` | no | no |
| 2.0 | 3 | 57 | `−0.71017551223` | 1 | `4.14691478×10⁻³²` | yes | `146.11354446381` | no | no |

`B` is positive definite at **all six** configurations — which is not a given, and is
inherited from the parent extension's Phase E measurement that `λ_min(R) > 0`, since
`B = R + 4pp^T ⪰ R`. `Q` has **exactly one** negative direction wherever it has any, so a
rank-one channel is capable in principle of accounting for it; had `#neg Q ≥ 2`, `B`
could not have been positive definite and the single-channel reduction would have been
inapplicable.

Validation gates, all passing: the criterion `η ≤ 1 ⟺ Q ⪰ 0` agrees with the direct
eigenvalue sign at 6/6; the inertia relation (number of `K`-eigenvalues above `1` equals
`#neg Q`) agrees at 6/6; reconstruction `‖Q − (B − uu^T)‖_max ≤ 4.2×10⁻⁸¹`; the identity
`‖B − (R + 4pp^T)‖_max ≤ 8.4×10⁻⁸¹`; linear-solve residuals `~10⁻⁷³` at 80 digits. **These
are gates, not findings.**

### 1.5 Precision ladder

At the crossing configuration `L = 1.8`, level 3 — the delicate one, with
`cond(B) = 8.45×10²⁵`:

| digits | `η` | solve residual | determinant-identity residual |
|---|---|---|---|
| 80 | `1.1236970969958460727` | `1.04×10⁻⁷³` | `3.02×10⁻⁶⁰` |
| 140 | `1.1236970969958460727` | `1.02×10⁻¹³³` | `3.79×10⁻¹²⁰` |
| 220 | `1.1236970969958460727` | `1.15×10⁻²¹³` | `1.30×10⁻¹⁹⁹` |

Identical to all printed digits, with residuals tracking the precision. The gap from the
threshold, `η − 1 = 0.1237`, is resolved by far more than the conditioning costs. **This
is stability under increased precision, not an interval certificate**, and is nowhere
described as one.

### 1.6 Coordinate invariance

`η` must not depend on the coordinates chosen on a fixed test space. Under
`c = Tc'` the form becomes `T^TQT`, the channel `T^Tu`, and
`η' = u^TT(T^TBT)^{-1}T^Tu = u^TB^{-1}u = η`. Checked numerically with a fixed,
deliberately badly scaled invertible `T` (diagonal entries spanning six orders of
magnitude plus off-diagonal coupling; no randomness), carrying the Gram matrix
`G → T^TGT` through the generalized eigenproblem.

The Gram deviation is `8.2×10⁻⁸⁰` at 80 digits, `1.2×10⁻¹³⁹` at 140, `1.1×10⁻²¹⁹` at 220
— the basis is orthonormal and the deviation is a precision floor. `η` is unchanged to
all printed digits under the change of coordinates (relative difference `9.46×10⁻⁶⁵`),
and the generalized `λ_min(Q', G')` reproduces the orthonormal `λ_min(Q)` to
`2.0×10⁻⁸⁰`.

---

## 2. What is automatic linear algebra, and what adds information

### Automatic — carries no information beyond the eigenvalue computation

1. **`Q ⪰ 0 ⟺ η ≤ 1`.** The Schur-complement criterion for a rank-one downdate. Built in.
2. **`#neg Q = #{eigenvalues of K above 1}`.** The inertia count. Built in.
3. **`η` increases under enrichment.** `K_W − K_V = EᵀS⁻¹E ⪰ 0` whenever `B_V, S ≻ 0`.
   Automatic from the identity, as the brief says.
4. **`η = 1 − det Q / det B`** — the matrix-determinant lemma. **This is the one that
   settles the question**, and I verified it rather than assuming it:

   | `L` | lvl | `η` | `1 − det Q/det B` | relative residual |
   |---|---|---|---|---|
   | 1.6 | 2 | `0.36560043476094143427` | `0.36560043476094143427` | `1.87×10⁻⁶⁷` |
   | 1.6 | 3 | `0.39016683249572650823` | `0.39016683249572650823` | `3.54×10⁻⁶⁵` |
   | 1.8 | 2 | `0.85006833019135937512` | `0.85006833019135937512` | `2.16×10⁻⁶⁵` |
   | 1.8 | 3 | `1.1236970969958460727` | `1.1236970969958460727` | `3.02×10⁻⁶⁰` |
   | 2.0 | 2 | `5.6022332216114890407` | `5.6022332216114890407` | `3.99×10⁻⁵⁹` |
   | 2.0 | 3 | `146.11354446381254239` | `146.11354446381254239` | `1.51×10⁻⁵¹` |

   So the "compatibility condition" `η ≤ 1` **is** the sign test `det Q ≥ 0`. The
   threshold `1` is not an independently derived constant; it is where a determinant
   changes sign, rewritten. This is the decisive disanalogy with the OU threshold, where
   `3 + 2√2` is an algebraic number that falls out of a Lyapunov solution and can be
   evaluated **without** solving the existence problem it decides.

5. **The suppression of channel overlaps along weak background directions is automatic.**
   This is worth spelling out because it looks like a mechanism and is not. For any unit
   eigenvector `v_j` of `B`, the term `|v_jᵀu|²/β_j` is one summand of `η`, so
   `|v_jᵀu|² ≤ η·β_j` *identically*. Overlaps along small `β_j` are therefore forced to
   be small in proportion, for no reason other than that `η` is finite — and `η` is
   finite precisely because `B ≻ 0`, which is inherited from the parent extension's
   Phase E measurement `λ_min(R) > 0`. The observed pattern
   `|v_jᵀu|²/β_j ∈ {0.851, 0.744, 0.409, …}` at `L = 2.0` is this inequality being
   roughly saturated, not an arithmetic coincidence.

### Adds information — three measurements, all negative or already known

6. **The near-null hypothesis is false, quantitatively.** The brief asks whether
   negativity arises because the channel aligns with weak directions of the background.
   It does not. At `L = 2.0` level 2, the smallest background eigenvalue
   `β = 2.05×10⁻²⁴` contributes `3.34×10⁻⁶` to `η = 5.602` — six millionths of one
   percent — because the channel's overlap with it is `|v̂ᵀû| = 2.51×10⁻¹⁵`. The brief's
   own caution applies exactly: a tiny background eigenvalue is not evidence of strong
   interaction when the channel barely couples to it.

7. **`η` is dominated by ordinary directions.** Grouping clustered eigenvalues and
   sorting by contribution:

   | configuration | `η` | largest single contribution | its `β` | top 3 | listed 8 |
   |---|---|---|---|---|---|
   | `L=1.6` lvl 2 | `0.3656` | `41.4%` | `1.5403` | `78.3%` | `92.9%` |
   | `L=1.8` lvl 3 | `1.1237` | `25.3%` | `0.9822` | `54.5%` | `93.8%` |
   | `L=2.0` lvl 2 | `5.6022` | `58.9%` | `0.2729` | `87.4%` | `99.1%` |

   In every case the dominant direction has `β` of order one. Summing, among the listed
   groups, every direction with `β < 10⁻⁶` gives `6.96%` (`L=1.6` lvl 2), `27.73%`
   (`L=1.8` lvl 3) and `36.76%` (`L=2.0` lvl 2) of `η`: substantial, never dominant, and
   never the explanation. So
   "a small spectral subspace accounts for most of the normalized interaction" is true in
   the weak sense that one to three directions carry most of it — but they are the
   *ordinary* ones, which is the opposite of the conjectured mechanism.

8. **What the channel actually aligns with is already recorded.** At `L = 2.0` level 2:

   | overlap | value |
   |---|---|
   | `\|⟨û, v₁(Q)⟩\|` — `Q`'s own **negative** eigenvector | `0.95090319` |
   | `\|⟨û, v₂(Q)⟩\|` — `Q`'s near-null eigenvector | `5.45×10⁻¹⁶` |
   | `\|⟨û, v_min(B)⟩\|` | `2.51×10⁻¹⁵` |
   | `\|⟨û, v_min(R)⟩\|` — the remaining form's near-null direction | `8.99×10⁻¹²` |
   | `\|⟨ŝ, t̂⟩\|`, `t_k = ∫f_k(u)·u cos(γ₁u) du` | `0.9998676112` |
   | `\|⟨t̂, v_min(R)⟩\|` | `9.08×10⁻¹²` |

   The last two lines are the point. The negative channel `s` **is** the
   derivative-evaluation extremiser direction `u cos(γ₁u)` to `1.3×10⁻⁴`, and that
   direction's orthogonality to the remaining form's near-null vector is
   `9.08×10⁻¹²` — which is the parent extension's Phase G reported overlap,
   `9.07×10⁻¹²`, recovered to three digits. The channel picture's one quantitative
   structural fact is a measurement the existing work already published. Likewise
   `|⟨û, v₁(Q)⟩| = 0.9509` sits beside Phase G's `0.9514` overlap of the frozen witness
   with `u cos(γ₁u)`: the same alignment, seen from a different side.

---

## 3. The strongest supported finding

> At `L = 1.8`, the step from the 37-dimensional level-2 space to the nested
> 57-dimensional level-3 space raises the background-normalized interaction from
> `K_V = 0.85006833019135937512` to `K_W = 1.1236970969958460727`, i.e. **across the
> threshold `1`**, and the increment is accounted for exactly by the Schur-complement
> identity
> `K_W − K_V = EᵀS⁻¹E = 0.27362876680448669759`,
> with `S = D − CᵀB_V⁻¹C` and `E = u_new − CᵀB_V⁻¹u_V`, to a relative residual of
> `4.3×10⁻⁶⁵` at 80 digits and `6.0×10⁻¹²⁵` at 140 digits.

Supporting consistency checks: the `V`-block of `B_W` equals the background built
intrinsically on `V` **exactly** (residual `0.0`), as does `u_V` (residual `0.0`), so the
partition is the right one. `‖E‖ = 0.3909046186`, and `S` has
`λ_min = 6.68×10⁻¹¹`, `λ_max = 6.186153548`.

Qualifications, all of which matter:

* **Finite space.** This is one nested pair at one window in one parity sector. It says
  nothing about any limit, and the parent extension's positive level-2 value at `L = 1.8`
  (`λ_min(Q) = +1.01×10⁻²⁰`) must **not** be read as positivity of the window — the
  negative witness lives in the larger space, and that is exactly what the increment
  above expresses.
* **Precision, not certification.** All values are stable at 80/140/220 digits; none is
  an interval enclosure.
* **The quartet is a numerical input.** `β₁`, `γ₁` come from the parent extension's
  census and are T2, not independently certified. Every statement here that involves the
  selected true zero orbit inherits that qualification.
* **Removing `O₁` does not leave only on-line zeros.** `R = Q − Q_O₁` still contains the
  other off-line quartets implicitly, since `Q` is the complete arithmetic form.
* **Nonnegativity of the increment is automatic**, as the brief states. What is reported
  here is its *size*, and that it is the quantity responsible for the crossing.

---

## 4. The principal obstruction

**The increment cannot be bounded from the quantities that describe it.** The identity
gives `K_W − K_V = EᵀS⁻¹E`, and the two bounds available from `‖E‖` and the extreme
eigenvalues of `S` are

```
‖E‖²/λ_max(S) = 0.1528/6.186     = 0.0247     (lower)
‖E‖²/λ_min(S) = 0.1528/6.68e-11  = 2.3×10⁹    (upper)
```

against a true value of `0.2736`. The bracket spans **eleven orders of magnitude**, and
the truth sits `11×` above the lower bound. So the increment is governed entirely by how
`E` distributes across `S`'s spectrum, and nothing short of that distribution determines
it. Knowing `‖E‖` and `cond(S)` is not enough.

That is the missing estimate: **an independent bound on `EᵀS⁻¹E` in terms of quantities
not derived from `Q` itself.** Without it the formulation cannot predict anything. And
this is where the OU analogy fails structurally rather than numerically: there, the
compatibility equation is solved in closed form, so the criterion is computable from the
model's parameters; here, `η` requires `B`, which requires `Q`, which is the entire
original object. **A criterion that needs the full arithmetic-side matrix is not an
independent explanation of that matrix's sign.**

Two further obstructions, recorded because they bound what a second investigation could
hope for:

* **No computational advantage.** Measured at `N = 37`, 80 digits: Cholesky on `Q`
  directly, which is the standard positivity test, `0.02 s`; the channel route (Cholesky
  on `B`, then one solve) `0.15 s`; full eigendecomposition `0.27 s`. All three return
  the same verdict. The channel route is cheaper than eigendecomposition and about
  seven times **more** expensive than the test it would replace, and it additionally
  requires `B ≻ 0` as a side condition that the direct test does not need.
* **No collective effect between channels.** With `O₁` and `O₂` on a common background at
  `L = 2.0` level 2: `K₁₁ = 5.6014586215603311`, `K₂₂ = 0.00180177135682622124`,
  `K₁₂ = 0.0278065531659852647`, `λ_max(K) = 5.601596698830836` against
  `max(K₁₁,K₂₂) = 5.6014586215603311`. The off-diagonal moves `λ_max` by
  `1.38×10⁻⁴`, i.e. `0.0025%`. `K` is symmetric to `1.03×10⁻⁶⁷` and the inertia relation
  holds (`one` eigenvalue above `1`, `one` negative direction of `Q`). The motivating
  idea — heterogeneous components satisfying a *collective* compatibility condition — is
  not realised here: the two channels are effectively independent, and the criterion
  degenerates to a maximum over individual channels. (`K₁₁ = 5.60146` differs from the
  single-channel `η = 5.60223` because the background is not the same object; the two
  are reported separately and not compared as though they belonged to one system.)

---

## 5. One highest-value next experiment — and whether it is warranted

**It is not warranted as a continuation of this line.** Classified against the brief's
four outcomes, this is **reformulation only**: the formulation reproduces the known
eigenvalue signs and standard rank-update behaviour, the one candidate mechanism is
measurably false, the one structural number it surfaces duplicates an existing Phase G
measurement, and it costs more than the test it would replace. Not *blocked* — every
background was positive definite, the conditioning was handled, and all gates passed; the
calculations are trustworthy. They simply do not add predictive information.

If a second investigation happens anyway, the single highest-value item is **not** in
this formulation but in the obstruction it exposes:

> **Candidate lemma.** Bound `EᵀS⁻¹E` for the nested enrichment in terms of the
> *analytic* data of the new modes — their frequencies `w_k = (k+1)π/L` relative to `γ₁`,
> and the decay of the arithmetic kernel `h(log n)` across them — rather than in terms of
> `S` and `E`, which are themselves assembled from `Q`.

That is the only route by which the increment could become *predictive* instead of
descriptive, and it is deliberately stated as a lemma about the enrichment increment and
**not** as a positivity statement, so it is not equivalent to assuming what one wants to
prove. I would rate its prospects as modest: the eleven-order bracket in §4 shows the
quantity is delicate, and a bound that does not already know how `E` sits inside `S`'s
spectrum will not be sharp enough to decide a crossing.

**What I would not do:** extend to all four quartets, add channels indefinitely, deform
orbit positions or weights to manufacture a collective crossing, or fit a Lyapunov
equation to `Q`. The last in particular would need independent arithmetic justification
that nothing here provides; the OU Lyapunov equation is derived from a stochastic model,
and there is no such derivation on the Weil side.

---

## 6. Scope statements, kept explicit

* **No finite calculation here establishes anything infinite-dimensional, and nothing
  here bears on RH or GRH.** Every matrix is a restriction to a stated finite space.
* **A positive finite restriction is not whole-window positivity.** In particular the
  positive level-2 value at `L = 1.8` is not the window being positive there; the
  negative witness is in the larger space.
* **Removing `O₁` does not leave only on-line zeros.** The other off-line quartets remain
  implicitly in `Q`.
* **The OU analogy is methodological only.** No OU parameter is identified with a zeta
  coordinate, no corridor is transferred, and `ou-threshold-lean` was read, not built.
  The analogy is *how* a compatibility identity yields a sign criterion — and §2 item 4
  and §4 record precisely where it breaks.
* **The channel split is not the sector split and not the Phase H split.** Positive/negative
  channels of one quartet in one parity sector are a different decomposition from
  even/odd and from prime-power/mixed-prime.
* **No novelty is claimed.** The Schur complement, the matrix-determinant lemma, the
  interlacing and inertia relations, and the block-inverse increment are textbook.
* **No manuscript was edited** and nothing was committed, pushed, tagged, or published.
  No wording concerns about the manuscripts arose during this exploration; had they, they
  would be recorded here and nowhere else.

### Corrections made during this exploration

1. **A check that could not fail.** The first version of `run_config` computed
   `gram_dev = max |norm2(k)/norm2(k) − 1|`, which is identically zero. It was replaced
   by a real computation of the Gram matrix through the parent extension's `phaseD.gram`,
   and `src/fix_gram_column.py` recomputed that column for the rows already written.
   This is sound only because `gram_dev` is an independent check that feeds no other
   quantity; a column that fed `η` would have required regenerating the assemblies. The
   recomputed values are `8.2×10⁻⁸⁰` (80 digits) to `1.1×10⁻²¹⁹` (220 digits).
2. **`mp.lu_solve` with a matrix right-hand side.** The invariance routine passed a
   matrix where mpmath expects a vector and raised `IndexError`. Replaced by column-wise
   triangular solves.

### Deviation from the requested file list

The brief asks for five files. I produced those five plus `data/twochannel_results.csv`,
`data/invariance_results.csv`, `src/run_all.sh`, `src/fix_gram_column.py` and the
per-run `data/*.log` transcripts. The two extra CSVs have different schemas from
`channel_results.csv` and deleting them would make §4's two-channel and invariance
results unreproducible; the two extra scripts are the run driver and the correction in
item 1 above. Everything else requested is present, and all ten same-schema channel rows
(the grid and both precision ladders) are consolidated into `data/channel_results.csv`.

---

## Closing answer

> **Does this exploration identify a concrete mechanism or advantage worth a second
> investigation, or does it only reorganize the existing positivity problem?**

**It only reorganizes it.** The background-normalized interaction `η` is exactly
`1 − det Q/det B`; its threshold `1` is a determinant sign test in disguise, not an
independently derived constant. The mechanism it was meant to expose — negativity through
alignment with weak background directions — is false by a factor of `10⁻⁶` at the
configuration where it should have been clearest, and `η` is instead dominated by
ordinary order-one directions. The one quantitative structural fact the picture yields,
the `9.08×10⁻¹²` orthogonality between the negative channel and the remaining form's
near-null direction, is the Phase G measurement `9.07×10⁻¹²` seen from a different angle,
because the negative channel *is* the derivative-evaluation extremiser to four digits.
Two channels do not interact collectively: the off-diagonal moves the criterion by
`0.0025%`. And the reduction is seven times more expensive than a Cholesky factorization
of `Q` itself.

The exploration is nonetheless worth having done, for a reason the brief anticipates: it
converts "maybe the negativity is an alignment effect" from an open suggestion into a
measured negative, and it localizes where a genuinely new estimate would have to live —
in a bound on `EᵀS⁻¹E` from analytic data, not from `S` and `E`. That is a narrower and
more honest target than the one this exploration started with.

**Stopping here.** The proposed lemma of §5 is not executed.
