# RESULT.md

## Verdict: **PARTIAL**

The known positive construction has been connected **exactly** to the repository's `Q_ζ`, and one
sharply stated unsolved estimate remains. That estimate is **not** a restatement of RH: it is a
specific finite-window extension whose blocking step is named in the sources themselves.

> **CORRECTION 2026-10-07 — see [`SEMILOCAL_LITERATURE_UPDATE.md`](SEMILOCAL_LITERATURE_UPDATE.md).**
> This audit's statement that the semi-local setting "is not available" was **too strong**.
> Semi-local Sonin spaces and the relevant structural identifications are established in later work
> (Connes–Consani–Moscovici, arXiv:2310.18423v2, 2024). The present audit did not locate the
> specific positive-trace comparison with sufficient discrepancy control needed for our
> application. **This is a gap in the verified transfer, not an absence of semi-local geometry.**
> The verdict **PARTIAL** and every statement below about the archimedean place are unaffected.

---

## 1. Strongest statement established (all inputs are Connes–Consani's)

For `g ∈ C_c^∞(R*_+)` supported in the **open** interval `(2^{-1/2}, 2^{1/2})`, with
`ĝ(i/2) = ĝ(−i/2) = 0` and `ĝ(0) = 0`:

```
      Q_zeta(g)  =  P_F(g)  +  R_F(g)
      P_F(g)  =  Pole(g) + Tr( theta(g) S theta(g)* )   >= 0     (independently constructed)
      R_F(g)  =  - E(g*g*)  -  sum_p W_p(g*g*)          >= 0     on this test space
```

because (i) `Tr(ϑ(g)Sϑ(g)*) = ‖ϑ(g)S‖²_HS ≥ 0` is a Gram square requiring **no** `ϑ`-invariance of
the Sonin space; (ii) `Pole(g) = 2F(i/2)² ≥ 0` in the even sector; (iii) `E` is identified
independently by **Theorem 4.6** and bounded by **Lemma 6.10 / Theorem 6.11**, giving
`E ≤ c|ĝ(0)|² = 0` with `c = 4γ/log2`; and (iv) `Σ_p W_p(g*g*) = 0` because **no prime power is
active** at that window.

**My contribution is (iv) and the placement, not the inequality.** The dictionary — that their
`I = [−½log2, ½log2]` is our `[−L,L]` with `L = ½log2`, that `W_∞` is our `Arch − log π` block and
*not* the whole form, and that `e^{2L} = 2` puts `p = 2` exactly on the boundary so the prime block
is empty — is derived and checked here. Everything else is theirs.

## 2. The precise gap

> **Open:** a positive trace with an independently identified, controlled discrepancy for the
> semi-local Weil distribution at `{∞, 2}` — or, equivalently, an extension of Lemma 6.10 from
> `I = [−½log2, ½log2]` to any strictly larger interval.

Both candidate routes fail at the source, not at an obstruction found here:

- **Enlarging the window** breaks *two* conditions simultaneously: Lemma 6.10 is proved on that
  specific `I`, and `W_2` switches on as soon as `L > ½log2`. These are separate conditions and
  both fail at the same point.
- **The semi-local setting** ~~is not available~~ **does not supply the needed estimate**
  *(corrected 2026-10-07)*. arXiv:2008.10974 proves `u = ρ_∞∏ρ_p` is quasi-inner and that the
  semi-local Sonin spaces are infinite-dimensional and filtering, but it *defines* them as
  `ker U₂₂`, calls itself "a **first test** pertaining to the general **strategy**", and
  **postpones to a future paper** the proof that the definition reproduces the semi-local analogue.
  **That postponed work has since appeared:** CCM arXiv:2310.18423v2 (2024) defines
  `S_λ(X_S,α)` directly (Def. 4.5), maps it from the archimedean space compatibly with Fourier
  (Prop. 4.6–4.7), and proves stability in `S` (its Theorem 4.6 / Theorem 2). What is still absent
  from the checked sources is the semi-local analogue of **CC2021 Theorem 4.6 and Theorem 6.11** —
  an independently defined discrepancy `E_S` plus a bound on it; CCM2024 calls that a *"more
  precise strategy"* it *"expect[s]"* to work. Note also that CCM2024's `θ_S` is a **hilbertian**
  isomorphism (their footnote 2: the topological vector space structure), i.e. bounded and
  invertible but **not** isometric, so archimedean orthogonality does not transport. Details and
  claim-by-claim status: [`SEMILOCAL_LITERATURE_UPDATE.md`](SEMILOCAL_LITERATURE_UPDATE.md).

Accordingly **no semi-local `P_F` is written down here.**

## 3. Minimal checks actually run (`checks.py`, exact rational unless stated)

| § | check | result |
|---|---|---|
| 1 | their `I` = our `[−L,L]`, `L = ½log2 = 0.34657359…` | PASS |
| 1 | prime activation vs `L`: none at `½log2` (`e^{2L} = 2` exactly, boundary); `{2,3}` at `L = 0.55`; `{2,3,4}` at `L = 0.8`; `{2,3,4,5,7}` at `L = 1.0` | PASS |
| 2 | **control**: `A S A*` is PSD (all leading minors `≥ 0`) for an explicit `A` that does **not** commute with `S`, because `A S A* = (AS)(AS)*` | PASS |
| 3 | Lemma 6.9 criterion `≡` (trace `≥ 0` and det `≥ 0`) on five cases via the 2×2 reduction, **including deliberately failing ones** — with `(a,b,c) = (3,1,1)`, alignment `1` passes and alignment `1/10` fails | PASS |
| 4 | `c = 4γ/log2 = 16.98658…` recomputed from the quoted `γ`; the source prints no decimal, so only the arithmetic is asserted | PASS |

No large sweeps, no new constant searches, no background jobs, no packages installed. No Cholesky
factorisation of a matrix whose positivity is the question. Numerical agreement is not treated as
certification anywhere.

## 4. Novelty

**None claimed.** The positive Sonin trace and the alignment criterion are prior art
(Connes–Consani, July 2021 PDF, Theorems 1, 4.6, 6.11 and Lemma 6.9). The dictionary and the
prime-activation computation are elementary bookkeeping. The Ramanujan comparison transfers
nothing — that repository's own README already states the shared lattice "is not a mechanism" —
and it is retired after the mechanism table.

## 5. Does our projection work contribute? Diagnostically only.

Stated plainly because the brief asks for it: **no positivity estimate in this audit comes from our
sensitivity law.** Our Theorem 7 is a *relative* budget, not a baseline; Theorem 8 stays
conditional; `K(L,γ)` remains a norm and not an attained coefficient; our `E = ker A_0` is a finite
exact kernel and the Sonin `S` is an infinite-rank projection, with no map between them and no
transfer of the quartic formula. What our work supplies is the *vocabulary* — magnitude, alignment,
response — that makes Lemma 6.9 legible, which is interpretation, not estimate.

One genuinely useful structural observation does fall out: the Sonin compression is a **third** kind
of projection, distinct from both already catalogued in
[`../second-moment-generalization-2026-10-07/PROJECTION_CONTEXT.md`](../second-moment-generalization-2026-10-07/PROJECTION_CONTEXT.md).
It is neither an exact reducing projection (commutes with the operator) nor a baseline-kernel
spectral projection (fixes a leading coefficient): it is a **Gram/compression** projection, positive
by construction and carrying **no invariance at all**. That is why `ran S` need not be
`ϑ`-invariant, and the finite-dimensional control in `checks.py` §2 shows exactly this.

## 6. Comparison with what the project already has

The intrinsic construction gives archimedean-only positivity at `L = ½log2 ≈ 0.3466`. The
manuscript already quotes Zhu's certified `λ*(0.8) ≥ 8.9×10^{-18}` at `L ≤ 0.8`, a **larger** window
that **does** contain active primes. These are different kinds of statement — a structural reason
and operator inequality versus a numerical certificate — and neither subsumes the other. **But as a
positivity range the intrinsic route does not currently reach as far as the certificate already in
use.** No improvement is claimed, and none was found.

---

## Plain language

Connes and Consani build a quantity that is positive for free: it is a sum of squares in a Hilbert
space, so it cannot be negative no matter what the primes do. They then prove an exact identity
saying how far that quantity sits from the archimedean part of the Weil form, and they bound the
difference. The result is genuine positivity with a reason behind it, rather than positivity
observed numerically.

Two things stop this from touching the arithmetic. First, their window is small — in our
coordinates, half-length `½log 2 ≈ 0.347` — and at exactly that size the prime `2` has not yet
switched on, so the theorem is about the archimedean place alone. Second, the version of the
construction that would include primes does not reach far enough yet. *(Corrected 2026-10-07: this
paragraph previously said the semi-local construction "does not exist yet". The spaces do exist —
the 2020 paper's postponed geometry was carried out in CCM arXiv:2310.18423v2, 2024, which defines
the semi-local Sonin spaces outright and proves they are stable as you add primes. What is missing
is the estimate: the map joining the two pictures is only an isomorphism of topological vector
spaces, not an isometry, so the archimedean "positive for free" does not carry across, and the
authors describe the semi-local positivity comparison as a strategy they expect to work rather than
a theorem.)*

Our own projection work does not help with the estimate. It supplies the language for reading their
alignment lemma — how big the perturbation is matters less than whether it points along the
directions that carry the answer — and it identifies their construction as a third kind of
projection we had not catalogued. That is worth recording, and it is not progress on RH.

**What independently supplied identity makes the positive structure constrain the actual arithmetic
Weil form, rather than merely resemble it?** Theorem 4.6: `Tr(ϑ(f)S) = W_∞(f) + E(f)`, with `E`
defined independently of `W_∞` rather than by subtraction, and then bounded above by a rank-one
operator in Lemma 6.10. That is the real bridge. It currently reaches one place, the archimedean
one, on one window. *(Corrected 2026-10-07: the semi-local **spaces** this would need are no
longer postponed — see [`SEMILOCAL_LITERATURE_UPDATE.md`](SEMILOCAL_LITERATURE_UPDATE.md). The
semi-local **identity and bound** are still only a proposed strategy.)*
