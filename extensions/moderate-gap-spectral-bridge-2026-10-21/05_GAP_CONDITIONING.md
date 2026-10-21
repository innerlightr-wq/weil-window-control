# 05 — Do moderate gaps improve detector conditioning?

Answers the brief's §6 (eight questions) and Final Question 2.
**Verdict: E — decisive no, on two independent grounds, one combinatorial and one
harmonic-analytic. Neither depends on the other.**

## 1. Ground one — anti-localisation (the gaps are guaranteed exactly where there
is nothing to detect)

> **Proposition 1 (anti-localisation).** `[IV, proved; verify.py J1–J5]`
> Let `j ∈ W_ℓ(T,c)` with `ℓ ≥ 2`. For every interior index `1 ≤ a ≤ ℓ−1`, the
> ordinate `γ_{j+a}` is attained by exactly one zero of `ζ`, which is therefore
> simple and on the critical line. The endpoints `γ_j` and `γ_{j+ℓ}` need not be.
>
> *Proof.* By the functional equation, `ρ` off the critical line forces the distinct
> zero `1 − ρ̄`, and `Im(1−ρ̄) = Im ρ`: **an off-line zero always has a partner at
> ordinate distance exactly 0.** A zero of multiplicity `m ≥ 2` likewise contributes
> `m` equal entries to the ordinate list. For an interior index both adjacent gaps
> are `≥ 2πc/log T > 0`, so no second zero shares `γ_{j+a}`; a zero not sharing its
> ordinate with any other is simple and cannot be off the line. For an endpoint only
> one adjacent gap is constrained, so an endpoint may be the upper (resp. lower)
> member of a zero-gap pair. ∎

Three consequences.

1. **The brief's caution is exactly right.** A run of four `c`-moderate gaps
   certifies **three** simple critical-line zeros, not five. Verified
   combinatorially on 300 random configurations: 0 interior violations, and
   off-line endpoints do occur (`verify.py` J1, J2).
2. **The theorem is structurally blind to off-line zeros.** The set of ordinates it
   certifies as `c`-separated and the set of ordinates of hypothetical off-line
   zeros are **disjoint by construction**, not merely unrelated. No strengthening of
   `c`, `ℓ`, or the proportion can change this: an off-line pair's internal ordinate
   gap is identically zero.
3. **The paper's own proof says the same thing.** ARTZ's Proposition 1 constructs
   runs by *discarding* every string whose interior meets `Z_T^m` (clusters
   containing a simple critical zero and something else) or which straddles a
   cluster of `Z_T^n` (clusters with no simple critical zero) — i.e. by removing
   exactly the neighbourhoods of the multiple and off-line zeros. The output set is
   the complement of the object of interest.

**Answer to the brief's §6 Q7 ("can a hypothetical off-line zero be guaranteed to
lie near one of the favourable blocks?") and Q8 ("does the theorem supply a
guarantee for the actual zero being tested, or only for some favourable zeros?"):
no, and only for some. It is worse than a localisation failure — it is an exclusion.**

## 2. Ground two — the guaranteed separation is sub-Nyquist

The kernel is `K_τ = \hat{η_τ²}` with `η_τ²` supported in `(−1/2,1/2)`. The nodes
`f_x(u) = η_τ(u)e^{−2πiux}` are therefore band-limited exponentials on an interval of
length 1, whose **Nyquist spacing is 1**: at spacing exactly 1 and `η_τ² ≡ 1` the
nodes are orthonormal (`K(j) = 0` for nonzero integer `j`, `verify.py` H5). ARTZ's
theorem delivers separation `c ≤ 1/4` — **four times oversampled.**

> **Proposition 2 (exponential decay of `λ_min` under fixed separation).**
> `[IV, certified; verify.py H1–H4]`
> Let `x_j = jc`, `0 ≤ j ≤ d−1`, `m = d−1`, and let `G_d = [K_τ(x_i − x_j)]`. Taking
> the test vector `a_j = (−1)^j\binom{m}{j}` gives
> `Σ_j a_j e^{−2πiujc} = (1 − e^{−2πiuc})^m`, hence `‖Σ a_j f_{x_j}‖² =
> ∫η_τ²(u)\,(2|\sin(πuc)|)^{2m}du ≤ (2\sin(πc/2))^{2m}` and `‖a‖² = \binom{2m}{m}`, so
> ```
>      lambda_min(G_d)  <=  ( 2 sin(pi c/2) )^{2m} / C(2m, m)
>                        ~  sqrt(pi m) * ( 2 sin^2(pi c/2) )^m .
> ```
> The rate `2\sin²(πc/2) < 1` for every `c < 1/2`; at `c = 1/4` it is
> `2\sin²(π/8) = 0.2928932…`.

Certified upper bounds on `λ_min(G_d)` (exact-rational interval arithmetic):

| `c` | `d=2` | `d=4` | `d=6` | `d=8` | `d=12` | `d=16` | `d=20` |
|---|---|---|---|---|---|---|---|
| `1/4` | `2.9e−1` | `1.0e−2` | `2.7e−4` | `6.9e−6` | `4.0e−9` | `2.1e−12` | `1.1e−15` |

(see `results/verify.txt`, block H, for the full table including `c = 1/3, 1/2, 3/4, 1`).

**So minimum separation at the guaranteed scale makes conditioning *worse*, not
better, as the block lengthens — exponentially so.** This is the exact opposite of
the hoped-for transfer.

### Why this is the relevant regime
In my detector problem the Gram matrix is over **all** ordinates in the window, not
just a run: `d ≈ (window length) × (density) ≍ L·log T/(2π)`. The large-`d` regime is
the only one that occurs. The moderate-gap theorem bounds `d` by `ℓ+1 ≤ 5`, which is
precisely the range in which conditioning is still harmless — so the theorem's own
ceiling and the only safe regime coincide. That coincidence is arithmetic, not
structural, and it is not usable: a 5-point block is far too small to annihilate the
competing ordinates in any finite window (my `INTERP` work needs `d ≳ 3N`; see 06 §4).

## 3. The eight questions of §6, answered

1. **Does minimum separation imply positive definiteness of the finite Gram?** Yes,
   trivially — the `f_x` are linearly independent for distinct `x` (Lamzouri §2,
   from the linear independence of `e^{zu}`). Positive *definiteness* is free and
   worthless. `[II]`
2. **Does it give an effective quantitative lower bound?** For `d = 2`, yes and
   certified: `λ_min = 1 − K_τ(c) ≥ 1 − (Â(c) + 2(1−A_τ)/A_τ) ≥ 1 − 0.9274 = 0.0726`
   at `c = 1/4`, `τ = 1/200`. For larger `d` an explicit lower bound is available but
   it necessarily decays — see 3. `[IV, certified from B3–B6]`
3. **How does the bound depend on `d`?** It must decay at least as fast as
   `√(πm)(2\sin²(πc/2))^m`, by Proposition 2. Empirically (diagnostic, `η²≡1`),
   `λ_min` at `c = 1/4` for `d = 2,4,6,8,10` runs `1.0e−1, 1.6e−3, 1.6e−5, 1.3e−7,
   9.3e−10` — a factor `≈ 10^{−2}` per two extra nodes. `[V]`
4. **How does it depend on support length?** Through the single dimensionless
   product `c × (support length)`. Doubling the support halves the effective
   oversampling factor; the paper's support is fixed at 1, so `c ≤ 1/4` is fixed at
   4× oversampled. `[IV]`
5. **Can conditioning deteriorate despite fixed nearest-neighbour separation?**
   **Yes — provably and exponentially.** This is Proposition 2. `[IV, certified]`
6. **What if the window contains additional zeros outside the selected block?**
   It is strictly worse: adding nodes can only lower `λ_min` (Cauchy interlacing for
   the Gram matrix of a larger family). And in the real problem the window *does*
   contain them — `d ≍ L log T/(2π)` of them. `[I]`
7. **Can a hypothetical off-line zero be guaranteed near a favourable block?**
   **No — it is guaranteed *not* to be in one.** Proposition 1. `[IV, proved]`
8. **Does the theorem guarantee anything for the zero being tested?** No; only for
   a positive proportion of *other*, provably on-line, ordinates. `[IV, proved]`

## 4. Survival under an adversarial surrounding configuration

Tested directly: a single off-line pair at depth `y = 1` surrounded by 24 on-line
zeros placed at the 24 most negative values of `Re[K(Δ − iy)²]`, subject to mutual
separation `≥ 0.30` (`verify.py` block F, case 4). The chain of 04 §1 still holds
(`‖G‖²_F = 3565.91` vs chain RHS `1767.68`), and `trΨ(P)` *rises* from 24 to 32.05 —
the adversarial placement costs the adversary budget rather than saving it. So the
adversary's leverage is **not** in the on-line surroundings. It is entirely in other
*off-line* pairs (04 §3.1), which is a different mechanism and is where the method
actually dies. `[V]`

## 5. What survives as a standalone statement

Proposition 2 is worth keeping on its own: it is an explicit, certified, dimension-
dependent ceiling on how well-conditioned a band-limited Gram matrix can be under a
*fixed* minimum separation below Nyquist. It is not RH-specific and it is the exact
quantitative form of the folklore statement "sub-Nyquist sampling is exponentially
ill-conditioned". **Verdict for Candidate T2: E as an RH transfer; C as a standalone
finite-configuration estimate.**
