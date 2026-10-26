# notes/normalization.md — the three conventions, and the closed forms

## 1. The three conventions side by side

| | Suzuki (arXiv:2610.32704v1) | Zhu / this repository |
|---|---|---|
| Fourier | `û(z) = ∫u(x)e^{izx}dx`, `‖Fu‖² = 2π‖u‖²` [(2.1)] | `F(t) = ∫f(x)e^{itx}dx` (`notes/stage3_assembly.md`) |
| form | `Q_W(u,f) = Σ_{γ∈Γ} m_γ û(γ) conj(f̂(γ))`, linear in the **first** variable [(2.3)] | `Q(f) = Σ_ρ \|F(γ_ρ)\|²` (under RH), `Σ` over nontrivial zeros **with multiplicity** |
| support | `u, f ∈ C_c^∞(−a,a)` | `supp f ⊂ [−L,L]` |
| coercivity | `Q_W(u) ≥ λ_a‖u‖²_{L²(−a,a)}`, `λ_a := inf σ(A_a)` [(2.4)] | `λ*(L) := inf{Q(f)/‖f‖²}` |
| prime cutoff | `p^r ≤ e^{2a}`, equality contributes 0 | `log n < 2L` |

**`Γ` is the set of *distinct* ordinates with multiplicity weights `m_γ`, so
`Σ_{γ∈Γ} m_γ (·)` = `Σ_ρ (·)` over zeros counted with multiplicity.** Both sums run over
`±γ`. The Fourier conventions are **identical**. The support conventions are identical.
Therefore:

> **Q0-a: `a = L`.  Q0-b: `λ_a = λ*(L)` exactly, constant `1`.**

This is an algebraic identification, not a fit. It is verified numerically in P0-T2
against `src/zeta_window.py` — code not written in this extension — and in P0-T4
structurally.

**Caveat carried forward.** Retraction 15 of the parent record notes that on ζ windows the
indefiniteness of the *recorder split* `A = Pole + Arch − LogPi` versus `Prime` depends on
where `−log π` is placed. That affects the **decomposition only**, never the total `Q`.
Nothing here uses the split.

## 2. The geometric side, specialized to the exponential basis

The repository's form (`notes/stage3_assembly.md`), written sesquilinearly:

```
Q(u,f) = h(i/2) + h(−i/2) + (1/2π)∫_ℝ h(t)[Re ψ(1/4+it/2) − log π]dt
         − Σ_{log n < 2L} 2Λ(n)n^{−1/2}·(1/2π)∫_ℝ h(t)cos(t log n)dt ,
         h(t) := û(t) conj(f̂(t)) .
```

`h` continues to the entire function `h(z) = û(z) conj(f̂(z̄))`, because with
`C(v) := ∫u(y+v)conj(f(y))dy` one has `∫C(v)e^{izv}dv = û(z)·conj(f̂(z̄))`. Hence

```
   h(i/2)  = φ_m(i/2)·conj(φ_n(−i/2)) = φ_m(i/2)·φ_n(i/2) ,
   h(−i/2) = φ_m(−i/2)·φ_n(−i/2) = conj( h(i/2) ) ,
   pole term = 2 Re[ φ_m(i/2) φ_n(i/2) ] .
```

> **Recorded error, fixed before any result (Retraction 21).** The first version of
> `src/suzuki.py` wrote the pole term as `φ_m(i/2)conj(φ_n(i/2)) + φ_m(−i/2)conj(φ_n(−i/2))`.
> That is **also real**, and for a *real* basis it coincides with the correct value, so the
> repository's own real-basis code cannot distinguish the two — only the derivation above
> does. It is wrong for the complex exponential basis. Cross-check: the general formula
> reproduces the repository's `pole_sign = −1` for the odd sector automatically (for odd
> real `f`, `F(i/2)` is purely imaginary, so `F_jF_k = −\|F_j\|\|F_k\|`), which the
> conjugated version does not.

Also by `(1/2π)∫h(t)e^{−itv}dt = C(v)`:

```
   (1/2π)∫h dt = C(0) ,      (1/2π)∫h(t)cos(t log n)dt = [C(log n) + C(−log n)]/2 ,
   prime term = − Σ_n Λ(n) n^{−1/2} · S(log n) ,     S(v) := C(v) + C(−v) .
```

> **Recorded error, fixed before any result (Retraction 22).** The first version used
> `2Λ(n)n^{−1/2}·S(log n)`, a factor **2 too large**. The repository's
> `zeta_window.build_parts_idx` uses `2Λ/√n · Csym(log n)` with
> `Csym = (C_jk + C_kj)/2 = S/2` for a real basis — i.e. `Λ/√n · S`, which is what the
> derivation gives.

## 3. The correlation, in closed form

With `χ_n(x) = e^{iπnx/a}/√(2a)` on `[−a,a]`, for `|v| ≤ 2a` and `d = m−n`:

```
   C_{mn}(v) = ∫χ_m(x) conj(χ_n(x−v)) dx
             = e^{iπ(m+n)v/(2a)} · (2a−|v|)/(2a) · sincpi( d(2a−|v|)/(2a) ) ,
```

(the overlap of `[−a,a]` with `[−a+v, a+v]` has length `2a−|v|` and centre `v/2`),
`C_{mn}(0) = δ_{mn}`, and **every correlation vanishes at `|v| = 2a`** — which is Suzuki's
remark that equality at the prime cutoff contributes zero. Symmetrizing, and using
`sin(πd ∓ κv) = ∓(−1)^d sin(κv)` with `κ = πd/(2a)`:

```
   d ≠ 0:   S(v) = −((−1)^d/(πd))·[ sin(πmv/a) − sin(πnv/a) ]
   d = 0:   S(v) = 2 cos(πmv/a)·(2a−|v|)/(2a) ,      S(0) = 2
```

**`S` is real in both cases**, which (with the real pole term) *proves* `G` is real — the
first clause of P0-T4, derived rather than tested.

## 4. The archimedean term in closed form — why truncation is not an option

```
   arch = −γ_E·C(0) + Σ_{k≥0}[ C(0)/(k+1) − ∫_0^{2a} S(v) e^{−c_k v} dv ] ,   c_k = 1/2+2k
```

The bracket is `O(1/k²)`, so reaching Suzuki's entrywise tail `< 10^{−188}` by truncation
would need `k ~ 10^{188}` terms. **Truncation is impossible; closed form is mandatory.**
Using `Σ_{k≥0} e^{−c_k v} = e^{−v/2}/(1−e^{−2v})`:

```
   d ≠ 0 (so C(0) = 0):
      arch = ((−1)^d/(πd)) ∫_0^{2a} [sin(πmv/a) − sin(πnv/a)] · e^{−v/2}/(1−e^{−2v}) dv

   d = 0:
      arch = ψ(1/4) + e^{−a}·Φ(e^{−4a}, 1, 1/4)
             − ∫_0^{2a} [S(v) − 2] · e^{−v/2}/(1−e^{−2v}) dv
```

using `Σ_{k≥0}[1/(k+1) − 1/(k+1/4)] = ψ(1/4) + γ_E` and
`Σ_{k≥0} e^{−2a c_k}·2/c_k = e^{−a}Φ(e^{−4a},1,1/4)`. Both integrands are `O(v)` against
the `1/(2v)` singularity, hence **analytic on `[0,2a]`**: one non-oscillatory
one-dimensional integral per entry, no series tail, no `10^{−188}` budget to defend.

This is a different route from Suzuki's (he truncates the archimedean series and bounds
the tail) and from the repository's (`arch_closed`, which uses digamma/Lerch/Hurwitz on a
`sin(αu+β)` decomposition). Agreement with the repository's route is P0-T2.
