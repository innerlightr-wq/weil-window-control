# SOURCES_AND_DEPENDENCIES.md

Dated **2026-10-07**. New files live only in this directory. All other tracked and untracked
content of `weil-window-control` is treated as read-only. The `thales-ramanujan-quarter-gamma`
repository was **not accessed** in this pass.

## Repository state at the start of this pass (= at the end)

| | value |
|---|---|
| HEAD | `13ad9dd` "PROJECTION_CONTEXT: add the interpretive principle (magnitude, alignment, response)" |
| branch | `main`, 0 ahead / 0 behind `origin/main` |
| remote | `https://github.com/innerlightr-wq/weil-window-control.git` |
| `git status --short` | `?? extensions/intrinsic-positivity-bridge-2026-10-07/` (pre-existing, untracked) |

Nothing was reset, fetched, pulled, switched, staged, committed or pushed. The local manuscript is
**not** ahead of the remote: `paper/weil_window_control.tex` sha256 `bd71ec14f4044455…`,
`paper/weil_window_control.pdf` sha256 `4ed8d7eef6d3c7d1…` (md5 `a29f9670844f341c7986cfc1d4382bbe`,
byte-identical to the Zenodo deposit `10.5281/zenodo.23212998`).

## Carried forward, with provenance (not re-audited)

| prior result | where | carried as |
|---|---|---|
| `W_∞` is the repo's `Arch − LogPi` block, **not** the whole `Q_ζ`; the repo's form is `M = Pole + Arch − LogPi − Prime` with `Prime` built from von Mangoldt terms `log n < 2L` | `src/zeta_window.py:148,214`; `intrinsic-positivity-bridge/MECHANISM_DICTIONARY.md` §3 | **independently re-checked** against `src/zeta_window.py` in this pass |
| CC2021 archimedean chain: positive Sonin trace, Thm 4.6 identity, Lemma 6.9/6.10, Thm 6.11, `c = 4γ/log2 = 16.98658…` | `intrinsic-positivity-bridge/POSITIVITY_DEFECT_AUDIT.md` §§1–3 | **carried**; only the four dependencies in §2B below were rechecked |
| Theorem 7 of the manuscript is a **relative** estimate; "`λ*→0`" was circular (Retraction 18) | `openai-quasi-rh-audit-2026-10-07/`; memory `weil-window-openai-quasirh-audit` | **circularity safeguard carried**: no step below normalises by an unproved gap of the Weil form or by a nonexistent bounded `L²` operator norm of it |
| `P_old = J_e S_e/m_e` commutes with the **full** operator and needs BOTH intertwinings; `P_E` onto `ker A_0` fixes only a leading coefficient; `η = 0` iff `ker(G)` holds a **nonzero-mean** vector | `second-moment-generalization-2026-10-07/PROJECTION_CONTEXT.md` §5 (as corrected) | **diagnostic motivation only.** No step of the derivation below uses a finite-matrix projected-moment or quartic result. Kept strictly separate, per §6 of the brief. |
| `pdftotext` drops `\left[`/`\right]`; `str(Decimal)` slicing drops exponents | `second-moment-generalization-2026-10-07/` | every displayed formula below was read in layout output **and** cross-checked against the surrounding derivation; all printed numbers use `'%.*e' % float(x)` |

Not re-run, per the brief: quasi-RH transfer, prime-comb, scalar coupling, nested-window
inheritance, Ramanujan analogy.

## A. Connes, *Trace formula in noncommutative geometry and the zeros of the Riemann zeta
function*, Selecta Math. (N.S.) **5** (1999), 29–106 — **the previously unchecked dependency**

**READ** from `arXiv:math/9811068v1`, 10 Nov 1998, 88 pp (519,632 B; `pdftotext -layout`).
Citations give **preprint pagination**. Sections read in full: **VII** (pp. 28–36, the `S`-local
trace formula) and **VIII** (pp. 37–46, the global case). Sections I–VI and Appendices I–III
**skimmed only — NOT CHECKED in detail**; nothing below depends on them.

Label convention: **C99-Thm VII.4**, **C99-Lem VIII.1**, **C99-Cor VIII.2**, **C99-Thm VIII.5**.
These are distinct from CC2021-Thm 4.6 and from CCM2024-Thm 4.6.

### Setup (C99 §VII, pp. 28–30)

`k` a global field, `S` a finite set of places containing **all infinite places**;
`O_S^* = {q ∈ k^* : |q_v| = 1, v ∉ S}` the `S`-units; `J_S = ∏_{v∈S} k_v^*`, `C_S = J_S/O_S^*`;
`A_S = ∏_{v∈S} k_v` and `X_S = A_S/O_S^*`. `L²(X_S)` is the completion of `S(A_S)` for
`‖f‖² = ∫ |Σ_{q∈O_S^*} f(qx)|² |x| d^*x` (eq. 5, p. 29). Scaling `(U(λ)ξ)(x) = ξ(λ^{-1}x)` (eq. 7);
`U(h) = ∫ h(g)U(g)dg` for `h ∈ S(C_S)` compactly supported.

**The two cutoffs (eqs. 12–13, p. 31), verbatim:**

> `P_Λ = {ξ ∈ L²(X_S) ; ξ(x) = 0  ∀x , |x| > Λ}` … "This gives an infrared cutoff and to get an
> ultraviolet cutoff we use `P̂_Λ = F P_Λ F^{-1}`" … `R_Λ = P̂_Λ P_Λ`.

`|x|` is the **total module** on `X_S`, not a per-place condition.

### C99-Thm VII.4 (printed p. 31) — the proved semi-local trace formula

> "Let `A_S` be as above, with basic character `α = ∏ α_v`. Let `h ∈ S(C_S)` have compact support.
> Then when `Λ → ∞`, one has
> `Trace(R_Λ U(h)) = 2h(1) log′Λ + Σ_{v∈S} ∫′_{k_v^*} h(u^{-1})/|1−u| d^*u + o(1)`"

with `2log′Λ = ∫_{λ∈C_S, |λ|∈[Λ^{-1},Λ]} d^*λ`, each `k_v^*` embedded in `C_S` by
`u ↦ (1,…,u,…,1)`, and the principal value `∫′` fixed by the unique distribution agreeing with
`du/|1−u|` off `u=1` whose `α_v`-Fourier transform vanishes at 1.

**Hypotheses and scope, stated exactly.** Valid for **any** global field `k` and any such `S`,
hence for `k = Q`, `S = {∞,2}`. It is an **asymptotic** statement in the cutoff `Λ`, with an
explicitly divergent leading term `2h(1)log′Λ`. The arithmetic side is exactly the sum of the
**local Weil distributions over the places of `S`**. **No convolution-square structure and no
positivity are asserted or available**: the operator is `R_Λ = P̂_Λ P_Λ`, a product of two
orthogonal projections, which is **neither self-adjoint nor positive** unless they commute.

### C99-Lem VIII.1 (printed pp. 38–39) and C99-Cor VIII.2 (printed p. 41)

**C99-Lem VIII.1:** "Let `χ_0` be a character of `C_{S,1}`, then **for `Λ` large enough `P̂_Λ` and
`P_Λ` commute** on the Hilbert space `L²_{χ_0}`."

**C99-Cor VIII.2:** "We can thus rewrite Theorem 4 **in the case of positive characteristic** as,
*Corollary 2.* Let `Q_Λ` be the **orthogonal projection** on the subspace of `L²(X_S)` spanned by
the `f ∈ S(A_S)` which vanish **as well as their Fourier transform for `|x| > Λ`**. … Then when
`Λ → ∞`, `Trace(Q_Λ U(h)) = 2h(1)log′Λ + Σ_{v∈S} ∫′_{k_v^*} h(u^{-1})/|1−u| d^*u + o(1)`."

**The restriction is essential and is not cosmetic.** The proof of C99-Lem VIII.1 is a function-field
counting argument throughout: it sets `Λ = q^N`, uses `2log′Λ = (2N+1)log q`,
`λ = l/log(q)`, the local principal value `∫′_{R_v^*} χ_v(u)/|1−u| d^*u = −f_v log(q_v)` with `f_v`
the ramification order, `|d| = ∏|d_v| = q^{2−2g}` with **`g` the genus of the curve**, and the
exact count `(2N+1)l − fl + (2−2g)l` of elements `g ∈ D_S` with `|g| ≤ Λ`, `|g^{-1}| ≤ q^{2−2g−f}Λ`.
It constructs explicit vectors `η_χ = ∏_S φ_v` ([W1] Prop. VII.13) and notes "Since `ξ` is locally
constant, its Fourier transform has compact support" — a non-archimedean phenomenon. The §VIII
hypothesis is also explicit: `S ⊃ S_0` finite "large enough so that `mod(C_S) = mod(C_k) = q^Z`."

### C99-Thm VIII.5 (printed pp. 41–42) and the global formula (16)

The global analogue (16) of the trace formula is **not proved**: "We can prove directly that (16)
holds when `h` is supported by `C_{k,1}` but **are not able to prove (16) directly for arbitrary
`h`**." C99-Thm VIII.5 then shows, **for `k` of positive characteristic**, that (16) `⟺` RH for all
`L`-functions with Grössencharakter on `k`.

### SOURCE GATE — what remains missing after reading A

1. **In characteristic zero the arithmetic identity is attached to a non-self-adjoint operator.**
   For `k = Q`, `S = {∞,2}`, the proved statement is C99-Thm VII.4 about `R_Λ = P̂_Λ P_Λ`. The
   upgrade to the **orthogonal** projection (C99-Cor VIII.2) is stated and proved **only in positive
   characteristic**, through a genus/degree count that has no characteristic-zero counterpart. **No
   checked source proves `[P̂_Λ, P_Λ] = 0` for `k = Q`.**
2. **The cutoff convention is complementary to the Sonin convention.** C99's `B_Λ` is
   `{f : f(x)=0 & f̂(x)=0 for |x| > Λ}` — doubly **low**-pass. The Sonin space of CC2021 and
   CCM2024 is `{f : f(x)=0 & F f(x)=0 for |x| < λ}` — doubly **high**-pass. In the two-projection
   decomposition of the same pair these are the opposite corners (`ran ∩ ran` versus `ker ∩ ker`).
   **C99-Thm VII.4 is therefore not an arithmetic identity for the Sonin projection**, and no
   checked source supplies the translation.
3. **No ordinary positive trace.** C99's trace is regularised by the `Λ`-cutoff and carries the
   divergent `2h(1)log′Λ`. Per the gate, no sign is assigned to it by analogy.
4. **Still NOT CHECKED:** Weil [W1] *Basic Number Theory* Prop. VII.13 and Appendix IV, used by
   C99-Lem VIII.1; C99 §§I–VI and Appendices; Yoshida. Used only as cited.

## B. Connes–Consani, *Weil positivity and Trace formula, the archimedean place*

**READ** (cached): author PDF dated **July 4, 2021**, `alainconnes.org/wp-content/uploads/Selecta.pdf`,
55 pp, 935,906 B. Related record arXiv:2006.13771 **downloaded, NOT used for numbering**. Labels
below are **CC2021-**.

Only the four dependencies the brief names were rechecked:

- **CC2021-Thm 4.6 (p. 24, eqs. 84–85).** `Tr(ϑ(f)S) = W_∞(f) + ∫ f(ρ^{-1})ε(ρ) d^*ρ` for **all**
  `f ∈ C_c^∞(R_+^*)` — **no support restriction**; `ε` given independently by (85), `ε(ρ^{-1})=ε(ρ)`.
  Write `E(f) := ∫ f(ρ^{-1})ε(ρ)d^*ρ`. **Independently defined, not by subtraction.**
- **Positive Sonin trace.** `S = Π_{S(1,1)}` is the orthogonal projection onto the Sonin space of
  even functions vanishing on `[−1,1]` together with their Fourier transform;
  `ϑ(λ)ξ(v) = λ^{-1/2}ξ(λ^{-1}v)` on `L²(R)_ev`. Positivity is the Gram identity
  `Tr(ϑ(g)Sϑ(g)^*) = ‖ϑ(g)S‖²_{HS} ≥ 0`, which needs **no** `ϑ`-invariance of `ran S`.
- **Support/vanishing hypotheses.** CC2021-Thm 6.11 (p. 45, eq. 142): `g` supported in
  `[2^{-1/2},2^{1/2}]`, `ĝ(−i/2) = 0` — **one** condition — gives
  `W_∞(g*g^*) ≥ Tr(ϑ(g)Sϑ(g)^*) − c|ĝ(0)|²`, `c = 4γ/log 2`. CC2021-Thm 1 (p. 2) adds `ĝ(0)=0`,
  i.e. **two** conditions.
- **Rank-one majorant / alignment.** CC2021-Lem 6.10 (p. 44, eq. 141): with
  `N_I = −2ℓ_1(1+ε_1)(Id − K_I)` on `H = L²(I)`, `I = [−½log2, ½log2]`, `γ ≈ 2.94355`,
  `⟨ξ|N_I ξ⟩ ≤ γ|⟨ξ_0|ξ⟩|²`. CC2021-Lem 6.9 (p. 44, eqs. 139–140) is the alignment criterion
  (`a + c ≥ b` **and** `b(a+c) ≤ a(b+c)|⟨φ|ψ⟩|²`). CC2021-Rem 6.12 quotes `λ_max = 1.05158`,
  `ε_1 ≈ 0.00122`, **not** independently derived here.

**Both of CC2021-Lem 6.10 and CC2021-Thm 6.11 are tied to `I = [−½log2, ½log2]`** and are not
extended to a larger interval anywhere in the checked sources. No step below extends them silently.

## C. Connes–Consani–Moscovici, *Zeta zeros and prolate wave operators*, arXiv:2310.18423v2

**READ** (cached PDF, 4 May 2024, 30 pp). Labels **CCM2024-**. §4 read in full; §§2–3, 5–6
**NOT CHECKED in detail**. The five load-bearing items, with printed pages:

- **CCM2024-Def 4.5 (p. 21).** `S_λ(X_S,α) := {f ∈ L²(X_S)^{K_S} | f(x) = 0 & F_S f(x) = 0  ∀x, |x| < λ}`.
- **CCM2024-Prop 4.6 (pp. 21–22).** (i) `θ_S(f) :=` class of `σ_S ⊗ f`; `θ_S(S_λ(R,e_∞)) ⊂ S_λ(X_S,α)`.
  (ii) eq. (57): `F_μ(w_S(θ_S f))(s) = F_μ(w_∞ f)(s) × ∏_{S∖{∞}} (1 − p^{−½−is})`.
  **From the proof, with `σ_p = ε_0 − (1/p)ε_1` and `g = w_∞(f)`:**
  `w_S(σ_S ⊗ f)(1×λ) = g(λ) − p^{−1/2} g(λ/p)` and `F_μ(g_1)(s) = p^{-is}F_μ(g)(s)` for
  `g_1(λ) = g(λ/p)`. So **multiplication by `p^{-is}` is scaling `λ ↦ λ/p`** — the fact that makes
  the multiplier an operator identity rather than a symbol.
- **CCM2024-Prop 4.7 (p. 22).** (i) `F_S ∘ θ_S = θ_S ∘ F̃_R`. (ii) the commutative diagram (59)
  with the **two different measures** `ds/|E_∞(s)|²` and `ds/|E_S(s)|²`.
  (iii) `⟨θ_S(f)|η_S(g)⟩ = ⟨f|g⟩`. **Its proof states `U_S = F_μ w_S : L²(X_S)^{K_S} → L²(R)` is
  unitary.** (iii) gives `θ_S^* η_S = Id`, i.e. `η_S = (θ_S^*)^{-1}` — **not** `θ_S^{-1}`.
- **CCM2024-Thm 4.6 (p. 23) = intro Theorem 2 (p. 5).** "the map `θ_S` is a **hilbertian
  isomorphism** of the Sonin spaces `θ_S : S_λ(R,e_∞) → S_λ(X_S,α)`." The proof establishes
  **surjectivity** ("To show that one has equality, let `h ∈ S_λ(X_S,α)` …"), so
  `S_λ(X_S,α) = θ_S(S_λ(R,e_∞))` exactly.
- **Footnote 2 (p. 5), verbatim.** "We use the term "hilbertian" to denote the underlying
  topological vector space structure of a Hilbert space." So bounded with bounded inverse
  (CCM2024-Prop 4.3(i)), **not isometric, not unitary**; and §4.8: "the choice of the finite set `S`
  plays a key role in **fixing the inner product**", `B_λ` "inherits **different inner products**".

**Not established in CCM2024:** the semi-local positive-trace comparison. §1 says only "We
**expect** that the use of such operator-theoretic tools in the semilocal case **opens a way** to
handle Weil's positivity as in [10]" and offers "a more precise **strategy**"; the archimedean
identification holds "up to a finite dimensional possible discrepancy"; the abstract defers a second
candidate prolate operator to "a forthcoming paper". **PROPOSED, not proved.**

## Later literature

Searched only for a directly relevant continuation. **Not found in the sources checked:** any
semi-local analogue of CC2021-Thm 4.6 or CC2021-Thm 6.11; any evaluation or bound of the principal
angles between `P_Λ` and `P̂_Λ` on `L²(X_S)` for `k = Q`; any characteristic-zero proof of
`[P̂_Λ,P_Λ]=0`. **NOT CHECKED / NOT LOCATED:** "New eigenfunctions for the negative part of the
Connes–Moscovici prolate spectrum" (2025 item, not retrieved); the forthcoming CCM paper. **This is
not an exhaustive literature review**; absence here means "not found in the sources checked".

## Added in the post-review amendment pass (2026-10-07)

| item | status |
|---|---|
| Poisson-kernel / Abel-summation identity `(1 − 2r\cosθ + r²)^{-1} = (1−r²)^{-1}(1 + 2Σ_{k≥1}r^k\cos kθ)`, `0<r<1` | **standard**; verified here at `r = 2^{-1/2}`, `θ = s\log2` to `9.0e-14` (`amend_checks.py` §A). It is what retracts the reach inference: `1/M` has **infinite** Laurent support where `M` has three terms |
| `\log M_σ(s) = −2Σ_{k≥1}k^{-1}2^{-kσ}\cos(ks\log2)` and `∂_σ\log M_σ = 2\log2 Σ_{k≥1}2^{-kσ}\cos(ks\log2)` | **standard** expansion of the local Euler factor (cf. DLMF §25.2 for the Euler product); verified to `3.3e-10` (`amend_checks.py` §B). It is what withdraws the parity reading. **No geometric operator family off `σ = ½` is constructed or claimed** |

Neither item is a new arithmetic result, and neither supplies the missing trace identity.

## Tooling

Python 3 standard library only (`fractions.Fraction`, `decimal.Decimal`, `math`). `mpmath` is
**not installable on this box** (no `pip`/`uv`/`pacman`), so the repository's `src/` regression
scripts were **not run**; nothing below depends on them. No background jobs, no installs, no scans.
