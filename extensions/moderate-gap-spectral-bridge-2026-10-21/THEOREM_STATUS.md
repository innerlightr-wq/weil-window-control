# THEOREM_STATUS — every statement in this package, with class and tier

**Provenance classes** (brief §13): **I** classical Fourier/matrix identity ·
**II** proved by ARTZ or an antecedent (Lamzouri, Alpöge–Furman, Goldston–Sahay,
BGSTB, Montgomery) · **III** previously established in my RH program ·
**IV** exact new consequence of combining the frameworks, proved here ·
**V** numerical observation made here · **VI** conditional or open.

**Tiers:** `[certified]` exact rational / rigorous interval arithmetic ·
`[proved]` pen-and-paper proof given here · `[diagnostic]` double precision ·
`[cited]` taken from the literature, not re-derived.

---

## A. The new paper, reproduced

| # | statement | class | tier | where |
|---|---|---|---|---|
| A1 | `C_MT = 1/2 + (1/√2)cot(1/√2) = 1.3274992963205883…` (width `3.3e−40`) | II | certified | 02 §5, verify A1 |
| A2 | `C_MT < 4/3`, hence the four-gap bound is positive | II | certified | verify A2 |
| A3 | the three proportions `0.3450014…, 0.1812518…, 0.0175021…` | II | certified | verify A3 |
| A4 | `ℓ = 5` gives `−0.146248`; the ceiling is `ℓ < (4−2C_MT)/(C_MT−1) = 4.1069` | II | certified | verify A4, A5 |
| A5 | the ceiling is set by `C_MT`, **not** by the kernel bound (`ℓ/(ℓ+2) < K²` holds to `ℓ=6`) | IV | certified | verify A6 |
| A6 | `A(x) = cos(√2 x)/(√2 sin(1/√2))` is a probability density on `(−1/2,1/2)` | II | certified | verify B1 |
| A7 | corrected `Â(c)` formula (pdftotext loses `1/` and the denominators of (70)); `Â(0)=1`, `Â(1/4)=0.907203020832` | IV (correction) | certified | 02 §2, verify B2, B3 |
| A8 | `A_τ ≥ 1−2τ ≥ 99/100` at `τ=1/200`; `K_τ(1/4) > √3/2`; `K_τ(c)² > 3/4` | II | certified | verify B4–B7 |
| A9 | `C_η(A) = (A∗A)(0)+2∫₀¹α(A∗A)(α)dα = C_MT` | II | diagnostic (1.7e−7) | verify B8 |
| A10 | Lemma 4: `(p−q)² ≥ Ψ(p)−4q`, with both case identities | II | certified | verify C1–C3 |
| A11 | Lemma 10: `trΨ = 2+2K²` exactly | II | certified | verify C4 |
| A12 | Lemma 11's step `4K−1 ≥ 3K²` on `[1/3,1]`, i.e. `−(3K−1)(K−1) ≥ 0` | II | certified | verify C5 |
| A13 | `tr G = N(T)` and `‖G‖²_F = Σ_{z,s}K_τ(z−s)²`, on 7 synthetic multisets | II | diagnostic | verify F1–F3 |
| A14 | Lemma 3.1 = `∫_{-1}^{1}f(u)F(T^u,T)du`; the `f(0)` term comes from `x` near 1, the `2∫₀¹αf` term from `log x` | IV | proved | 02 §5 |

## B. The Hilbert-space structure and the depth scalar

| # | statement | class | tier | where |
|---|---|---|---|---|
| B1 | `Im z_ρ = δ_ρ·log T/(2π)`: in ARTZ's coordinate the displacement is the imaginary part | IV | proved | 02 §1 |
| B2 | closed forms for `Re K`, `Im K` at complex argument (`η²≡1`), checked against Simpson to `5.3e−15` | I | diagnostic | verify D1 |
| B3 | `κ_ρ := K_τ(2i Im z_ρ) = ∫η_τ²(u)cosh(2δ_ρ u log T)du ≥ 1` | IV | proved | 02 §3 |
| B4 | `‖g_z‖² = (1+κ)/2`, `‖h_z‖² = (κ−1)/2`, `⟨g_z,h_z⟩ = 0`, `‖g‖²−‖h‖² = 1` | IV (last part II) | certified | verify D4–D7 |
| B5 | the full real Gram table for `{f_x, g_z, h_z}` at arbitrary depths and separations | IV | certified | 02 §3 |
| B6 | the bilinear pairing gives `K(z−s)`, the Hermitian one `K(z−s̄)`; they differ off the line and the totals agree only because `Z_T = Z̄_T` | IV | proved | 03 §2 Q8 |
| B7 | off-line block spectrum `{m(1+κ), −m(κ−1)}`; `tr = 2m` for every `κ`; `‖·‖²_F = 2m²(1+κ²)` | IV | certified | verify E1, E2 |
| B8 | **the trace and the inertia are exactly the `δ`-blind invariants** | IV | certified | 04 §2 |

## C. The new inequality and its scope

| # | statement | class | tier | where |
|---|---|---|---|---|
| C1 | `tr Q + K_τ(c)²\|R̃_1\| + Σ_{q_j>0}(q_j−2)² ≤ ‖G‖²_F − N(T)` | **IV** | proved; checked on 7 multisets | 04 §1, verify E4, F4 |
| C2 | the discarded slack for a simple off-line pair is exactly `(κ−1)²` | IV | certified | verify E3 |
| C3 | `κ−1 = 2δ²V_μ + O(δ⁴)`, `V_μ` = centred second moment of the window measure | IV | certified | verify E5, I1–I3 |
| C4 | `(κ−1)² = 4V_μ²δ⁴ + O(δ⁶)` — exponent **4** | IV | certified | verify E6 |
| C5 | `‖G‖²_F ≥ ‖Q_+‖²_F` (elementary: `‖P−Q_-‖²_F ≥ 0`, `tr(PQ_+) ≥ 0`) | I | diagnostic check | verify F5 |
| C6 | single-pair bound `2m + (m(1+κ)−2)² ≤ (C_{η_τ}−1+o(1))N(T)` | IV | proved | 04 §3 B1 |
| C7 | `κ ≤ cosh(δ log T)` for **every** admissible window; `δ ≤ 1/2` from `β ≤ 1` | IV | certified | verify G1, G2 |
| C8 | pointwise deficit `= (2/π)(C_MT−1)log(T/2π) = 0.2084930…·log(T/2π) > 1` for `T > 761` | IV | certified | verify G3–G5 |
| C9 | exponent identity `1−2δ = 2(1−σ)`: the summed form **is** the Density Hypothesis | IV | certified | verify G6 |
| C10 | that exponent beats Selberg 1946 and Ingham for every `σ ∈ (1/2,3/4]` | IV | certified | verify G7 |

## D. Conditioning and localisation

| # | statement | class | tier | where |
|---|---|---|---|---|
| D1 | **anti-localisation**: every interior ordinate of a run of `ℓ` positive gaps is a simple on-line zero; endpoints need not be | **IV** | proved; 8704 runs, 0 violations | 05 §1, verify J1–J5 |
| D2 | a 4-gap run certifies **three**, not five, simple critical-line zeros | IV | proved | 05 §1 |
| D3 | off-line endpoints occur (2685 of 8704 runs in the random test) | V | diagnostic | verify J2 |
| D4 | `λ_min(G_d) ≤ (2sin(πc/2))^{2(d−1)}/\binom{2d−2}{d−1}` | **IV** | certified | 05 §2, verify H1–H4 |
| D5 | decay rate `2sin²(πc/2) = 0.2928932…` at `c = 1/4`; decay for every `c < 1/2` | IV | certified | verify H1, H3 |
| D6 | at `c = 1` the nodes are orthonormal (Nyquist spacing of the band-limited kernel) | I | diagnostic | verify H5 |
| D7 | `λ_min` for `d ≤ 2` certified: `≥ 1 − 0.9274 = 0.0726` at `c = 1/4`, `τ = 1/200` | IV | certified | 05 §3 Q2 |
| D8 | adversarial **on-line** placements do not help the adversary (they raise `trΨ(P)` from 24 to 32.05) | V | diagnostic | 05 §4, verify F case 4 |

## E. The non-additivity law

| # | statement | class | tier | where |
|---|---|---|---|---|
| E1 | the slack is **not** additive over off-line pairs | **V** | diagnostic | 04 §3.1 |
| E2 | at depth `a = 2π`, separation `2`: `total/(κ−1)² = 0.9799 + 1.3315 ln g`, max residual `0.0880` for `g = 2..32` | V | diagnostic | results/additivity.txt |
| E3 | best linear fit at the same data has max residual `1.0402` — 12× worse | V | diagnostic | results/additivity.txt |
| E4 | at separations `0.5–0.7`, 16 pairs cost `0.589–0.657` of **one** isolated pair | V | diagnostic | results/cluster_smooth.txt |
| E5 | mechanism: `g`-direction correlation `−0.751` at separation 1, depth `a=2π`; `−0.884` at `a=3π` | IV (mechanism) + V (numbers) | proved + diagnostic | 04 §3.2 |
| E6 | the fraction of the slack carried by the top eigenvalue falls `1.000 → 0.312` as `g: 1 → 32` | V | diagnostic | results/additivity.txt |
| E7 | on-line neighbours at integer spacing contribute only `4(1−μ²)`, via `Σ_{d∈Z}(y²−d²)/(y²+d²)² = π²/sinh²(πy)` | I | proved, diag-corroborated | verify F8, F9 |
| E8 | **had the slack been additive, C1 would prove the Density Hypothesis** — so non-additivity is forced | IV | certified (from C9, C10) | 04 §3, 09 §5 |

## F. Cross-framework statements about my own program

| # | statement | class | tier | where |
|---|---|---|---|---|
| F1 | the `δ²` coefficient is "twice the centred second moment of the window" in both frameworks | IV | certified | 07 §4, verify I2, I3 |
| F2 | `2Σ_{i=0}^{R}(i−R/2)² = \binom{R+2}{3}` exactly for `R = 1..200` | III | certified | verify I4 |
| F3 | the coefficient scales as `L^{s+2}` for window mass `~L^s`: ARTZ is `s=0` (`L²`), my note is `s=1` (`L³`) | IV | certified | verify I5, I6 |
| F4 | exponent **4** agrees with the independently certified `M(δ) = δ⁴D₂`, `D₂ = 2911.050006543213352847841125028374273138 ± 1e−37` | IV | certified (here) + certified (prior) | 03 §3 |
| F5 | **correction**: `NG-2'` (margin exactly 0 at multiplicity `ν ≥ 2`) does **not** apply verbatim to an off-line pair; detectability is `δ⁴D₂ > 0`, not 0 | III, corrected | proved | 06 §2 |
| F6 | ARTZ's multiplicity weights are **scalar weights on independent vectors, not a jet** | III (preserved) | proved | 06 §7 |
| F7 | the method's extremal configuration (off-line pairs of depth `→ 0`) **is** the detector's blind spot | IV | proved from II + III | 06 §5 |

## G. Conditional and open

| # | statement | class | where |
|---|---|---|---|
| G1 | moderate gaps give `κ(Γ) ≳ c^{2n}D_n/(n!)²` at a positive proportion of indices — **conditional on the unfinished Hermite–Genocchi lower bound** | VI | 06 §1 |
| G2 | `(H-add)` ⟹ `N(σ,T) ≪_ε T^{2(1−σ)+ε}` — hypothesis measured false and DH-equivalent | VI | 10 T6 |
| G3 | the exact growth order in `g` of `sup Σ_j(q_j−2)²` under separation `≥ s` — **open** | VI | 11, next step 2 |
| G4 | whether BGSTB Lemma 5's error term is uniform over configurations with large symmetric-diagonal terms — **not checked** | VI | ASSUMPTIONS D1 |

## Refuted here

| statement | why |
|---|---|
| "moderate gaps give a *uniform in `d`* Gram lower bound under separation `c`" | D4/D5: exponential decay |
| "the slack `Σ_j(q_j−2)²` is additive over off-line pairs" | E1–E4 |
| "an off-line pair is a near-collision, so `NG-2'` gives margin 0" | F5 |
| "the gap theorem localises a run near a prescribed ordinate" | D1 |
| "the `δ²–L³` law is universal; the paper exhibits `L³`" | F3: the paper exhibits `L²` |
| "a run of four moderate gaps gives five simple critical-line zeros" | D2 |
