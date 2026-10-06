# Prior art — verified against full texts (Round 3, Task 1)

**Method.** Unlike the Round-2 novelty check (abstracts only, and flagged as such), every
item below was checked against the **full text**, downloaded and read locally. Quotes are
verbatim and ≤25 words. Where a claim could not be confirmed, that is stated.

**Headline: the prior art is stronger than Round 2 found, and it changes the claim list.**
Three of the results I had treated as contributions are published theorems.

---

## 1. Hallouin & Perret (2019) — the Toeplitz Gram matrix *is* my Stage 1 object

**Full citation (verified from the published PDF's own header):**

> Emmanuel Hallouin and Marc Perret, *A unified viewpoint for upper bounds for the number
> of points of curves over finite fields via Euclidean geometry and semi-definite symmetric
> Toeplitz matrices*, Trans. Amer. Math. Soc. **372** (2019), no. 8, 5409–5451.
> DOI: 10.1090/tran/7813. Published electronically June 17, 2019.

*Obtained from* `https://www.math.univ-toulouse.fr/~hallouin/Documents/eh-Weil.pdf`
(author's copy of the AMS offprint). Not on arXiv under this title; the 2014 preprint
arXiv:1409.2357, *From Hodge Index Theorem to the number of points of curves over finite
fields*, is an earlier, differently-titled and differently-numbered version.

> **Citation caveat, verified.** The authors' own later paper (arXiv:2005.12190,
> bibliography entry `[HP19]`) prints the volume as "312". The published header, and
> CvS's bibliography entry [12], both give **372**. 372 is correct; 312 is a typo. Use 372.

### (a) Proposition 5 — CONFIRMED, and stronger than stated in the task

§1.2, pp. 5415–5416. With `γ^k = p(Γ_k)/√(q^k)` (their Definition 4):

> "⟨γ^i, γ^j⟩ = 2g if i = j, x_{|i−j|} if i ≠ j"

and

> "x_i = ⟨γ^0, γ^i⟩ = ((q^i + 1) − |X(F_{q^i})|)/√(q^i)"

their equation (5) then displays `Gram(γ^0,…,γ^n)` as a Toeplitz matrix with `2g` on the
diagonal. **This is my `T_R` entry for entry — not merely up to scale.** My `t(0) = 2g` and
`t(n) = (q^n + 1 − N_n)/q^{n/2}` are their `2g` and `x_n`. The task statement omitted the
`2g` diagonal; with it, the match is exact.

*(The 2014 preprint used an extra `1/√(2g)`, normalising the diagonal to 1; the published
version does not.)*

### (b) Weil bound from the 2×2 corner — CONFIRMED

§2.5.1, p. 5423, titled "Weil bound as a bound of order 1":

> "the usual Weil bound comes from Cauchy-Schwartz inequality applied to γ^0 and γ^1"

from `Gram(γ^0,γ^1) ≥ 0`, i.e. `|x_1| ≤ 2g`, giving `||X(F_q)| − (q+1)| ≤ 2g√q`.

### (c) Theorem 6 / §1.3 — CONFIRMED verbatim, including the kernel generator

§1.3, "Toeplitz interpretation of Riemann hypothesis for curves over finite fields",
pp. 5417–5418. **Theorem 6 (Toeplitz version of Riemann hypothesis for curves).** Its proof
states:

> "the dimension d of F corresponds to the minimal integer such that Gram(γ^0,…,γ^d) is
> singular"

and

> "Let (a_0,…,a_d) ∈ Q^{d+1} be a generator of the kernel of this Gram matrix"

**This is my H2, restatement included.** Their `d` = rank of the Frobenius space
= degree of the minimal polynomial of Frobenius = my `d` = number of distinct `θ_j`.

### (d) Appendix A.2 — this is my H1, and my Round-2 rank method

- **Theorem 36(ii)**, p. 5448, for a real symmetric PSD Toeplitz matrix of size `n+1` and
  rank `n` with kernel vector `(a_0,…,a_n)` and `P = Σ a_j X^j`:
  > "The polynomial P(X) has n distinct complex roots ε_1,…,ε_n ∈ C all of norm 1."

  **That is H1 exactly** (Carathéodory–Fejér / Pisarenko).
- **Theorem 36(i)**: the shift `γ_i → γ_{i+1}` is an isometry with `P` as **minimal and
  characteristic** polynomial — from which my "kernel is the ideal `(P*)`" is immediate.
- **Theorem 36(iv)**: the factorisation `x_k = Σ_i λ_i ε_i^k` with `λ_i > 0` — the atomic
  measure representation, i.e. my `T_R = Σ_m w_m w_m^H` positivity.
- **Lemma 33**, p. 5445: rank of a PSD Toeplitz matrix = size of the largest non-zero
  **leading** minor. This is exactly the leading-minor/LDL rule I used in Round 2.
- **Lemma 34**, p. 5446: the `(γ_i ± γ_{n−i})/√2` decomposition — the parity split I used.

### (e) What is NOT in Hallouin–Perret

Searched the full text: no occurrence of perturbation, planted/fictitious spectra,
off-circle eigenvalues, or any `δ`-expansion. Their object is the **Weil domains** `W_n`
(the convex set where the Toeplitz matrix is positive definite) and upper bounds for point
counts. **They never run the form on data violating `|α| = √q`.** That is the gap this
project occupies.

---

## 2. Connes–van Suijlekom, arXiv:2511.23257 — both task claims CONFIRMED

Introduction, immediately after Corollary 1.1 (the Carathéodory–Fejér corollary):

> "resonating with the analogue of the Riemann Hypothesis for function fields; see [12]
> for a further discussion of this connection"

and, in the next sentence:

> "The key difficulty in this context, then, becomes the verification that zero is indeed
> the (simple) minimal eigenvalue of T."

Their **[12] is Hallouin–Perret**, cited with volume **372**, no. 8, 5409–5451.

**Two further passages that bear directly on my claims:**

- **Remark 2.3**: when the kernel is more than one-dimensional, kernel polynomials vanish
  on the Carathéodory–Fejér atoms "but which otherwise can have arbitrary other zeros",
  and "The correct formulation of the theorem is that if you take the intersection of the
  zeros of the various eigenfunctions, then they are all on the unit circle."
  **This is my Round-2 §3 result, stated in words.** My contribution reduces to making it
  quantitative (the dimension formula `R+1−d` and the parity split).
- **Appendix A**: "an example where the truncation matrices admit simple maximal and
  minimal eigenvalues but this property fails in the limit". CvS already exhibit simplicity
  failing.

---

## 3. Bombieri (2000) — the planted off-line zero experiment is ALREADY THERE

**Full citation (verified from the digitised article):**

> Enrico Bombieri, *Remarks on Weil's quadratic functional in the theory of prime numbers,
> I*, Atti Accad. Naz. Lincei Cl. Sci. Fis. Mat. Natur. Rend. Lincei (9) Mat. Appl.
> **11** (2000), no. 3, 183–233.
> Digitised at `http://www.bdim.eu/item?id=RLIN_2000_9_11_3_183_0` (full text obtained).

**Accessible and checked — not "not checked".**

### (i) Negative-eigenvalue count — FOUND

Abstract:
> "we show that the number of negative eigenvalues is precisely one-half of the number of
> zeros failing to [lie on the critical line]"

**Theorem 8** (§8, "The eigenvalues of finite approximations", p. ~213):
> "The number of negative eigenvalues of the matrix H(Γ;t) equals the number of distinct
> complex conjugate pairs (γ, γ̄) in Γ."

restated for `K_E(Γ)` and separately for the even and odd parts in §9 (pp. ~217–218).
**This is the Stage 2 signature count.** Round 2 guessed it was known; it is Bombieri's
Theorem 8.

### (ii) Perturbative / δ-expansion — NOT FOUND

Searched the full text for any expansion in the distance of a zero from the critical line.
None. The depth is **fixed**, never varied.

### (iii) But Bombieri DOES run the planted experiment — this was missed in Round 2

**§13, "Some numerical experiments", pp. ~227–228:**

> "the first N zeros of ζ(s) with positive imaginary parts and a fictitious zero ρ_0 off
> the critical line, together with their images by complex conjugation and reflection"

— i.e. **exactly my planted quartet `{ρ, 1−ρ, ρ̄, 1−ρ̄}`** — with

> "The fictitious zero off the line has been arbitrarily set at ρ_0 = 0.52 + i3.14"

(so `δ = 0.02`), `N` up to 160, and `E = [−t,t]` — **my window**. He reports

> "the existence of a critical value t_c± > 0 such that the unique negative eigenvalue
> λ_N±(t) tends to 0 if t < t_c±, as N → ∞"

and a natural even/odd division of eigenfunctions.

**Consequence for my claims.** My Stage 3e "detection onset" is a rediscovery of Bombieri's
critical value `t_c`, and my planted-quartet construction is his. What remains is
*quantitative*: the `δ²` scaling (he fixes one `δ`), the `L³` coefficient and its
Cauchy–Schwarz saturation, and the identification of the onset with Zhu's bandwidth
`T*(L) = 2πe^{2L}`.

---

## Revised verdict on novelty

| Round-2 claim | Round-3 verdict after reading full texts |
|---|---|
| Function-field control with planted spectra | **Partly new as an experiment.** The framework is Hallouin–Perret Prop. 5 + Thm 6 + Thm 36; what is not in the literature is running it on *planted* spectra end-to-end against the CvS pipeline. |
| `ker T_R = {P* q}`, simplicity failing past `2g` | **Not new.** HP Thm 36(i) + Thm 6; CvS Remark 2.3 states it in words. Keep only as a quantitative corollary. |
| Toeplitz–Rosati bridge | **Not new.** HP Proposition 5 *is* a Gram matrix in the positive form from the curve's geometry. My `Tr(φψ†)` version is an equivalent realisation, not a new identification. |
| `δ²`–`L³` sensitivity law | **Still appears new.** No `δ`-expansion in HP, CvS or Bombieri; Bombieri fixes `δ = 0.02`. This is now the main quantitative claim. |
| Block F inertia on ζ windows | **Appears new, weakly.** Not found; but it is a computational observation with a stated convention, not a theorem. |
| Detection onset vs `L_pred` | **Downgraded to a refinement.** Bombieri §13 already has the critical window `t_c`; my contribution is tying it to `T*(L)` and measuring the basis dependence. |
| Hurwitz / witness discussion | **Not new.** CvS Step 5 and Remark 2.3. Present as making their point concrete. |

**Nothing in this project should be described as new without the qualifications above.**
