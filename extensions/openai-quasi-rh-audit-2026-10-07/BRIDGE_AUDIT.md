# BRIDGE_AUDIT.md — dependency chain, proposed translation, cost ledger, controls

All OpenAI citations are to `openai/math` at **`adc7f1241b42`** with file and line.

## 1. Dependency chain ending at the zero-free half-plane

### 1a. The alternate manuscript (`paper2.tex`, Oct 5) — the 11/12 route

Reconstructed from §§1–3 and `sec:reduction`:

1. **Target object.** `paper2.tex:215–226`: for a fixed finite-order Hecke character $\nu$ of
   $K=\mathbb Q(\sqrt{-3})$, the **smoothed $\nu$-twisted ideal Möbius sum**
   $A_1(D)=\sum_{\mathfrak n}\mu(\mathfrak n)\nu(\mathfrak n)W(\mathrm N(\mathfrak n)/D)$,
   $W\in C_c^\infty((0,\infty))$, ideals prime to $2,3$ and the conductor.
2. **Reduction.** `paper2.tex:229–239`: a power saving
   $A_1(D)\ll_{\nu,W,\varepsilon}D^{1-\delta+\varepsilon}$ implies a zero-free half-plane — the
   *smoothed Hecke version of Littlewood's* $\sum_{n\le x}\mu(n)\ll x^{1/2+\varepsilon}\iff$ RH
   relation (they cite Lit12, MM09 §1).
3. **Family embedding.** `paper2.tex:266–268`: $A_u(D)=\sum_{\mathfrak n}\mu\,\nu\,\chi_{\mathfrak n}(u)W(\cdot)$
   with $\chi$ sextic residue characters; $u$ = **row**, $\mathfrak n$ = **column**.
4. **The crucial estimate** (`Proposition thm:ms`, display `eq:intro-ms`, `paper2.tex:272–276`):
   $$\sum_{0<\mathrm N(u)\le H}|A_u(D)|^2\ \ll_{\nu,W,\vartheta,\varepsilon}\ D^{1+\varepsilon}H,\qquad H=D^{1+\vartheta},\ 0<\vartheta\le1/10 .$$
   Square-root cancellation on average over the row.
5. **Mean-square extraction** (`paper2.tex:284–298`): $\chi_{\mathfrak n}(p^6)=\mathbf 1_{p\nmid\mathfrak n}$
   gives $A_{p^6}(D)=A_1(D)+O_W(D/Y)$; Landau's prime ideal theorem supplies $\asymp Y/\log Y$
   primary primes with $Y/2<\mathrm N(p)\le Y$; taking $Y=H^{1/6}$,
   $$|A_1(D)|^2\ \ll\ D^{1+\varepsilon}H^{5/6}+D^2H^{-1/3},$$
   and $H=D^{1+\vartheta}$ yields $A_1(D)\ll D^{11/12+5\vartheta/12+\varepsilon}$, hence $11/12$.
6. **Proof of the mean square** (`paper2.tex:303`): Poisson in the row $u$ → sextic Gauss sums →
   (Gauss–Jacobi) absorb $\mu$ → **cubic** Gauss sums → coefficients of **Kubota's cubic theta
   function** → automorphy/theta transformation (§`sec:reflection`). The decisive sieve input is
   **Goldmakher–Louvel's quadratic large sieve** (`paper2.tex:412`), which gives the factor
   $M+L$ and *avoids* the extra $(ML)^{2/3}$ of the general higher-order sieve of
   Blomer–Goldmakher–Louvel.

**Where 11/12 comes from, exactly.** $11/12=\tfrac12(1+\tfrac56)$: the $\tfrac56$ is
$1-\tfrac16$ with the $\tfrac16$ forced by $Y=H^{1/6}$ — and the sixth power is forced by the
*sextic* character structure ($\chi_{\mathfrak n}(p^6)$ is the trivialising power), balanced
against the second term $D^2H^{-1/3}$. **It is not an unoptimised parameter**; it is the balance
point of the extraction with a sixth-power trivialiser.

### 1b. The stronger manuscript (`paper.tex`, Sep 30) — the 7/8 route

`paper.tex` is a single 16,677-line document containing **both** parts.

1. **Quantity continued.** `paper.tex:175–179`: $\beta_*$ := the supremum of $1/2$ *and* the real
   parts of zeros in $1/2\le\Re s\le1$ for the family of primitive finite-order Hecke characters
   over $F$. **Note $\beta_*\ge1/2$ by definition.**
2. **Common continuation criterion.** `paper.tex:181–186`: for each target $\eta$ and large real
   scale $Z$, compare a normalized character sum against a Mellin integral containing
   $1/L_F(s,\eta)$, after deleting finitely many Euler factors and multiplying by a holomorphic
   factor bounded away from zero. A direct bound for the sum **plus** a power-saving bound for its
   difference from that integral gives the continuation.
3. **The two parts differ in the affine power** (`paper.tex:187–192`):
   $C_{\mathrm I}(s)=s-\tfrac23$ (Part I → 11/12) and $C_{\mathrm{II}}(s)=s-\tfrac{11}{16}$
   (Part II → 7/8). *"Thus the common analytic principle does not identify the two normalized
   sums."* Only the positive power margins must be independent of $\eta$.
4. **Completed sum.** `paper.tex:198–209`: smoothed cubic-theta Fourier coefficients averaged
   against sextic residue characters and target $\eta$; completed indices $cn^3$ ($c$ squarefree,
   $c,n\equiv1\bmod3$); the character acts on the whole product.
5. **Two exact representations.** `prop:completed-reflection` (completed cubic reflection) and
   Poisson summation, whose nonzero frequencies are $ua^6$ with $u$ sixth-power-free
   (`paper.tex:210–219`). For $u=1$ the row quotient **contains the reciprocal of the target
   $L$-function**.
6. **Estimates.** Quadratic large sieve (reflected-row mean square) + planar additive large sieve
   (reduced fractions in $\mathbb C/\mathcal O$) + Cauchy–Schwarz (`prop:balanced-low`); a **zero
   detector with saturated witnesses** (`sec:detector`) assigning zero-free rectangles; and a
   **sextic large sieve** specialization (`lem:sextic-large-sieve`) for the row count
   (`prop:sextic-row-count`).
7. **Principal row** residues + the local Euler identity give a nonzero scalar multiple of the
   target Mellin integral; normalization + the two bounds verify the criterion.
8. **Part II adds** prime compensation (selected primes cancel an unwanted Euler contribution),
   asymmetric scales, and two further moment estimates (`paper.tex:258–262`, `sec:compensation`).

**Decisive gain:** a power saving in $Z$ with exponent margins chosen independently of the target
character. **Conversion to nonvanishing:** the continuation criterion — if $\beta_*$ exceeded the
proposed boundary, $1/L_F(s,\eta)$ would continue holomorphically across a common positive
distance to the left of $\beta_*$, contradicting the definition of $\beta_*$. **This is the
"nonzero pole forced into a holomorphic region" contradiction, not a positivity argument.**

**The two mechanisms are not the same.** Both end at a reciprocal-$L$/Möbius power saving, but
`paper2` extracts $A_1$ from a row mean square via sixth powers, whereas `paper.tex` Part I runs
the *same analytic principle on a different normalized sum* with $C_{\mathrm I}(s)=s-2/3$, and
Part II changes the base sum itself. The manuscript says so explicitly at line 190.

## 2. Transferable-estimate record

| Field | `eq:intro-ms` (the one candidate worth recording) |
|---|---|
| Exact statement | $\sum_{0<\mathrm N(u)\le H}|A_u(D)|^2\ll D^{1+\varepsilon}H$, $H=D^{1+\vartheta}$ |
| Source | `paper2.tex:272–276`, Proposition `thm:ms` |
| Scale | norm scales $D$ (column), $H=D^{1+\vartheta}$ (row), $\vartheta\in(0,1/10]$ |
| Character / conductor | sextic residue characters over $\mathcal O_{\mathbb Q(\sqrt{-3})}$, twisted by a **fixed** finite-order Hecke $\nu$; ideals prime to $2,3,\mathrm{cond}(\nu)$ |
| Height range | none — this is a *norm*-scale estimate, not an estimate in the analytic conductor/height $t$ |
| Weight hypotheses | $W\in C_c^\infty((0,\infty);\mathbb C)$, arbitrary; constants depend on $W$ |
| Effective constants? | **No.** $\ll_{\nu,W,\vartheta,\varepsilon}$; `003.md` separately states no explicit $c$ for the Siegel corollary |
| Fixed vs uniform | **fixed-character**: constants and lower thresholds may depend on $\eta$/$\nu$; only the *power margins* are uniform (`paper.tex:191–192`) |
| Pointwise vs averaged | **averaged** (mean square over the row $u$); pointwise information about $A_1$ is recovered only through the sixth-power extraction, at the cost $H^{5/6}$ |
| Exponent losses | $D^\varepsilon$ throughout; $H^{5/6}$ and $H^{-1/3}$ in the extraction; $5\vartheta/12$ in the final exponent; deletion of finitely many Euler factors; $Y=H^{1/6}$ |
| Auxiliary normalizations | primary-generator convention; characters extended by zero on nonunits; holomorphic factor bounded away from zero |

**Every field above is in the wrong currency for the window Weil form**, for the reason in §4.

## 3. The proposed translation, written out, and where it breaks

The template from the Sturmian work is *exact correspondence → quantitative gain → full cost
accounting → independent contradiction*. Its native substitutes here:

| Sturmian part | Native RH substitute attempted |
|---|---|
| exact correspondence | the explicit formula: geometric side $\equiv$ zero side on the window |
| quantitative gain | $|\delta_\rho|\le3/8$ from Theorem 1.1 instead of the trivial $|\delta_\rho|\le1/2$ |
| full cost accounting | §4 ledger below |
| independent contradiction | Zhu's certified $\lambda^*(0.8)\ge8.9\times10^{-18}$ — an unconditional lower bound on the *actual* form |

*No Sturmian object is transplanted. No zeta zero is assumed algebraic; Ridout and the Subspace
Theorem are not invoked anywhere.*

**Route A, stated precisely.** Define the remaining-zero contribution: for unit $f$ even with
$\operatorname{supp}f\subset[-L,L]$,
$$Q_{\mathrm{rest}}(f)\;=\;\sum_{\rho\ \text{off-line}}4\operatorname{Re}F(\gamma_\rho+i\delta_\rho)^2 ,$$
so that $Q_{\mathrm{true}}=\sum_{\rho\ \text{on-line}}F(\gamma_\rho)^2+Q_{\mathrm{rest}}$ with the
first sum $\ge0$. **A lower bound $Q_{\mathrm{rest}}(f)\ge-\epsilon(L)\|f\|_2^2$ would give
$\lambda^*(L)\ge-\epsilon(L)$ and discharge exactly the hypothesis missing in
`BASELINE_AUDIT.md` §5.** A one-sided bound suffices; an absolute bound is not needed.

**What Theorem 1.1 contributes to it.** By Cauchy–Schwarz on $[-L,L]$,
$$|F(\gamma+i\delta)|^2=\Bigl|\int f(u)e^{i\gamma u}e^{-\delta u}du\Bigr|^2\le\int_{-L}^{L}e^{-2\delta u}du=\frac{\sinh(2L\delta)}{\delta},$$
so each off-line quartet satisfies
$4\operatorname{Re}F^2\ge-4\sinh(2L\delta_\rho)/\delta_\rho$, and $|\delta_\rho|\le3/8$ caps that
at $4\sinh(3L/4)/(3/8)$ instead of $8\sinh(L)$. **That is the entire contribution: a cap on the
depth of each individual quartet.**

**Where it breaks.** $Q_{\mathrm{rest}}$ is a sum over *all* off-line zeros. The bound above is
**uniform in $\gamma$** — it does not decay with the ordinate — so summing it requires a bound on
the **number and band-weight** of off-line zeros:
$$\sum_{\rho\ \text{off-line}}|F(\gamma_\rho+i\delta_\rho)|^2\ \le\ \epsilon(L)\|f\|_2^2\quad\text{uniformly over the unit ball.}$$
**Theorem 1.1 bounds depth, not count.** It places every zero in $\Re s\in[1/8,7/8]$ and says
nothing about how many lie off the line or where. The manuscript draws this distinction itself,
`paper.tex:100–104`: *"Zero-density estimates address a different question: they bound the number
of zeros to the right of a given vertical line. … these zero-density bounds do not exclude every
zero."* The currency Route A needs is precisely the one the source says it is not providing.

**The obstruction is not about the exponent 7/8.** For any $\theta\in(1/2,1)$, a zero-free
half-plane $\Re s>\theta$ leaves every zero with $|\delta_\rho|<\theta-\tfrac12$ unconstrained,
and those near-line zeros are exactly where the dangerous mass sits: $\lambda^*(L)$ is
super-exponentially small (Zhu's Landau–Widom law; measured $\lambda^*(0.8)\approx2.3\times10^{-17}$),
so a single quartet of depth $\delta>\sqrt{\lambda^*/c(L)}\approx1.8\times10^{-9}$ at $L=0.8$
already exhausts the whole budget. Improving $7/8\to11/16\to\tfrac12+\varepsilon$ never reaches
$\varepsilon=0$, and only $\varepsilon=0$ — RH itself — removes the near-line zeros. **Scope: this
argument rules out the specific bridge "zero-free half-plane $\Rightarrow$ lower bound on
$\lambda^*(L)$". It does not show that no bridge from this machinery exists.**

**And the logical direction runs the other way.** The geometric side is archimedean + pole
$-\log\pi$ minus a *finite* von Mangoldt sum over $n\le e^{2L}$ (`BASELINE_AUDIT.md` §1): it
contains **no zero data**. So window positivity is an unconditional statement about primes that
*constrains* zeros — it is an input to zero-exclusion, not an output of it. A theorem about zero
location cannot supply it.

**Route B, stated precisely.** Can a lemma of mine improve an estimate in the OpenAI proof, or
build a better arithmetic witness? **No, and the reason is structural, not quantitative.** Their
estimates are power savings for mean squares of sextic-character-twisted ideal Möbius sums over
$\mathbb Q(\sqrt{-3})$, proved by Poisson summation, cubic-theta automorphy and large sieves. My
tools are positivity and Carathéodory–Fejér/Toeplitz facts about the finite-window Weil form
(Hallouin–Perret Prop. 5, Thm. 36(ii); the $\delta^2$–$L^3$ law). **There is no object in their
proof satisfying the hypotheses of my lemmas**: no Toeplitz Gram matrix of a normalized
Frobenius-type family, no moment-positivity step, no $\lambda_{\min}$ of a window form. Nor can I
write the contradiction: my Theorem 8 would deliver $Q(f)<0$ for a hypothetical off-line zero, but
pairing it against Zhu's certificate needs an **upper** bound on the contribution of all *other*
zeros, which is unavailable — which is exactly why Theorem 8 is conditional. **The implication
$4K(L,\gamma)\delta^2>\text{error}$ proves nothing on its own, and a lower bound
$\lambda_{\min}\ge-C\delta^2$ does not exhibit a negative direction.** Both traps are live here and
both are avoided.

## 4. Cost ledger (`bridge_cost.py`, reproducible)

At the cap $\delta=3/8$ supplied by Theorem 1.1, comparing the two valid bounds on the drop
$Q_\delta-Q_0$ — (P) Theorem 7's $c(L)\delta^2+\frac{16}{3}L^5\delta^4e^{2L\delta}$ and
(D) the elementary $4\sinh(2L\delta)/\delta+8L$:

| $L$ | (P) main $c(L)\delta^2$ | (P) remainder | (P) total | (D) total | sharper |
|---|---|---|---|---|---|
| 0.8 | 0.4496 | 0.0630 | 0.5126 | 13.19 | P |
| 1.5 | 2.9636 | 2.4670 | 5.4306 | 26.70 | P |
| 2.0 | 7.0249 | 15.126 | 22.151 | 38.71 | P |
| 2.5 | 13.721 | 67.162 | 80.883 | 53.96 | **D** |
| 3.0 | 23.709 | 243.16 | 266.87 | 74.04 | **D** |
| 6.0 | 189.67 | 73825 | 74015 | 528.0 | **D** |

- **Remainder overtakes main term at $L=1.5897$** ($\delta=3/8$).
- **(D) overtakes (P) at $L=2.2848$** ($\delta=3/8$).
- What $7/8$ buys over the trivial $|\delta|\le1/2$ in (D): factor 1.04 at $L=1$, 1.41 at $L=3$,
  5.41 at $L=8$, 15.0 at $L=12$ — growing, but multiplying an **uncontrolled sum**.

**Charged costs.** Conductor: none transferable (their estimate is in norm scale $D$, mine in
window length $L$ — no dictionary exists between them). Ordinate: their estimate carries **no**
height range; mine is governed by the band $T^*(L)=2\pi e^{2L}$. Test-function seminorms: their
$W$-dependence is unquantified ($\ll_W$); my $f$ must be real, even, unit-norm, supported in
$[-L,L]$. Support: their ideals are prime to $2,3,\mathrm{cond}\,\nu$; my prime sum is the *full*
von Mangoldt sum to $e^{2L}$ — deleting Euler factors is not available to me. Smoothing: theirs
smooth and compactly supported in norm; mine compactly supported in $u$. Contour shifts: theirs a
Mellin integral with affine $C(s)$; mine a shift $\gamma\mapsto\gamma+i\delta$, cost $e^{2L\delta}$,
**charged in full above**. Exceptional terms: their principal row is extracted separately, my pole
term is $2F(i/2)^2\ge0$ (even sector). Distance to the boundary: $3/8$, which is *not small* — the
cost of that is the whole of the table above.

## 5. Controls

- **Control that the ledger discriminates.** At small $\delta$, bound (D) does **not** vanish
  (it tends to $16L$ as $\delta\to0$) while the true drop tends to $0$. So (D) must lose for small
  $\delta$ — and it does, for all $L<2.2848$ at $\delta=3/8$ and for all $L$ as $\delta\to0$. A
  "bound" that won in every regime would indicate an error.
- **Control on $K$.** The defining ODE $dK/dL=2L^2\sin^2(\gamma L)$ is satisfied by the
  manuscript's closed form and violated by the bracket-dropped variant; this is what caught my own
  misreading.
- **Negative control on provenance.** `ComparatorChallenges/QuasiRiemannHypothesis.lean` is a
  `sorry` stub. Treating a comparator file as a proof would have produced a false "formally
  verified" claim; the check that distinguishes them is `formalization.yaml`, which maps the
  comparator to `OAI/NumberTheory/DirichletL/Nonvanishing.lean`.
- **Routes deliberately not rerun** (§7 of the task): the prime-comb resonance search
  (`extensions/comb-dips/NULL.md`, verdict NULL), naive orthogonal inheritance, and the scalar
  coupling-constant search. **No new analytic lemma changed their premises**, so their outcomes
  stand as recorded over their tested ranges and were not re-executed.
