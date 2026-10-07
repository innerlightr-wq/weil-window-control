# BASELINE_AUDIT.md — my own statements, re-derived before any transfer

All line numbers are `paper/weil_window_control.tex` at HEAD `f6033d7`. Formulas were read from
the **TeX source**, not from `pdftotext` (see the retraction in §6 — this matters).

## 1. Conventions, recovered exactly

| Item | Value in the project |
|---|---|
| Fourier convention | $F=\hat f$ with $F(\xi)=\int f(u)e^{i\xi u}\,du$ (no $2\pi$); fixed by `src/zeta_window.py:F` and consistent with $K$ below |
| Form domain | $f$ real, **even**, $\operatorname{supp}f\subset[-L,L]$, $\|f\|_2=1$; §4.2 preamble (tex 263–264) |
| Parity sector | even sector throughout §4.2; the odd sector's pole sign flip $-2P_jP_k$ is derived separately in §5 |
| Window form | $M=\text{Pole}+\text{Arch}-\log\pi\cdot I-\text{Prime}$ (`src/zeta_window.py:214`) |
| Prime term | `von_mangoldt_terms(2*L)` — a **finite** von Mangoldt sum over $n\le e^{2L}$ |
| Zero multiplicities | ζ zeros assumed simple (true for every computed zero); the counterfactual is count-preserving against a **hypothetical double** on-line zero |
| $\lambda^*(L)$ | used as $\lambda_{\min}$ of the baseline window form; **only ever introduced as a macro (`\newcommand` tex 25), never formally defined.** See §5 |

**Structural fact, load-bearing below.** The geometric side is archimedean + pole + $-\log\pi$
minus a *finite* von Mangoldt sum. **It contains no zero data at all.** So $\lambda^*(L)\ge0$ is a
positivity question about an explicit, finite-data operator — which is why Zhu can certify it at
$L=0.8$ by rigorous numerics.

## 2. The count-preserving counterfactual — VERIFIED independently

Baseline (tex 272, eq. `baseline`): $Q_0=Q_{\mathrm{true}}+2F(\gamma)^2$. **Recorded as the task
asks: the baseline sits $+2F(\gamma)^2$ above the actual zeta Weil form.**

I re-derived eq. `pert` (tex 285) from scratch. With $\rho=\tfrac12+iγ_\rho$, the quartet
$\{\tfrac12\pm\delta\pm i\gamma\}$ has $γ_\rho\in\{\gamma\mp i\delta,\,-\gamma\mp i\delta\}$;
since $F$ is even, the $\pm\gamma$ pairs coincide, so the quartet contributes
$2\bigl[F(\gamma-i\delta)^2+F(\gamma+i\delta)^2\bigr]=4\operatorname{Re}F(\gamma+i\delta)^2$.
With $G(\delta)=F(\gamma+i\delta)=F+i\delta F'-\tfrac{\delta^2}{2}F''+O(\delta^3)$,
$\operatorname{Re}G^2=F^2-\delta^2(F'^2+FF'')+O(\delta^4)$. The double on-line baseline
contributes $4F(\gamma)^2$. Hence

$$Q_\delta(f)-Q_0(f)=-4\delta^2\bigl(F'(\gamma)^2+F(\gamma)F''(\gamma)\bigr)+R_4 .$$

**Matches the manuscript exactly, including the factor 4.** The three distinct operations the task
asks me to separate are kept separate in the manuscript and here:
*adding* zeros (changes the count), *replacing a pair by a quartet* (changes the count),
*preserving multiplicity* (the double-zero baseline, the only convention with
$Q_\delta\to Q_0$ as $\delta\to0$ — tex 279–281). The cost of that continuity is exactly the
$+2F(\gamma)^2$ baseline shift, which the manuscript states rather than hides.

## 3. The remainder and its $L|\delta|$ dependence — VERIFIED, and it is *not* an asymptotic expansion

`notes/proofs.md` §B.3: with $|m_k|\le\sqrt{2L}\,L^k$ one gets $G(\delta)^2=\sum c_n\delta^n$,
$|c_n|\le 2L(2L)^n/n!$; $\operatorname{Re}G^2$ is **even** in $\delta$ so odd terms vanish, and

$$|R_4|\;\le\;8L\sum_{n\ge4}\frac{(2L\delta)^n}{n!}\;\le\;\tfrac{16}{3}L^5\delta^4e^{2L\delta}.$$

I confirm the arithmetic: $8L\cdot\frac{(2L\delta)^4}{4!}e^{2L\delta}=\frac{128}{24}L^5\delta^4e^{2L\delta}=\frac{16}{3}L^5\delta^4e^{2L\delta}$.

**Important correction to how this is usually read.** This is a convergent **series-tail** bound,
so Theorem 7 is *rigorously valid at every* $\delta$, including $\delta=3/8$ — it is not a
small-$\delta$ asymptotic that becomes false. What fails at large $L\delta$ is its *usefulness*:
the $\delta^4$ allowance overtakes the $\delta^2$ main term. Measured (`bridge_cost.py`):

> at $\delta=3/8$ the remainder equals the main term at $L=1.5897$, and beyond that **Theorem 7's
> bound is dominated by its own $O(\delta^4)$ allowance.**

So the task's warning is justified in substance — one must not *use* the $\delta^2$ law at
$\delta=3/8$ for $L\gtrsim1.6$ — but the reason is loss of content, not loss of validity.

## 4. $K(L,\gamma)$ — VERIFIED, with my own false alarm retracted

$K(L,\gamma)=\sup\{|F'(\gamma)|^2\}$ over the unit ball of the form domain. With the convention of
§1, $F'(\gamma)=-\int_{-L}^{L}u\sin(\gamma u)f(u)\,du$, and $u\mapsto u\sin(\gamma u)$ is **even**,
so it lies in the domain; Cauchy–Schwarz gives
$K(L,\gamma)=\|u\sin(\gamma u)\|_2^2=\int_{-L}^{L}u^2\sin^2(\gamma u)\,du$, attained at
$f\propto u\sin(\gamma u)$. **Exactly as stated** (tex 308–313).

The closed form (tex 315–316) is $K=\frac{L^3}{3}-\Bigl[\frac{L^2\sin2\gamma L}{2\gamma}+\frac{L\cos2\gamma L}{2\gamma^2}-\frac{\sin2\gamma L}{4\gamma^3}\Bigr]$.
Checked four ways (`check_K.py`): it satisfies the defining ODE $dK/dL=2L^2\sin^2(\gamma L)$ to
10 digits at six values of $L$; it agrees with quadrature to $10^{-15}$; and it reproduces the
manuscript's own tables of $4K(L,\gamma_1)$ (0.7419, 1.3427, 3.1158, 5.1130, 10.652) and of
$|K/(L^3/3)-1|$ (17.3%, 8.68%, 0.699%, 6.37%, 6.38%, 0.136%) **exactly**. **CORRECT.**

> **Retraction of my own intermediate finding.** I first reported two sign errors in this closed
> form. That was wrong. It was an artifact of reading the `pdftotext` extraction, where the tall
> `\left[ \right]` delimiters are emitted as stray `"` and `#` characters on a separate line, so
> the bracket vanishes and two signs appear flipped. The TeX, the code
> (`src/task_final_checks.py:11`, which brackets the group explicitly) and every printed number
> agree. **No defect. Formulas in this audit were re-read from TeX after this.**

**The $K$-versus-$C(L)$ distinction the task asks about is already correctly drawn in the
manuscript.** $4K$ is the sharp *norm* of the derivative-evaluation functional and hence a valid
$\delta^2$ *coefficient in a lower bound* (Remark 10); it is **not** claimed to be the coefficient
of the re-optimised minimum eigenvalue. §4.3 measures the gap honestly: the re-optimised minimiser
reaches 41.9% of $4K$ at $L=0.8$ and 96.3% at $L=1.6$. `notes/proofs.md` §B.5 states explicitly
that $C(L)$ is a property of the *re-optimised* minimiser and that
$4(F'^2+FF'')$ at the *unperturbed* ground state is 3–6 orders of magnitude too small. Attaining
the derivative norm does not minimise the whole form, and the manuscript does not say it does.

## 5. BASELINE ISSUE — one sentence is not justified as written

**Theorem 7, final sentence (tex 294–296), verbatim:**

> "This is unconditional exactly where Weil positivity on the window has been certified: for
> $L\le0.8$, since Zhu certifies $\lambda^*(0.8)\ge8.9\times10^{-18}$ and $\lambda^*$ is
> non-increasing in $L$. **For $L>0.8$ it holds with $\lambda^*(L)$ replaced by $0$, using no
> positivity input at all.**"

The same step appears in `notes/proofs.md` §B.4, whose Theorem B chain ends
`… ≥ λ*(L) − c(L)δ² − (16/3)L⁵δ⁴e^{2Lδ} ≥ −c(L)δ² − (16/3)L⁵δ⁴e^{2Lδ}`.

**The final inequality in that chain requires $\lambda^*(L)\ge0$.** That is precisely
finite-window Weil positivity at $L$, which is *unknown* for $L>0.8$ — it is the entire subject of
the paper, and as $L\to\infty$ it is RH.

**Applying the task's test.** Set $\delta=0$. The substituted statement reads
$\lambda_{\min}(Q_0)\ge0$ — i.e. **at $\delta=0$ the conclusion already implies the positivity
input it claims not to use.** The claim is circular, not merely loose.

**What is genuinely unconditional** is the *relative* bound, which needs no sign for $\lambda^*$:

$$\boxed{\;\lambda_{\min}(Q_\delta)\;\ge\;\lambda^*(L)-c(L)\delta^2-\tfrac{16}{3}L^5\delta^4e^{2L\delta}\;}$$

equivalently $\lambda_{\min}(Q_\delta)-\lambda^*(L)\ge-c(L)\delta^2-\tfrac{16}{3}L^5\delta^4e^{2L\delta}$:
*the perturbation is at most $c(L)\delta^2+$ remainder, whatever the sign of the baseline.*

**Recommended correction** (not applied — the task forbids editing the manuscript): delete "using
no positivity input at all" and replace the sentence with

> "For $L>0.8$ the inequality continues to hold as stated, but $\lambda^*(L)$ is then an unknown
> quantity of unknown sign; what remains unconditional is the relative bound
> $\lambda_{\min}(Q_\delta)-\lambda^*(L)\ge-c(L)\delta^2-\frac{16}{3}L^5\delta^4e^{2L\delta}$.
> Replacing $\lambda^*(L)$ by $0$ would assume window positivity at $L$."

**Blast radius — deliberately small.** This does **not** affect: Theorem 6 (function field,
self-contained); Theorem 7 for $L\le0.8$ (Zhu supplies $\lambda^*>0$); Theorem 8 (explicitly
conditional on all other zeros being on the line); Corollary 9; the dictionary; §§3, 5, 6, 7; or
any numeric. It affects one clause about $L>0.8$, and `Limitations` item 1 (tex 557) already says
the honest thing — *"No rigorous lower bound on $\lambda^*(L)$ of our own"* — so the paper
contradicts itself in the reader's favour elsewhere.

## 6. Theorem 8's hypothesis — stated exactly, and it is strong

$\widetilde Q:=Q_0-4F(\gamma)^2\ge0$, which via the baseline shift reads
$Q_{\mathrm{true}}-2F(\gamma)^2\ge0$, i.e. **"every nontrivial zero other than the one being moved
lies on the critical line"** (tex 276–279). The manuscript says this in those words. It is
RH-minus-one-zero: far stronger than any quasi-RH statement, and in particular **not implied by a
zero-free half-plane at any $\theta<1$.** This is the hypothesis a transfer would have to replace.

## 7. Four things kept logically distinct, as the task requires

| | In this project |
|---|---|
| (a) bound on a perturbation | Theorems 7, 8, Corollary 9 — all of them. This is all the $\delta^2$–$L^3$ law is. |
| (b) negative trial direction for a fabricated spectrum | §3.2 planted spectra, §6 detection onset. Explicitly fabricated; Non-claim 3 says so. |
| (c) positivity of the actual Weil form | **Not proved here.** Zhu's certificate at $L=0.8$ only; Non-claim 5 and Limitation 1 say so. |
| (d) exclusion of an actual zeta zero | **Nowhere.** Non-claim 1: *"No proof or advance of RH for ζ."* |

The manuscript does not conflate these. The audited defect in §5 is the one place where (c) is
quietly borrowed.
