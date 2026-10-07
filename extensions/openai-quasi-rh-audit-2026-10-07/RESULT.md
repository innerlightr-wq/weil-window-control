# RESULT.md

## Verdict: **NO BRIDGE FOUND**, with one **required correction** to my own manuscript

Two separate outputs. The second is the answer to the research question; the first is a defect
found on the way and must be fixed regardless.

---

## A. Required correction (my manuscript, not the source)

**Theorem 7's closing sentence is not justified as written.** Verbatim (`paper/weil_window_control.tex:294–296`):

> "For $L>0.8$ it holds with $\lambda^*(L)$ replaced by $0$, using no positivity input at all."

Substituting $0$ for $\lambda^*(L)$ in a lower bound requires $\lambda^*(L)\ge0$, which *is*
finite-window Weil positivity at $L$ — unknown for $L>0.8$. At $\delta=0$ the substituted
statement reads $\lambda_{\min}(Q_0)\ge0$: **it asserts the input it disclaims.** The same step is
the last link of the chain in `notes/proofs.md` §B.4.

What is unconditional is the **relative** bound:
$$\lambda_{\min}(Q_\delta)\ \ge\ \lambda^*(L)-c(L)\delta^2-\tfrac{16}{3}L^5\delta^4e^{2L\delta},\qquad c(L)=8L^3\bigl(\tfrac13+\tfrac1{\sqrt5}\bigr),$$
i.e. *the perturbation is at most $c(L)\delta^2+$ remainder, whatever the sign of the baseline.*
Suggested replacement text is in `BASELINE_AUDIT.md` §5. **Nothing else is affected**: Theorem 6,
Theorem 7 for $L\le0.8$, Theorem 8, Corollary 9, the dictionary and every numeric stand, and
`Limitations` item 1 already states the honest position.

**Everything else I audited is correct.** Eq. (4) re-derived from scratch with the factor 4
confirmed; the baseline shift is exactly $Q_0=Q_{\mathrm{true}}+2F(\gamma)^2$; the remainder bound
is a convergent series tail (hence valid at *every* $\delta$, not an asymptotic); $K(L,\gamma)$,
its closed form, its extremal function, Corollary 9 and both printed tables verified to
$10^{-15}$. **I withdraw an intermediate claim of mine that eq. (5) had two sign errors** — that
was an artifact of `pdftotext` rendering `\left[ \right]` as stray characters; the TeX and the code
are right (`BASELINE_AUDIT.md` §4).

## B. The research question: no transfer

**Strongest statement actually proved here.** None about ζ. The only new rigorous content is
negative and ledger-shaped:

> **Observation.** Let $f$ be real, even, unit-norm, supported in $[-L,L]$. For a zero at
> $\tfrac12+\delta+i\gamma$, Cauchy–Schwarz gives
> $|F(\gamma+i\delta)|^2\le\sinh(2L\delta)/\delta$, so each off-line quartet depresses the window
> form by at most $4\sinh(2L\delta)/\delta$. OpenAI Theorem 1.1
> (`paper.tex:106`, pinned `adc7f1241b42`) caps $|\delta|\le3/8$, hence caps this at
> $4\sinh(3L/4)/(3/8)$ in place of the trivial $8\sinh(L)$ — a gain of 1.04 at $L=1$, 1.41 at
> $L=3$, 15.0 at $L=12$. **This bound is uniform in $\gamma$ and therefore not summable over
> zeros; it yields no bound on $\lambda^*(L)$.**

**Why the route fails, localized.** What my Theorem 7 is missing for $L>0.8$ is a one-sided bound
$Q_{\mathrm{rest}}(f)\ge-\epsilon(L)\|f\|_2^2$ on the off-line zeros' contribution. That requires
control of the **number and band-weight** of off-line zeros. Theorem 1.1 controls their **depth**
and nothing else, and the manuscript itself separates the two (`paper.tex:100–104`: zero-density
estimates "address a different question … do not exclude every zero"). The needed currency is
exactly what the source states it does not produce.

Three further reasons the gap is not a matter of tuning:

1. **No exponent helps.** For any $\theta\in(\tfrac12,1)$, $\Re s>\theta$ leaves every zero with
   $|\delta_\rho|<\theta-\tfrac12$ free, and those carry the mass: $\lambda^*(L)$ falls
   super-exponentially (measured $\lambda^*(0.8)\approx2.3\times10^{-17}$), so one quartet of depth
   $\approx1.8\times10^{-9}$ already exhausts the budget at $L=0.8$. Only $\theta=\tfrac12$ — RH —
   removes them.
2. **The logical direction is reversed.** The geometric side is archimedean + pole $-\log\pi$ minus
   a *finite* von Mangoldt sum over $n\le e^{2L}$ and contains **no zero data**. Window positivity
   is an unconditional statement about primes that *constrains* zeros; a theorem about zero
   location cannot supply it.
3. **No shared object (Route B).** Their estimates are power savings for mean squares of
   sextic-twisted ideal Möbius sums over $\mathbb Q(\sqrt{-3})$ (Poisson, cubic-theta automorphy,
   quadratic/sextic large sieves). Nothing in that proof satisfies the hypotheses of my
   Toeplitz/Carathéodory–Fejér lemmas, and it contains no moment-positivity step to sharpen.

**Secondary, genuinely useful but elementary.** At the cap $\delta=3/8$, Theorem 7's own
$O(\delta^4)$ allowance overtakes its $\delta^2$ main term at $L=1.5897$, and the elementary bound
$4\sinh(2L\delta)/\delta+8L$ is sharper for $L>2.2848$. So the right bound to quote in the regime
that $7/8$ makes relevant is the minimum of the two. This is Cauchy–Schwarz, not new mathematics;
its only content is the crossover, and it does not advance RH. I do **not** label it a theorem.

## C. Prior-art comparison and what is mine

| | Whose |
|---|---|
| $7/8$ half-plane, the continuation criterion, cubic-theta/large-sieve machinery | **OpenAI**, `adc7f1241b42`. Not verified by me beyond the dependency chain; **no Lean build was run** |
| Window Toeplitz form, step (i), Hurwitz step, parity split | Hallouin–Perret; Connes–van Suijlekom |
| Certified $\lambda^*(0.8)\in[8.9\times10^{-18},2.27\times10^{-17}]$, bandwidth $T^*(L)=2\pi e^{2L}$, Landau–Widom decay | **Zhu** |
| Planted-ζ precedent, negative-eigenvalue count | Bombieri; window form itself Yoshida |
| $\delta^2$–$L^3$ law, Theorems 6–8, the dictionary | **mine** (prior to this audit) |
| The $\lambda^*\to0$ circularity; re-derivation of eq. (4); verification of $K$ and both tables; the $\delta=3/8$ cost ledger and its two crossovers; the depth-versus-count localization | **mine, this audit** |

**Novelty check.** I claim no new theorem, so no primary-prior-art search for a new result was
needed. The two substantive points are an internal audit finding and an elementary ledger.

## D. Remaining gap

A lower bound on $\lambda^*(L)$ for $L>0.8$. The estimate that would close it:

> for every $f$ real, even, unit-norm, supported in $[-L,L]$,
> $\displaystyle\sum_{\rho:\,\Re\rho\neq1/2}|F(\gamma_\rho+i\delta_\rho)|^2\le\epsilon(L)$
> with $\epsilon(L)<\lambda^*_{\text{certified}}$ — a **window-band-weighted zero-density
> estimate valid down to the critical line**.

Every known zero-density estimate degenerates to the full zero count as $\sigma\to\tfrac12$, which
is why this is hard and why no fixed zero-free half-plane substitutes for it. Nothing in family 003
is of this type.

Realistic alternatives, not attempted here: extend Zhu's certification to $L>0.8$ (a numerical
analysis problem about an explicit finite-data operator, needing no new arithmetic); or keep
Theorem 8's conditional form and state it as conditional.

---

## Plain-language summary

I could not find a local copy of the project, so I cloned it fresh and left it untouched; the PDF
in Downloads turned out to be the committed PDF re-saved by Preview, with identical text, so there
was no version ambiguity to worry about. OpenAI's repository has a single commit, which I pinned.
Their 7/8 result is a real and clearly written theorem, and their Lean scope page is honest, but
the file their catalogue points a reader to first is a statement stub ending in `sorry` — the
actual proof lives elsewhere, and I did not build it, so I am claiming nothing about formal
verification.

On the mathematics: their theorem says every zero of ζ sits in the strip
$\tfrac18\le\Re s\le\tfrac78$. My window result needs something different in kind — not *how far*
a stray zero can be from the critical line, but *how much total damage all stray zeros can do*.
A depth cap divided by nothing is still infinite: their bound on a single zero's effect does not
shrink as the zero climbs, so adding up over infinitely many possible stray zeros gives nothing.
Worse, the margin I would be protecting is about $10^{-17}$, and a single stray zero only
$2\times10^{-9}$ off the line would eat it. That stays true no matter how much the 7/8 is improved,
right up to — but not including — the Riemann hypothesis itself. And the direction of information
runs backwards for this purpose: the quantity I am missing is built from primes, not zeros, so it
is the sort of thing that constrains zeros rather than being constrained by them.

The audit did earn its keep in one place. My Theorem 7 says that for windows larger than 0.8 it
holds "with $\lambda^*(L)$ replaced by $0$, using no positivity input at all". That is not right:
setting that term to zero assumes the window is positive, and at $\delta=0$ the sentence asserts
precisely the fact it says it is not using. The fix is small — state it as a relative bound — and
nothing else in the paper moves. Separately, I talked myself into believing one of my own formulas
had two sign errors, then found the error was in how I was reading the PDF rather than in the
paper; I have recorded that retraction rather than quietly dropping it.

**Verdict: NO BRIDGE FOUND** — these particular estimates do not transfer usefully — together with
one required correction to the manuscript's Theorem 7.
