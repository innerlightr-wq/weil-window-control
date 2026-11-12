> # **CERTIFIED: `1.4330e-17 <= lambda_Q(4/5) <= 1.674e-17`**
>
> **Updated 2026-11-12 (`REPORT.md` §§12-13).** With **no reassembly** — the
> frozen `T=400` blocks had more in them — the lower endpoint rose from
> `1.41e-17` to `1.43300048827856e-17`, and the odd-sector constant from
> `2.3e-17` to `1.3857e-14` (`577x`; now `1.689x` Zhu's `8.206e-15`). The
> `T=400` even route is exhausted: certified floor and measured ceiling agree
> to `5.0e-5`.
>
> Separately certified, **per parity sector**: conditioning on four explicit
> even directions gives `Q_zeta >= 0.617553710936882 ||f||^2`, and on three odd
> directions `>= 0.106805801391493` (the paper's transitional eigenvalue),
> rising to `0.784130096434765` with a fourth. Over **both** sectors, `7` conditions give `>= 0.1068` and `8` give
> `>= 0.6176`. Four conditions on the combined space cannot: min-max forbids it,
> since the odd sector carries three minuscule eigenvalues of its own. This does
> **not** raise `lambda_Q(4/5)`, which stays in the bracket above.
>
> Counts below `10^{-2}`: **exactly four** even (both bounds certified) and
> **at most three** odd, against the paper's Table 3 entries `4` and `3`.
>
> A second window is now **bracketed on both sides**:
> `8.27343868254787e-18 <= eps <= 9.475e-18` at
> `lambda_-^2 = 4.999999999999999927` (§§13-14). By monotonicity the upper
> bound holds for every `lambda >= lambda_-`, including `lambda^2 = 5` exactly,
> where the paper's `9.1663e-18` falls **inside** our bracket.
>
> And the normalization dictionary to arXiv:2611.13129v1 is **proved**
> (`lambda = e^L`, verified to ~30 digits), giving a certified
> `eps_lambda(1)/delta_4 in [6.3999, 7.4723]` at `L = 4/5` — a **numerical
> comparison** until `delta_4` is interval-enclosed, which it is **not yet**. New drivers must
> go through `phase5/src/constants.py`: see `RETRACTIONS.md` **R15**.
>
> Superseded headline:
