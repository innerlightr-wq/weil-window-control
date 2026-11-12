(The grouping matters: the terms are `O(1/m^2)`. The variant
`sum_{m>=0}[C(0)/(m+1) - 2 int C e^{-(m+1/2)u}]` printed in earlier revisions of
this file **diverges**, its `m`-th term behaving like `-C(0)/m`; the constant
there also dropped the `log 4` of `int_0^inf (1-e^{-s})/sinh s ds = 2 log 2`.
Nothing computed ever used that expression: the pipeline uses the frequency form
below, via `baseline.H`, and the time-domain evaluator uses
`zeta_window.arch_closed`.)

