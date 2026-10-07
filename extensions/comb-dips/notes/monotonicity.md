# H is strictly increasing — used as an exact bound in the branch-and-bound (T1)

`H(t) = Re psi(1/4 + i t/2) - log pi`.

```
H'(t) = d/dt Re psi(1/4 + it/2) = Re[ psi'(1/4 + it/2) * (i/2) ] = -(1/2) Im psi'(1/4+it/2).
```

With `z = 1/4 + it/2`, `t > 0`, the series `psi'(z) = sum_{k>=0} (z+k)^{-2}` gives

```
Im (z+k)^{-2} = Im( (conj z + k)^2 ) / |z+k|^4 = -2 (1/4 + k)(t/2) / |z+k|^4  <  0
```

for every `k >= 0`, so `Im psi'(z) < 0` and hence **`H'(t) > 0` for all `t > 0`**.

**Consequence used in Stage 2.** On any cell `[a, b]`, `min H = H(a)` *exactly* — no
Lipschitz slack is needed for the smooth trend. The rigorous cell bound is therefore

```
min_{[a,b]} Psi_L  >=  H(a) - min( A_L ,  P_L(m) + Lip_P * (b-a)/2 ),
m = (a+b)/2,  Lip_P = sum |w_k| * log n_k ,
```

using also the unconditional `P_L <= A_L`, which is what makes the envelope cut free.
