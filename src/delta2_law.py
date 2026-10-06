#!/usr/bin/env python3
r"""Tasks 2 and 3: derive and check the delta^2 coefficient and the L^3 (Cauchy-Schwarz)
bound, on both the zeta side and the function-field side.

TASK 2 -- the expansion (T1).  For real even f, F = f-hat is real on R with real Taylor
coefficients there, so F(conj z) = conj F(z).  With w = gamma + i delta,

    F(gamma + i delta) = F + i d F' - (d^2/2) F'' + O(d^3)          (F, F', F'' at gamma)
    F(w)^2             = F^2 + 2 i d F F' - d^2 (F'^2 + F F'') + O(d^3)
    Re F(w)^2          = F^2 - d^2 (F'^2 + F F'') + O(d^4)          (odd powers are pure
                                                                     imaginary)
NORMALISATION (checked against the code's zeros-side convention in Gate A):
    on-line pair {+-gamma}            contributes   2 F(gamma)^2
    off-line quartet {+-gamma +- i d} contributes   4 Re F(w)^2
so the brief's "2 Re F(gamma - i delta)^2" is the PAIR-normalised half of the quartet, and
F(conj w) = conj F(w) makes the two signs of delta agree.  Moving a zero off the line
(remove the pair, insert the quartet) therefore changes the form by

    Delta = -2 F^2 + 4 Re F(w)^2 = 2 F(gamma)^2 - 4 d^2 (F'^2 + F F'') + O(d^4),

which is negative exactly when the trial function nearly annihilates gamma.  At the window
minimiser F(gamma_n) ~ 0 for the leading zeros, leaving  Delta ~ -4 d^2 F'(gamma_n)^2.

TASK 2c -- NOT the Laguerre expression.  The coefficient here is  F'^2 + F F'', whereas
the Laguerre/Laguerre-Polya quantity is  L_1 = F'^2 - F F''  (opposite sign on F F'').
They agree only to leading order where F(gamma) ~ 0 -- which is exactly what the minimiser
arranges, so the agreement is a property of the minimiser, not an identity.  Any link to
the Partition note's Remark 1 is therefore LEADING ORDER ONLY.

TASK 3 -- the L^3 bound (T1).  For f supported in [-L,L], real and even, with ||f||_2 = 1,
    F'(gamma) = -\int_{-L}^{L} u f(u) sin(gamma u) du,
so by Cauchy-Schwarz
    |F'(gamma)|^2 <= \int_{-L}^{L} u^2 sin^2(gamma u) du = L^3/3 - (1/2)\int u^2 cos(2 gamma u) du
                  -> L^3/3   as gamma -> infinity.
(The cruder bound that drops the sin weight gives 2L^3/3, a factor 2 worse.)  Hence
    C(L) := -Delta/delta^2 = 4 |F'|^2  <=  4L^3/3 + o(L^3).

FUNCTION-FIELD ANALOGUE.  For the off-line real pair {rho, 1/rho}, rho = 1+eps, Stage 2
gives the exact lam_min = (R+1) - |u||w|, u_i = rho^i, w_i = rho^-i.  Expanding,
    lam_min = -2 eps^2 [ S2 - S1^2/(R+1) ] + O(eps^3),   S1 = sum i, S2 = sum i^2,
and  S2 - S1^2/(R+1) = sum_{i=0}^{R} (i - R/2)^2 = R(R+1)(R+2)/12,  so
    lam_min = -eps^2 R(R+1)(R+2)/6 = -C(R+2,3) eps^2,
reproducing the Stage 2 law EXACTLY -- and exhibiting it as 2 eps^2 times the CENTRED
second moment of the window.  The centring is not a convention: it falls out of the
optimisation (the minimiser removes the mean).  The naive uncentred Cauchy-Schwarz bound
sum k^2 = R(R+1)(2R+1)/6 overshoots by (2R+1)/(R+2) -> 2.

SAME BOUND UP TO CONVENTION.  Identifying the window LENGTH (R+1 lattice points <-> an
interval of length 2L, i.e. R <-> 2L) and the pair-vs-quartet factor 2:
    function field:  C_R   = R^3/6  + O(R^2)
    zeta         :  C(L)  = 4L^3/3 = (2L)^3/6 .
They are the same law, (window length)^3 / 6, in the two settings.
"""
import sys, os, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from zeta_window import build_matrix, F, norm2, EVEN_D
from weilform import eigsym

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def Fv_and_derivs(v, L, g, sector=EVEN_D, h=None):
    """F_v(gamma), F_v'(gamma), F_v''(gamma) for f = sum v_k f_k/||f_k||."""
    L = mp.mpf(L)
    nrm = [mp.sqrt(norm2(k, L, sector)) for k in range(len(v))]
    fun = lambda t: sum(v[k] * F(k, L, t, sector) / nrm[k] for k in range(len(v)))
    return fun(g), mp.diff(fun, g), mp.diff(fun, g, 2)


def quartet_value(v, L, g, d, sector=EVEN_D):
    """4 Re F(w)^2 with w = g + i d, evaluated exactly (no expansion)."""
    L = mp.mpf(L)
    nrm = [mp.sqrt(norm2(k, L, sector)) for k in range(len(v))]
    w = mp.mpc(mp.mpf(g), mp.mpf(d))
    Sw = sum(v[k] * F(k, L, w, sector) / nrm[k] for k in range(len(v)))
    Sm = sum(v[k] * F(k, L, -w, sector) / nrm[k] for k in range(len(v)))
    return 4 * mp.re(Sw * Sm)


def cs_bound(L, g):
    """\\int_{-L}^{L} u^2 sin^2(g u) du  -- the refined Cauchy-Schwarz bound on |F'|^2."""
    L, g = mp.mpf(L), mp.mpf(g)
    return mp.quad(lambda u: u ** 2 * mp.sin(g * u) ** 2, [-L, 0, L])


def main():
    mp.mp.dps = 60
    g1 = mp.im(mp.zetazero(1))
    rows = []
    print("TASK 2a -- exact quartet value vs the delta^2 expansion, at the window "
          "minimiser\n")
    print(f"{'L':>6s} {'N':>4s} {'F(g1)':>13s} {'F prime':>13s} "
          f"{'exact 4ReF(w)^2':>18s} {'expansion':>18s} {'rel err':>11s}")
    for Lx, N in [('0.8', 16), ('1.0', 16), ('1.3', 20), ('1.6', 20), ('2.0', 24)]:
        L = mp.mpf(Lx)
        M = build_matrix(L, N, 60, EVEN_D)
        vals, vecs = eigsym(M)
        v = vecs[0]
        nv = mp.sqrt(sum(x ** 2 for x in v))
        v = [x / nv for x in v]
        F0, F1, F2 = Fv_and_derivs(v, L, g1)
        d = mp.mpf('1e-3')
        exact = quartet_value(v, L, g1, d)
        expan = 4 * (F0 ** 2 - d ** 2 * (F1 ** 2 + F0 * F2))
        rel = abs(exact - expan) / max(abs(exact), mp.mpf(10) ** -60)
        print(f"{Lx:>6s} {N:>4d} {mp.nstr(F0,6):>13s} {mp.nstr(F1,6):>13s} "
              f"{mp.nstr(exact,10):>18s} {mp.nstr(expan,10):>18s} {mp.nstr(rel,4):>11s}")
        # Task 2b / Task 3b
        C_pred = 4 * F1 ** 2
        B = cs_bound(L, g1)
        rows.append(dict(L=Lx, N=N, lam_star=mp.nstr(vals[0], 10),
                         F_at_g1=mp.nstr(F0, 10), Fp_at_g1=mp.nstr(F1, 10),
                         Fpp_at_g1=mp.nstr(F2, 10),
                         C_pred_4Fp2=mp.nstr(C_pred, 10),
                         cs_bound_Fp2=mp.nstr(B, 10),
                         C_cs_bound_4B=mp.nstr(4 * B, 10),
                         ratio_C_over_CS=mp.nstr(C_pred / (4 * B), 8),
                         C_over_L3=mp.nstr(C_pred / L ** 3, 8),
                         four_thirds=mp.nstr(mp.mpf(4) / 3, 8),
                         laguerre_L1=mp.nstr(F1 ** 2 - F0 * F2, 10),
                         ours_Fp2_plus_FFpp=mp.nstr(F1 ** 2 + F0 * F2, 10)))
    print("\nTASK 2b / 3b -- predicted C = 4 F'(g1)^2 vs the Cauchy-Schwarz bound 4 L^3/3\n")
    print(f"{'L':>6s} {'C=4F prime^2':>15s} {'CS bound 4B':>14s} {'C/CS':>9s} "
          f"{'C/L^3':>10s} {'4/3':>8s} {'L1=Fp2-FFpp':>14s} {'ours Fp2+FFpp':>15s}")
    for r in rows:
        print(f"{r['L']:>6s} {r['C_pred_4Fp2']:>15s} {r['C_cs_bound_4B']:>14s} "
              f"{r['ratio_C_over_CS']:>9s} {r['C_over_L3']:>10s} {r['four_thirds']:>8s} "
              f"{r['laguerre_L1']:>14s} {r['ours_Fp2_plus_FFpp']:>15s}")
    with open(os.path.join(ROOT, 'results', 'stage3_delta2_law.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("\nwrote results/stage3_delta2_law.csv")

    print("\nTASK 3c -- function-field law: exact vs centred/uncentred second moment")
    print(f"{'R':>4s} {'exact C(R+2,3)':>16s} {'2*centred moment':>18s} "
          f"{'uncentred sum k^2':>19s} {'uncentred/exact':>16s}")
    for R in (1, 2, 4, 8, 16, 32, 64):
        exact = mp.mpf(R * (R + 1) * (R + 2)) / 6
        cent = 2 * mp.mpf(R * (R + 1) * (R + 2)) / 12
        unc = mp.mpf(R * (R + 1) * (2 * R + 1)) / 6
        print(f"{R:>4d} {mp.nstr(exact,10):>16s} {mp.nstr(cent,10):>18s} "
              f"{mp.nstr(unc,10):>19s} {mp.nstr(unc/exact,6):>16s}")
    print("\n  R <-> 2L identification:  C_R = R^3/6  vs  C(L) = 4L^3/3 = (2L)^3/6  "
          "-- the same law.")


if __name__ == '__main__':
    main()


def perturbed_C(Lx, N, dps, gs, delta, sector=EVEN_D):
    r"""C(L) measured the way Stage 3e measures it: from lam_min of the PERTURBED form,
    together with 4 F'^2 evaluated at the PERTURBED minimiser (not at the unperturbed
    ground state).  These are the two things Task 2b asks to compare."""
    import sys as _s, os as _o
    _s.path.insert(0, _o.path.dirname(_o.path.abspath(__file__)))
    from stage3e_certify import move_block
    mp.mp.dps = dps
    L = mp.mpf(Lx)
    M = build_matrix(L, N, dps, sector)
    Z = M + move_block(L, range(N), gs, delta, sector)
    vals, vecs = eigsym(Z)
    v = vecs[0]
    nv = mp.sqrt(sum(x ** 2 for x in v))
    v = [x / nv for x in v]
    F0, F1, F2 = Fv_and_derivs(v, L, gs, sector)
    g0 = eigsy_min_vec(M)
    F0g, F1g, F2g = Fv_and_derivs(g0, L, gs, sector)
    return dict(lam_min=vals[0], C_measured=-vals[0] / mp.mpf(delta) ** 2,
                C_at_perturbed_min=4 * (F1 ** 2 + F0 * F2),
                C_at_unperturbed_gs=4 * (F1g ** 2 + F0g * F2g),
                F_at_g=F0, Fp_at_g=F1)


def eigsy_min_vec(M):
    vals, vecs = eigsym(M)
    v = vecs[0]
    nv = mp.sqrt(sum(x ** 2 for x in v))
    return [x / nv for x in v]
