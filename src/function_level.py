#!/usr/bin/env python3
"""Items 1 and 2: the kernel is an ideal, and FUNCTION-level re-scoring of step (iv).

THE STRUCTURAL FACT (T1).  Let P*(z) = prod over DISTINCT beta_j of (z - beta_j),
degree d.  For c = (c_0..c_R), Q(c) = 0 iff chat(beta_j) = 0 for all j iff P* | chat.
Hence

    ker T_R = { P*(z) q(z) : deg q <= R - d },     dim ker T_R = max(0, R + 1 - d).

Consequences.
 (2a) At R = d the kernel is 1-dimensional: lam_min = 0 is SIMPLE, and the ground state
      IS P* itself.  At R = d+1 the kernel is 2-dimensional, so lam_min = 0 is NOT
      simple and the CvS hypothesis (simple isolated lowest eigenvalue) FAILS for every
      R > d.  A finite spectrum simply runs out of content past its own degree.
 (2b) P* is self-reciprocal up to sign (the spectrum is closed under beta -> 1/beta by
      the functional equation), so J(P* q) = +- P* q^rev: the kernel is J-invariant and
      splits into even and odd parts of dimensions ceil((m+1)/2), floor((m+1)/2) with
      m = R - d.  For m >= 1 BOTH parts are nonzero, so a generic kernel vector -- which
      is what any eigensolver returns -- has no definite parity.  The "evenness failure"
      observed at R > 2g in Stage 1 is therefore an artefact of degeneracy, not a
      property of the curve.

RE-SCORING (iv) AT THE FUNCTION LEVEL (T1 + T2).  In CvS, (iv) is Hurwitz transfer:
the NORMALIZED MINIMIZER FUNCTION converges to the target, and zero convergence is a
corollary.  The converse fails, and that is the whole point:

   THEOREM (T1).  Suppose lam_min(T_R) is simple for infinitely many R, and the
   normalized ground-state polynomials converge locally uniformly on C to some F != 0.
   By H1 every zero of every chat_R lies on |z| = 1; by Hurwitz every zero of F is a
   limit of zeros of chat_R, hence |z| = 1.  So if the target has ANY zero off the unit
   circle, the minimizer functions CANNOT converge to it.

So H1 -- the free witness -- is exactly the obstruction to function-level convergence in
the off-line case.  (i) and (iv) are two sides of one coin, and (iv) is NOT free.
Measured below with three distances.
"""
import sys, os, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from weilform import toeplitz, eigsym, parity_of, poly_zeros, inertia
from analysis import analyse_window


def target_poly(betas, tol=None):
    """Monic P*(z) = prod over DISTINCT beta of (z - beta), little-endian coefficients."""
    tol = tol or mp.mpf(10) ** (-(mp.mp.dps // 2))
    uniq = []
    for b in betas:
        if not any(abs(b - u) < tol for u in uniq):
            uniq.append(b)
    p = [mp.mpc(1)]
    for b in uniq:
        p = [-b * p[0]] + [p[k - 1] - b * p[k] for k in range(1, len(p))] + [p[-1]]
    return p, uniq


def _mgs(basis):
    """Modified Gram-Schmidt on a list of mp vectors (lists); returns orthonormal set."""
    Q = []
    for v in basis:
        w = list(v)
        for u in Q:
            c = sum(mp.conj(u[k]) * w[k] for k in range(len(w)))
            w = [w[k] - c * u[k] for k in range(len(w))]
        n = mp.sqrt(sum(abs(x) ** 2 for x in w))
        if n > mp.mpf(10) ** (-(mp.mp.dps - 10)):
            Q.append([x / n for x in w])
    return Q


def dist_to_ideal(v, Pstar, R):
    """Relative distance from v (length R+1) to span{z^k P*: k=0..R-d}, in [0,1]."""
    d = len(Pstar) - 1
    if R < d:
        return mp.mpf(1)
    basis = []
    for k in range(R - d + 1):
        e = [mp.mpc(0)] * (R + 1)
        for i, c in enumerate(Pstar):
            e[i + k] = c
        basis.append(e)
    Q = _mgs(basis)
    nv = mp.sqrt(sum(abs(x) ** 2 for x in v))
    w = [mp.mpc(x) / nv for x in v]
    proj = [mp.mpc(0)] * (R + 1)
    for u in Q:
        c = sum(mp.conj(u[k]) * w[k] for k in range(R + 1))
        proj = [proj[k] + c * u[k] for k in range(R + 1)]
    return mp.sqrt(sum(abs(w[k] - proj[k]) ** 2 for k in range(R + 1)))


def coeff_distance(v, Pstar):
    """l2 distance between unit-normalized v and unit-normalized P*, best sign.
    Only meaningful when len(v) == len(Pstar)."""
    a = list(v); b = list(Pstar)
    na = mp.sqrt(sum(abs(x) ** 2 for x in a)); nb = mp.sqrt(sum(abs(x) ** 2 for x in b))
    a = [x / na for x in a]; b = [x / nb for x in b]
    d1 = mp.sqrt(sum(abs(a[k] - b[k]) ** 2 for k in range(len(a))))
    d2 = mp.sqrt(sum(abs(a[k] + b[k]) ** 2 for k in range(len(a))))
    return min(d1, d2)


def normalized_value_at(v, beta):
    """|chat(beta)| / (||v||_2 * (sum |beta|^{2k})^{1/2})  in [0,1] by Cauchy-Schwarz.
    Zero iff beta is a zero of chat. This is the FUNCTION-level 'does the minimizer
    see the off-circle point' statistic."""
    R = len(v) - 1
    val = sum(mp.mpc(v[k]) * beta ** k for k in range(R + 1))
    nv = mp.sqrt(sum(abs(x) ** 2 for x in v))
    r = abs(beta)
    bound = mp.sqrt(sum(r ** (2 * k) for k in range(R + 1)))
    return abs(val) / (nv * bound)


def kernel_parity_dims(T, R, tol_exp=None):
    """Dimensions of the even and odd parts of ker T_R, from the numerical kernel."""
    vals, vecs = eigsym(T)
    scale = max([abs(v) for v in vals] + [mp.mpf(1)])
    tol = scale * mp.mpf(10) ** (-(tol_exp or (mp.mp.dps - 12)))
    K = [vecs[m] for m in range(R + 1) if abs(vals[m]) <= tol]
    if not K:
        return 0, 0, 0
    ev, od = [], []
    for v in K:                                   # split each basis vector by J
        rev = list(reversed(v))
        ev.append([(a + b) / 2 for a, b in zip(v, rev)])
        od.append([(a - b) / 2 for a, b in zip(v, rev)])
    de = len(_mgs([[mp.mpc(x) for x in w] for w in ev]))
    do = len(_mgs([[mp.mpc(x) for x in w] for w in od]))
    return len(K), de, do


def run(tag, desc, t, betas, Rmax, dps, rows):
    Pstar, uniq = target_poly(betas)
    d = len(Pstar) - 1
    offline = [b for b in uniq if abs(abs(b) - 1) > mp.mpf(10) ** -20]
    print(f"\n=== {tag} ===\n  {desc}")
    print(f"  target P*: degree d = {d} (#distinct beta); off-circle beta: "
          f"{[mp.nstr(b, 6) for b in offline] or 'none'}")
    for R in range(1, Rmax + 1):
        T = toeplitz([mp.mpf(x) for x in t], R)
        row, v, zs = analyse_window(t, R, [mp.arg(b) for b in uniq], dps)
        nk, de, do = kernel_parity_dims(T, R)
        di = dist_to_ideal(v, Pstar, R) if R >= d else None
        cd = coeff_distance(v, Pstar) if R == d else None
        gaps = [normalized_value_at(v, b) for b in offline]
        radial = [abs(abs(b) - 1) for b in offline]
        rec = dict(config=tag, R=R, d=d, lam_min=row['lam_min'],
                   lam_min_sign=row['lam_min_sign'],
                   lam_min_simple=row['lam_min_simple'], parity=row['parity'],
                   kernel_dim=nk, kernel_dim_pred=max(0, R + 1 - d),
                   kernel_even_dim=de, kernel_odd_dim=do,
                   M1_coeff_dist_at_R_eq_d=(mp.nstr(cd, 6) if cd is not None else ''),
                   M1p_dist_to_ideal=(mp.nstr(di, 6) if di is not None else ''),
                   M2_min_normalized_value_at_offline_beta=(
                       mp.nstr(min(gaps), 6) if gaps else ''),
                   radial_gap=(mp.nstr(min(radial), 6) if radial else ''),
                   H1_max_abs_zero_minus_1=row['H1_max_abs_zero_minus_1'])
        rows.append(rec)
        print(f"   R={R:2d} lam={row['lam_min']:>15s}[{row['lam_min_sign']:>4s}]"
              f" simple={str(row['lam_min_simple']):5s} par={row['parity']:6s}"
              f" ker={nk}(pred {max(0,R+1-d)}) even/odd={de}/{do}"
              f"  M1'(dist to ideal)={rec['M1p_dist_to_ideal']:>10s}"
              f"  M2(|chat(beta)| norm)={rec['M2_min_normalized_value_at_offline_beta']:>10s}"
              + (f"  M1(coeff dist at R=d)={rec['M1_coeff_dist_at_R_eq_d']}" if cd is not None else ""))


def main():
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument('--dps', type=int, default=50)
    ap.add_argument('--Rmax', type=int, default=14)
    a = ap.parse_args(); mp.mp.dps = a.dps
    from curves import count_points_hyperelliptic, newton_p_to_A, A_to_power_sums
    from stage1 import CURVES
    from planted import t_from_blocks, betas as planted_betas
    rows = []
    print(f"mpmath dps = {mp.mp.dps}\n")
    print("############ HONEST CURVES (item 2: kernel dimension for R > 2g) ############")
    for tag, desc, f, q, g in CURVES:
        p = [None] + [1 + q ** n - count_points_hyperelliptic(f, q, n)
                      for n in range(1, g + 1)]
        A = newton_p_to_A(p, q, g)
        pall = A_to_power_sums(A, a.Rmax)
        sq = mp.sqrt(mp.mpf(q))
        t = [mp.mpf(2 * g)] + [mp.mpf(pall[n]) / sq ** n for n in range(1, a.Rmax + 1)]
        al = mp.polyroots([mp.mpf(A[i]) for i in range(2 * g + 1)],
                          maxsteps=300, extraprec=300)
        run(tag, desc, t, [x / sq for x in al], min(a.Rmax, 2 * g + 4), a.dps, rows)
    # the repeated-angle counterexample
    f, q, g = [0, 1, 0, 0, 0, 1], 5, 2
    p = [None] + [1 + q ** n - count_points_hyperelliptic(f, q, n) for n in (1, 2)]
    A = newton_p_to_A(p, q, g); pall = A_to_power_sums(A, a.Rmax); sq = mp.sqrt(mp.mpf(q))
    t = [mp.mpf(2 * g)] + [mp.mpf(pall[n]) / sq ** n for n in range(1, a.Rmax + 1)]
    al = mp.polyroots([mp.mpf(A[i]) for i in range(2 * g + 1)], maxsteps=300, extraprec=300)
    run('E5x-repeated', 'y^2=x^5+x / F_5, P=(1+5T^2)^2, repeated angles',
        t, [x / sq for x in al], 8, a.dps, rows)

    print("\n############ PLANTED (item 1: function-level re-scoring of (iv)) ############")
    for tag, desc, blocks in [
            ('C0-online', 'CONTROL: on-line, +-0.7, +-2.0', [('online', 0.7), ('online', 2.0)]),
            ('C1-g1real', 'off-line real pair rho=1.3', [('realpair', '1.3')]),
            ('C2-quartet', 'off-line quartet rho=1.3, phi=0.7', [('quartet', '1.3', 0.7)]),
            ('C4-mixed', 'off-line quartet + on-line pair', [('quartet', '1.3', 0.7), ('online', 2.0)])]:
        t, _ = t_from_blocks(blocks, a.Rmax)
        run(tag, desc, t, planted_betas(blocks), a.Rmax, a.dps, rows)

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'results', 'stage2_function_level.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {out}")


if __name__ == '__main__':
    main()
