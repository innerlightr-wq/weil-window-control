#!/usr/bin/env python3
"""Stage 1: function-field control for the CvS/Connes window Weil form.

Everything upstream of the eigensolve is EXACT integer arithmetic (T1):
point counts -> L-polynomial -> integer power sums p_n.
Only t(n) = p_n q^{-n/2} and the eigensolve are numerical (T2), at --dps digits.
"""
import sys, os, csv, argparse, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from curves import (squarefree_over_Fq, count_points_hyperelliptic,
                    newton_p_to_A, A_to_power_sums)
from weilform import toeplitz, eigsym, parity_of, poly_zeros, inertia
from analysis import analyse_window

CURVES = [
    # tag, description, f (little-endian), q, g
    ('E5',  'y^2 = x^3 + x  / F_5',            [0, 1, 0, 1],             5, 1),
    ('H2F3','y^2 = x^5 + x^2 + 1  / F_3',      [1, 0, 1, 0, 0, 1],       3, 2),
    ('G3F5','y^2 = x^7 + x + 1  / F_5',        [1, 1, 0, 0, 0, 0, 0, 1], 5, 3),
    ('G3F3','y^2 = x^7 + x^3 + x + 1  / F_3',  [1, 1, 0, 1, 0, 0, 0, 1], 3, 3),
]


def build_curve(tag, desc, f, q, g, brute_budget=200000):
    assert squarefree_over_Fq(f, q), f"{tag}: f not squarefree over F_{q}"
    assert (len(f) - 1) == 2 * g + 1, f"{tag}: deg f must be 2g+1"
    nb = 1
    while q ** (nb + 1) <= brute_budget:
        nb += 1
    N_brute = {n: count_points_hyperelliptic(f, q, n) for n in range(1, nb + 1)}
    p_brute = {n: 1 + q ** n - N_brute[n] for n in N_brute}

    # H3: build P from N_1..N_g only, then predict the rest.
    A = newton_p_to_A([None] + [p_brute[n] for n in range(1, g + 1)], q, g)
    p_all = A_to_power_sums(A, max(4 * g + 8, nb))
    N_pred = {n: 1 + q ** n - p_all[n] for n in range(1, nb + 1)}
    h3_ok = all(N_pred[n] == N_brute[n] for n in N_brute)
    h3_detail = {n: (N_brute[n], N_pred[n]) for n in sorted(N_brute)}

    # cross-check: L-polynomial from ALL brute-forced p_n (no functional equation)
    if nb >= 2 * g:
        from fractions import Fraction
        e = [Fraction(1)] + [Fraction(0)] * (2 * g)
        for k in range(1, 2 * g + 1):
            s = Fraction(0)
            for i in range(1, k + 1):
                s += Fraction(-1) ** (i - 1) * e[k - i] * p_brute[i]
            e[k] = s / k
        A_full = [int(Fraction(-1) ** i * e[i]) for i in range(2 * g + 1)]
        A_consistent = (A_full == A)
    else:
        A_full, A_consistent = None, None

    return dict(tag=tag, desc=desc, f=f, q=q, g=g, nb=nb, N_brute=N_brute,
                A=A, A_full=A_full, A_consistent=A_consistent,
                p=p_all, h3_ok=h3_ok, h3_detail=h3_detail)


def frobenius_angles(A, q, g):
    """theta_j from the L-polynomial, high precision. Also returns |alpha_j|/sqrt(q)."""
    # P(T) = sum A_i T^i has roots 1/alpha_j.  The alpha_j are the roots of the
    # reversed polynomial  Q(T) = T^{2g} P(1/T) = sum_i A_i T^{2g-i},
    # whose highest-degree-first coefficient list is exactly [A_0, A_1, ..., A_2g].
    rev = [mp.mpf(A[i]) for i in range(2 * g + 1)]
    roots = mp.polyroots(rev, maxsteps=300, extraprec=300)
    sq = mp.sqrt(q)
    out = []
    for r in roots:
        out.append((mp.arg(r), abs(r) / sq, r))
    out.sort(key=lambda z: z[0])
    return out


def t_values(p, q, nmax):
    sq = mp.sqrt(mp.mpf(q))
    t = [mp.mpf(0)] * (nmax + 1)
    for n in range(0, nmax + 1):
        t[n] = mp.mpf(p[n]) / sq ** n if n > 0 else mp.mpf(p[0])
    return t


def analyse(cur, dps, rows, zrows):
    q, g, tag = cur['q'], cur['g'], cur['tag']
    Rmax = 4 * g
    p = dict(cur['p']) if isinstance(cur['p'], dict) else {n: v for n, v in enumerate(cur['p'])}
    p[0] = 2 * g
    t = t_values(p, q, Rmax)
    ang = frobenius_angles(cur['A'], q, g)
    thetas = [a[0] for a in ang]
    maxdev = max(abs(a[1] - 1) for a in ang)     # |alpha|/sqrt(q) - 1 : Weil check

    n_distinct = len(set(round(float(th), 12) for th in thetas))
    print(f"\n=== {tag}: {cur['desc']}  (g={g}, q={q}) ===")
    print(f"  L-poly P(T) coeffs A_0..A_2g = {cur['A']}")
    print(f"  H3 (N_1..N_g determine all N_n up to n={cur['nb']}): "
          f"{'CONFIRMED' if cur['h3_ok'] else 'REFUTED'}   "
          f"(independent integer predictions: {cur['nb'] - g})")
    if cur['A_consistent'] is not None:
        print(f"  L-poly from all p_1..p_2g equals functional-equation build: {cur['A_consistent']}")
    print(f"  max | |alpha_j|/sqrt(q) - 1 | = {mp.nstr(maxdev, 5)}   (Weil/RH-for-curves check)")
    print(f"  theta_j = {[mp.nstr(th, 10) for th in thetas]}  (#distinct={n_distinct})")
    print(f"  t(0..{Rmax}) = {[mp.nstr(x, 8) for x in t]}")

    # ---- exact prediction for the R=2g kernel vector (T1 derivation) ----------
    # kernel poly = prod_j (z - alpha_j/sqrt q) = sum_k A_{2g-k} q^{-(2g-k)/2} z^k,
    # and the functional equation A_{2g-k}=q^{g-k}A_k collapses this to c_k ~ A_k q^{-k/2},
    # which is palindromic -- so the R=2g ground state is forced EVEN.
    sq = mp.sqrt(mp.mpf(q))
    c_pred = [mp.mpf(cur['A'][k]) / sq ** k for k in range(2 * g + 1)]
    npred = mp.sqrt(sum(x ** 2 for x in c_pred))
    c_pred = [x / npred for x in c_pred]

    lam_seq = []
    def vals_lam(row):
        return mp.mpf(row['lam_min'])
    for R in range(0, Rmax + 1):
        row, v, zs = analyse_window(t, R, thetas, dps)
        row.update(curve=tag, q=q, g=g, rank_pred_theory=min(R + 1, n_distinct))
        row['H2_kernel_vec_err'] = ''
        if R == 2 * g:
            nv = mp.sqrt(sum(x ** 2 for x in v))
            w = [x / nv for x in v]
            e1 = max(abs(a - b) for a, b in zip(w, c_pred))
            e2 = max(abs(a + b) for a, b in zip(w, c_pred))
            row['H2_kernel_vec_err'] = mp.nstr(min(e1, e2), 6)
        rows.append(row)
        zang = sorted(mp.arg(z) for z in zs)
        for za, zz in zip(zang, sorted(zs, key=lambda z: mp.arg(z))):
            zrows.append(dict(curve=tag, R=R, zero_angle=mp.nstr(za, 14),
                              abs_zero=mp.nstr(abs(zz), 14)))
        lam_seq.append(vals_lam(row))
        print(f"  R={R:2d} n={R+1:2d} lam_min={row['lam_min']:>15s}[{row['lam_min_sign']:>4s}]"
              f" f64={row['lam_min_float64']:>12s}[{row['float64_sign']}]"
              f" gap={row['gap_lam2_lam1']:>10s} simple={str(row['lam_min_simple']):5s}"
              f" par={row['parity']:6s} H1:max||z|-1|={row['H1_max_abs_zero_minus_1']:>10s}"
              f" H4 z->th={row['H4_zero_to_theta']:>9s} th->z={row['H4_theta_to_zero']:>9s}"
              f" inert={row['inertia']} rk_th={min(R+1,n_distinct)}"
              + (f" H2err={row['H2_kernel_vec_err']}" if R == 2 * g else ""))
    # ---- variational (Cauchy-interlacing) monotonicity: lam_min(R+1) <= lam_min(R) ----
    tolm = mp.mpf(10) ** (-(dps - 10))
    viol = [(R, mp.nstr(lam_seq[R + 1] - lam_seq[R], 4))
            for R in range(len(lam_seq) - 1) if lam_seq[R + 1] - lam_seq[R] > tolm]
    print(f"  variational monotonicity lam_min(R+1) <= lam_min(R): "
          f"{'HOLDS' if not viol else 'VIOLATED ' + str(viol)}"
          f"   (tol {mp.nstr(tolm, 3)})")
    print(f"  lam_min(R) sequence: {[mp.nstr(x, 6) for x in lam_seq]}")
    # pole-part rank-2 identity: A(n)=q^{n/2}+q^{-n/2} => A = u v^T + v u^T  (T1)
    sqq = mp.sqrt(mp.mpf(q))
    u = [sqq ** i for i in range(Rmax + 1)]
    w = [sqq ** (-i) for i in range(Rmax + 1)]
    resid = max(abs(u[i] * w[j] + w[i] * u[j] - (sqq ** abs(i - j) + sqq ** (-abs(i - j))))
                for i in range(Rmax + 1) for j in range(Rmax + 1))
    print(f"  pole-part rank-2 identity residual (A = u w^T + w u^T): {mp.nstr(resid, 4)}")
    return thetas


def recorder_split(cur, dps, srows):
    """T = A - 2M with A = 2g*I + pole part, 2M = prime part (brief's convention).
    The n=0 assignment is a CONVENTION (documented); both choices reported."""
    q, g, tag = cur['q'], cur['g'], cur['tag']
    Rmax = 4 * g
    sq = mp.sqrt(mp.mpf(q))
    p = {n: v for n, v in enumerate(cur['p'])}
    p[0] = 2 * g
    t = t_values(p, q, Rmax)
    pole = [sq ** n + sq ** (-n) for n in range(Rmax + 1)]              # n=0 -> 2
    Nn = {n: 1 + q ** n - p[n] for n in range(1, Rmax + 1)}
    for conv in ('brief', 'uniform'):
        if conv == 'brief':      # A(0)=2g+2, 2M(0)=2
            Avec = [pole[0] + 2 * g] + pole[1:]
            Mvec = [mp.mpf(2)] + [mp.mpf(Nn[n]) / sq ** n for n in range(1, Rmax + 1)]
        else:                    # A(0)=2, 2M(0)=N_0:=2-2g
            Avec = list(pole)
            Mvec = [mp.mpf(2 - 2 * g)] + [mp.mpf(Nn[n]) / sq ** n for n in range(1, Rmax + 1)]
        for R in range(0, Rmax + 1):
            TA = toeplitz(Avec, R)
            TM = toeplitz([x / 2 for x in Mvec], R)
            TS = toeplitz(t, R)
            resid = max(abs(TA[i, j] - 2 * TM[i, j] - TS[i, j])
                        for i in range(R + 1) for j in range(R + 1))
            iA, _, _ = inertia(TA); iM, _, _ = inertia(TM); iS, _, _ = inertia(TS)
            srows.append(dict(curve=tag, convention=conv, R=R, size=R + 1,
                              split_residual=mp.nstr(resid, 4),
                              A_inertia=f"({iA[0]},{iA[1]},{iA[2]})",
                              M_inertia=f"({iM[0]},{iM[1]},{iM[2]})",
                              S_inertia=f"({iS[0]},{iS[1]},{iS[2]})",
                              A_indefinite=(iA[0] > 0 and iA[2] > 0),
                              M_indefinite=(iM[0] > 0 and iM[2] > 0),
                              S_psd=(iS[0] == 0)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dps', type=int, default=60)
    ap.add_argument('--out', default=os.path.join(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))), 'data'))
    a = ap.parse_args()
    mp.mp.dps = a.dps
    print(f"mpmath working precision: {mp.mp.dps} decimal digits "
          f"(backend: {'gmpy' if mp.libmp.BACKEND == 'gmpy' else mp.libmp.BACKEND})")
    rows, zrows, srows = [], [], []
    built = []
    for spec in CURVES:
        cur = build_curve(*spec)
        built.append(cur)
        analyse(cur, a.dps, rows, zrows)
        recorder_split(cur, a.dps, srows)
    os.makedirs(a.out, exist_ok=True)
    for name, data in (('stage1_spectra.csv', rows), ('stage1_zeros.csv', zrows),
                       ('stage1_recorder_split.csv', srows)):
        with open(os.path.join(a.out, name), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0].keys()))
            w.writeheader(); w.writerows(data)
        print(f"wrote {os.path.join(a.out, name)}  ({len(data)} rows)")
    with open(os.path.join(a.out, 'stage1_curves.json'), 'w') as fh:
        json.dump([{k: (v if k not in ('N_brute', 'h3_detail', 'p') else
                        {str(kk): vv for kk, vv in (v.items() if isinstance(v, dict) else enumerate(v))})
                    for k, v in c.items()} for c in built], fh, indent=1)


if __name__ == '__main__':
    main()
