#!/usr/bin/env python3
r"""Function-field cross-check of the parity sign structure.  EXACT over Q.  [T1]

Setting (Theorem 6 / notes/proofs.md A.1-A.2 of the parent repository): an off-circle real
pair of normalised Frobenius eigenvalues {rho, 1/rho} gives the window Toeplitz matrix

    T_R = u w^T + w u^T,     u_i = rho^i,  w_i = rho^{-i},   0 <= i <= R,

with the two nonzero eigenvalues (R+1) +- ||u|| ||w||.  This script verifies, exactly in
Fraction arithmetic for rational rho, the SIGNED RANK-ONE SPLIT and its parity structure,
which is the function-field mirror of the Phase-5 result (+2A^2 even / -2S^2 odd).

Let J be the flip (J x)_i = x_{R-i}.  Since u_{R-i} = rho^R w_i and w_{R-i} = rho^{-R} u_i,

    J u = rho^R w,      J w = rho^{-R} u,

so J does NOT fix u + w; it fixes the rho^R-WEIGHTED combination.  Put

    v_+ = u + rho^R w        (J-palindromic:      J v_+ = + v_+)
    v_- = u - rho^R w        (J-anti-palindromic: J v_- = - v_-)

Then, exactly over Q,

    T_R          = (1 / (2 rho^R)) [ v_+ v_+^T  -  v_- v_-^T ]
    T_R v_+      = lam_+ v_+,    lam_+ = (R+1) + rho^{-R} ||u||^2   > 0
    T_R v_-      = lam_- v_-,    lam_- = (R+1) - rho^{-R} ||u||^2   <= 0
    rho^{-R} ||u||^2 = ||u|| ||w||   (because ||w||^2 = rho^{-2R} ||u||^2 identically)

so the J-EVEN direction carries the POSITIVE rank-one piece and the J-ODD direction the
NEGATIVE one -- the same signed split as the zeta/Dirichlet side.

TWO CORRECTIONS to the commissioning brief, both checked below rather than argued:
  * the brief's identity  T_R = (1/2)[(u+w)(u+w)^T - (u-w)(u-w)^T]  is TRUE, exactly;
  * but u+w and u-w are NOT the eigenvectors and NOT J-parity vectors.  They are both
    iff rho^R = 1.  The rho^R weight is exactly what J requires.
"""
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAIL = []
OUT = {}


def check(name, ok, detail=''):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ''))
    if not ok:
        FAIL.append(name)
    return ok


def outer(a, b):
    return [[a[i] * b[j] for j in range(len(b))] for i in range(len(a))]


def madd(A, B, ca=Fr(1), cb=Fr(1)):
    return [[ca * A[i][j] + cb * B[i][j] for j in range(len(A[0]))]
            for i in range(len(A))]


def mscale(A, c):
    return [[c * A[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def mzero(A):
    return all(A[i][j] == 0 for i in range(len(A)) for j in range(len(A[0])))


def matvec(A, x):
    return [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(A))]


def dot(a, b):
    return sum(a[i] * b[i] for i in range(len(a)))


def flip(x):
    return list(reversed(x))


def parallel(a, b):
    """exact: are a, b parallel?  (cross products all vanish)"""
    n = len(a)
    return all(a[i] * b[j] == a[j] * b[i] for i in range(n) for j in range(n))


def ratio(a, b):
    """exact lambda with a = lambda b, assuming parallel and b nonzero."""
    for i in range(len(b)):
        if b[i] != 0:
            return Fr(a[i], 1) / b[i]
    return None


def main():
    RHOS = [Fr(2), Fr(3, 2), Fr(5, 4), Fr(7, 3), Fr(1)]
    RS = [1, 2, 4, 8, 16]
    rows = []
    print("=" * 100)
    print("Exact (Fraction) verification of the signed parity split of T_R")
    print("=" * 100)
    hdr = (f"  {'rho':>7s} {'R':>3s} {'split v+-':>10s} {'eigvec v+':>10s} "
           f"{'eigvec v-':>10s} {'J-par v+':>9s} {'J-par v-':>9s} "
           f"{'lam_+ > 0':>10s} {'lam_- <= 0':>11s} {'brief split':>12s} "
           f"{'u+w eigvec':>11s} {'u+w J-even':>11s}")
    print(hdr)
    for rho in RHOS:
        for R in RS:
            n = R + 1
            u = [rho ** i for i in range(n)]
            w = [rho ** (-i) for i in range(n)]
            rhoR = rho ** R
            T = madd(outer(u, w), outer(w, u))

            # --- the rho^R-weighted parity vectors
            vp = [u[i] + rhoR * w[i] for i in range(n)]
            vm = [u[i] - rhoR * w[i] for i in range(n)]

            # (1) the signed rank-one split, exactly
            split = madd(outer(vp, vp), outer(vm, vm), Fr(1), Fr(-1))
            ok_split = mzero(madd(T, mscale(split, Fr(1, 2) / rhoR), Fr(1), Fr(-1)))

            # (2) eigenvector property + eigenvalues, exactly
            Tvp, Tvm = matvec(T, vp), matvec(T, vm)
            ok_evp = parallel(Tvp, vp)
            ok_evm = parallel(Tvm, vm) if any(x != 0 for x in vm) else True
            nu2 = dot(u, u)
            lam_p_pred = Fr(R + 1) + nu2 / rhoR
            lam_m_pred = Fr(R + 1) - nu2 / rhoR
            lam_p = ratio(Tvp, vp) if ok_evp else None
            lam_m = ratio(Tvm, vm) if (ok_evm and any(x != 0 for x in vm)) else None
            ok_lamp = (lam_p == lam_p_pred)
            ok_lamm = (lam_m == lam_m_pred) if lam_m is not None else (lam_m_pred == 0)
            # ||w||^2 = rho^{-2R} ||u||^2 identically  =>  rho^{-R}||u||^2 = ||u|| ||w||
            ok_norm = (dot(w, w) * rhoR ** 2 == nu2)

            # (3) parity, exactly
            ok_jp = (flip(vp) == vp)
            ok_jm = (flip(vm) == [-x for x in vm])
            ok_signp = (lam_p_pred > 0)
            ok_signm = (lam_m_pred <= 0)

            # (4) the brief's un-normalised claims
            upw = [u[i] + w[i] for i in range(n)]
            umw = [u[i] - w[i] for i in range(n)]
            brief_split = madd(outer(upw, upw), outer(umw, umw), Fr(1), Fr(-1))
            ok_brief = mzero(madd(T, mscale(brief_split, Fr(1, 2)), Fr(1), Fr(-1)))
            ok_upw_eig = parallel(matvec(T, upw), upw)
            ok_upw_par = (flip(upw) == upw)

            for nm, o in (('split', ok_split), ('eig v+', ok_evp), ('eig v-', ok_evm),
                          ('lam+', ok_lamp), ('lam-', ok_lamm), ('normid', ok_norm),
                          ('J v+', ok_jp), ('J v-', ok_jm),
                          ('sign+', ok_signp), ('sign-', ok_signm),
                          ('brief split', ok_brief)):
                if not o:
                    FAIL.append(f"rho={rho} R={R}: {nm}")
            print(f"  {str(rho):>7s} {R:>3d} {str(ok_split):>10s} {str(ok_evp):>10s} "
                  f"{str(ok_evm):>10s} {str(ok_jp):>9s} {str(ok_jm):>9s} "
                  f"{str(ok_signp):>10s} {str(ok_signm):>11s} {str(ok_brief):>12s} "
                  f"{str(ok_upw_eig):>11s} {str(ok_upw_par):>11s}")
            rows.append(dict(rho=str(rho), R=R, rho_R=str(rhoR),
                             lam_plus=str(lam_p_pred), lam_minus=str(lam_m_pred),
                             norm_u_sq=str(nu2), norm_w_sq=str(dot(w, w)),
                             uu_ww_product=str(nu2 / rhoR),
                             split_exact=ok_split, vplus_is_eigvec=ok_evp,
                             vminus_is_eigvec=ok_evm,
                             vplus_palindromic=ok_jp, vminus_antipalindromic=ok_jm,
                             lam_plus_positive=ok_signp, lam_minus_nonpositive=ok_signm,
                             brief_split_identity=ok_brief,
                             brief_upw_is_eigvec=ok_upw_eig,
                             brief_upw_palindromic=ok_upw_par))
    print()
    nontrivial = [r for r in rows if r['rho'] != '1']
    check("signed rank-one split T_R = (1/2rho^R)[v+v+^T - v-v-^T], exactly, all cells",
          all(r['split_exact'] for r in rows))
    check("v+ and v- are eigenvectors with lam_+ = (R+1) + rho^-R||u||^2 and "
          "lam_- = (R+1) - rho^-R||u||^2, exactly",
          all(r['vplus_is_eigvec'] and r['vminus_is_eigvec'] for r in rows))
    check("||w||^2 = rho^-2R ||u||^2 identically, so rho^-R||u||^2 = ||u||.||w|| in Q",
          all(r['split_exact'] for r in rows))
    check("v+ palindromic and v- anti-palindromic under i -> R-i, exactly, all cells",
          all(r['vplus_palindromic'] and r['vminus_antipalindromic'] for r in rows))
    check("lam_+ > 0 and lam_- <= 0 in every cell (Cauchy-Schwarz, equality iff rho = 1)",
          all(r['lam_plus_positive'] and r['lam_minus_nonpositive'] for r in rows))
    check("the brief's identity T_R = (1/2)[(u+w)(u+w)^T - (u-w)(u-w)^T] is also exact",
          all(r['brief_split_identity'] for r in rows))
    bad_eig = [r for r in nontrivial if r['brief_upw_is_eigvec']]
    bad_par = [r for r in nontrivial if r['brief_upw_palindromic']]
    check("CORRECTION: u+w is NOT an eigenvector for any rho != 1 "
          "(the rho^R weight is required)", len(bad_eig) == 0,
          f"{len(bad_eig)} counterexamples to the refutation" if bad_eig
          else f"refuted in all {len(nontrivial)} cells with rho != 1")
    check("CORRECTION: u+w is NOT palindromic for any rho != 1", len(bad_par) == 0,
          f"refuted in all {len(nontrivial)} cells with rho != 1")
    triv = [r for r in rows if r['rho'] == '1']
    check("at rho = 1 (the degenerate on-circle case) u+w DOES work, and lam_- = 0",
          all(r['brief_upw_is_eigvec'] and r['brief_upw_palindromic']
              and r['lam_minus'] == '0' for r in triv),
          f"{len(triv)} cells, consistent with Corollary A3 (strictness fails iff rho = 1)")
    with open(os.path.join(HERE, 'results', 'ff_parity_split.json'), 'w') as fh:
        json.dump(dict(rows=rows, failures=FAIL), fh, indent=1)
    print()
    if FAIL:
        print(f"FAILED: {len(FAIL)}")
        for f in FAIL[:10]:
            print(f"   - {f}")
        return 1
    print("FUNCTION-FIELD PARITY SPLIT: ALL EXACT OVER Q")
    return 0


if __name__ == '__main__':
    sys.exit(main())
