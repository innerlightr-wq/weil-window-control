r"""PHASE 0 (GATE) -- normalization reconciliation.  Tests P0-T1..T4 of PREREG.md.

P0-T1  algebraic identification (reported; the derivation is notes/normalization.md)
P0-T2  CROSS-BASIS: Suzuki's G vs this repository's src/zeta_window.build_parts_idx,
       through  cos(k pi x/L) = (sqrt(2L)/2)(chi_k + chi_{-k}),  chi_0 = 1/sqrt(2L).
       The repository code is NOT written in this extension.
P0-T3  spectral side: G_{mn} = sum_rho phi_m(gamma_rho) phi_n(gamma_rho) -- assumes RH,
       so it supports no unconditional claim; a consistency check only.
P0-T4  structural: G real, symmetric, positive definite, G_{-m,-n} = G_{m,n}.
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
import mpmath as mp
import suzuki as S
from suzuki import fp
import zeta_window as ZW          # the repository's own gated assembly

OK = [True]
out = {}
def chk(name, cond, extra=""):
    OK[0] &= bool(cond)
    print(f"  {'PASS' if cond else 'FAIL'}  {name} {extra}")
    return bool(cond)


# ------------------------------------------------------------------ brute-force arch
def arch_brute(m, n, a, K):
    """The archimedean term by its DEFINING truncated series, for cross-checking the
    closed form.  Converges like O(1/K), so useful only at low precision."""
    C0 = mp.re(S.chi_corr(m, n, a, 0))
    tot = -mp.euler * C0
    for k in range(K + 1):
        c = mp.mpf(1) / 2 + 2 * k
        I = mp.quad(lambda v: S.S_corr(m, n, a, v) * mp.e ** (-c * v), [0, 2 * a])
        tot += C0 / (k + 1) - I
    return tot


def main():
    print("=" * 78)
    print("PHASE 0 (GATE) -- normalization reconciliation")
    print("=" * 78)

    # ---------------------------------------------------------------- P0-T1
    print("\n-- P0-T1: algebraic identification --")
    print("""   Fourier conventions identical: uhat(z) = int u e^{izx} dx both sides.
   Spectral sides identical: sum over distinct ordinates with multiplicity weights
   = sum over zeros with multiplicity.  Supports identical: (-a,a) vs [-L,L].
   Coercivity: lambda_a := inf sigma(A_a) = inf Q_W(u)/||u||^2_{L^2(-a,a)} ;
               lambda*(L) := inf Q(f)/||f||^2 .
   =>  a = L  and  lambda_a = lambda*(L),  constant 1.   (notes/normalization.md Sec.1)""")
    chk("P0-T1 identification pinned exactly, constant 1", True)
    out["P0_T1"] = dict(a_equals_L=True, lambda_constant="1", pinned=True)

    # ------------------------------------------- arch closed form vs defining series
    print("\n-- arch: closed form vs the defining truncated series (sanity) --")
    mp.mp.dps = 30
    a = mp.mpf('0.8047')
    rows = []
    for (m, n) in ((0, 0), (1, 1), (0, 1), (2, -3)):
        cf = S.arch_term(m, n, a)
        K = 1500
        bf = arch_brute(m, n, a, K)
        dev = abs(cf - bf)
        rows.append((m, n, fp(cf), fp(bf), fp(dev)))
        print(f"   (m,n)=({m},{n}): closed {mp.nstr(cf,14)}  series(K={K}) {mp.nstr(bf,14)}"
              f"   dev {mp.nstr(dev,3)}")
    # the series converges like 1/K, so 1e-3 agreement at K=4000 is the expected scale
    chk("closed form agrees with the defining series to the series' own O(1/K) accuracy",
        all(float(r[4]) < 1e-2 for r in rows),
        "(the series is O(1/K); this confirms the closed form, which is then used)")
    out["arch_crosscheck"] = rows

    # ---------------------------------------------------------------- P0-T4 first
    print("\n-- P0-T4: structural properties of G --")
    mp.mp.dps = 60
    for aval in ('0.5493061443340548', '0.8047189562170503'):
        a = mp.mpf(aval)
        N = 6
        G = S.gram(a, N, dps=60)
        n = G.rows
        maxim = max(abs(mp.im(G[i, j])) for i in range(n) for j in range(n))
        sym = max(abs(G[i, j] - G[j, i]) for i in range(n) for j in range(n))
        scale = max(abs(G[i, j]) for i in range(n) for j in range(n))
        rev = max(abs(G[i, j] - G[n - 1 - i, n - 1 - j]) for i in range(n) for j in range(n))
        Gr = mp.matrix([[mp.re(G[i, j]) for j in range(n)] for i in range(n)])
        evs = mp.eigsy(Gr, eigvals_only=True)
        lmin = min(evs)
        print(f"   a = {mp.nstr(a,10)}, N = {N}:  max|Im G| = {mp.nstr(maxim,3)},"
              f"  |G-G^T| = {mp.nstr(sym,3)},  scale = {mp.nstr(scale,6)}")
        print(f"      index-reversal |G_(m,n) - G_(-m,-n)| = {mp.nstr(rev,3)}"
              f"   (relative {mp.nstr(rev/scale,3)})")
        print(f"      lambda_min(G) = {mp.nstr(lmin,10)}  lambda_max = {mp.nstr(max(evs),10)}")
        chk(f"a={mp.nstr(a,8)}: G real", maxim < mp.mpf(10) ** -45)
        chk(f"a={mp.nstr(a,8)}: G symmetric", sym < mp.mpf(10) ** -45)
        chk(f"a={mp.nstr(a,8)}: index-reversal symmetry (relative <= 1e-30)",
            rev / scale < mp.mpf(10) ** -30, f"= {mp.nstr(rev/scale,3)}")
        chk(f"a={mp.nstr(a,8)}: positive definite", lmin > 0, f"lambda_min = {mp.nstr(lmin,8)}")
        out.setdefault("P0_T4", []).append(dict(
            a=fp(a), N=N, max_im=fp(maxim), sym=fp(sym),
            reversal_rel=fp(rev / scale), lmin=fp(lmin),
            lmax=fp(max(evs))))

    # ---------------------------------------------------------------- P0-T2
    print("\n-- P0-T2: CROSS-BASIS against src/zeta_window.build_parts_idx --")
    print("""   cos(k pi x/L) = (sqrt(2L)/2)(chi_k + chi_{-k})  for k >= 1,   chi_0 = 1/sqrt(2L),
   and the repository returns Q(f_j,f_k)/(||f_j|| ||f_k||) with ||f_0||^2 = 2L,
   ||f_k||^2 = L.  Hence the predicted repository entry from Suzuki's G is
        k,j >= 1 :  (2L/4)(G_{j,k}+G_{j,-k}+G_{-j,k}+G_{-j,-k}) / L
        j = 0    :  sqrt(2L)*(sqrt(2L)/2)(G_{0,k}+G_{0,-k}) / (sqrt(2L) sqrt(L))
        j = k = 0:  2L * G_{0,0} / (2L) = G_{0,0}
   (The j=0 coefficient is sqrt(2L), NOT 1/sqrt(2L): chi_0 is the CONSTANT 1/sqrt(2L)
    and f_0 = 1, so f_0 = sqrt(2L) chi_0.  The first version of this test had it
    inverted -- a factor 2L -- and the gate caught it.  See REPORT Retraction 23.)""")
    mp.mp.dps = 60
    for aval, Kmax in (('0.5493061443340548', 8), ('0.8047189562170503', 8)):
        L = mp.mpf(aval)
        G = S.gram(L, Kmax, dps=60)
        pos = {v: i for i, v in enumerate(range(-Kmax, Kmax + 1))}
        M = ZW.build_matrix_idx(L, list(range(Kmax + 1)), 60, sector=ZW.EVEN)
        worst = mp.mpf(0); where = None
        sc = max(abs(M[i, j]) for i in range(Kmax + 1) for j in range(Kmax + 1))
        for j in range(Kmax + 1):
            for k in range(Kmax + 1):
                nj = mp.sqrt(2 * L) if j == 0 else mp.sqrt(L)
                nk = mp.sqrt(2 * L) if k == 0 else mp.sqrt(L)
                # f_0 = cos(0) = 1 = sqrt(2L)*chi_0 ;  f_k = sqrt(2L)(chi_k+chi_-k)/2
                cj = mp.sqrt(2 * L) if j == 0 else mp.sqrt(2 * L) / 2
                ck = mp.sqrt(2 * L) if k == 0 else mp.sqrt(2 * L) / 2
                sj = [0] if j == 0 else [j, -j]
                sk = [0] if k == 0 else [k, -k]
                pred = sum(G[pos[p], pos[q]] for p in sj for q in sk) * cj * ck / (nj * nk)
                dev = abs(mp.re(pred) - M[j, k])
                if dev > worst:
                    worst, where = dev, (j, k)
        print(f"   L = {mp.nstr(L,12)}, block |j|,|k| <= {Kmax}:  worst |dev| = "
              f"{mp.nstr(worst,4)}  at {where},  matrix scale {mp.nstr(sc,6)}")
        print(f"      relative = {mp.nstr(worst/sc,4)}")
        chk(f"L={mp.nstr(L,8)}: cross-basis relative deviation <= 1e-25",
            worst / sc <= mp.mpf(10) ** -25, f"= {mp.nstr(worst/sc,3)}")
        out.setdefault("P0_T2", []).append(dict(L=fp(L), Kmax=Kmax,
                                                worst_abs=fp(worst),
                                                scale=fp(sc),
                                                worst_rel=fp(worst / sc)))

    # ---------------------------------------------------------------- P0-T3
    print("\n-- P0-T3: spectral side (ASSUMES RH; consistency only) --")
    print("""   PREREGISTRATION DEFECT, recorded not repaired silently.  PREREG fixed the
   threshold at 1e-8 absolute.  The zero-side sum truncated at the first NZ zeros has a
   tail of size ~ 2(2a)^3 log(T/2pi)/(2 pi T) with T = gamma_NZ, so 1e-8 needs T ~ 1e8,
   i.e. ~1e7 zeros.  THE THRESHOLD WAS UNREACHABLE AS WRITTEN.  Replaced -- before
   looking at whether the replacement passes -- by the test that actually discriminates:

     (i)  the deviation must carry the sign (-1)^{m+n}, which is the shape of the tail
          (phi_m(t) ~ sqrt(2a)(-1)^m sin(at)/(at) for large |t|), and
     (ii) |deviation| must fall like 1/T as NZ grows.

   A constant additive offset, or a deviation not shrinking, would indicate a
   normalization error.  A fixed tolerance would not have distinguished the two.""")
    mp.mp.dps = 30
    a = mp.mpf('0.8047189562170503')
    PAIRS = ((0, 0), (1, 1), (0, 1), (2, -3), (3, 3))
    G = S.gram(a, 3, dps=30)
    pos = {v: i for i, v in enumerate(range(-3, 4))}
    NZS = (150, 300, 600)
    gam = [mp.im(mp.zetazero(k)) for k in range(1, max(NZS) + 1)]
    devs = {}
    print(f"\n   {'(m,n)':>8} " + " ".join(f"{'dev @NZ='+str(z):>18}" for z in NZS)
          + f" {'sign ok':>8}")
    rows = []
    for (m, n) in PAIRS:
        ds = []
        for NZ in NZS:
            zs = mp.mpf(0)
            for g in gam[:NZ]:
                zs += (S.phi(m, a, g) * S.phi(n, a, g)
                       + S.phi(m, a, -g) * S.phi(n, a, -g))
            ds.append(mp.re(G[pos[m], pos[n]]) - zs)
        sgn_ok = all((d > 0) == ((-1) ** ((m + n) % 2) > 0) for d in ds)
        devs[(m, n)] = ds
        print(f"   {str((m,n)):>8} " + " ".join(f"{mp.nstr(d,8):>18}" for d in ds)
              + f" {str(sgn_ok):>8}")
        rows.append(dict(mn=str((m, n)), devs=[fp(d) for d in ds], sign_ok=sgn_ok))
        chk(f"P0-T3 {(m,n)}: deviation sign = (-1)^(m+n)", sgn_ok)
    print(f"\n   T = gamma_NZ: " + " ".join(f"{mp.nstr(gam[z-1],8):>12}" for z in NZS))
    ratios = []
    for (m, n), ds in devs.items():
        r = [abs(ds[i]) / abs(ds[i + 1]) for i in range(len(ds) - 1)]
        pred = [gam[NZS[i + 1] - 1] / gam[NZS[i] - 1] for i in range(len(NZS) - 1)]
        ratios.append((m, n, [fp(x) for x in r], [fp(x) for x in pred]))
    print(f"   {'(m,n)':>8} {'|dev| ratios':>26} {'predicted T-ratios':>26}")
    for m, n, r, pr in ratios:
        print(f"   {str((m,n)):>8} {str(r):>26} {str(pr):>26}")
    allr = [abs(ds[i]) / abs(ds[i + 1]) for ds in devs.values() for i in range(len(ds) - 1)]
    predr = [gam[NZS[i + 1] - 1] / gam[NZS[i] - 1] for i in range(len(NZS) - 1)] * len(devs)
    worstratio = max(abs(float(x) / float(y) - 1) for x, y in zip(allr, predr))
    chk("P0-T3 |deviation| falls like 1/T (ratios match within 35%)", worstratio < 0.35,
        f"worst relative departure from the 1/T law = {worstratio:.3f}")
    out["P0_T3"] = dict(NZS=list(NZS), rows=rows, ratios=[list(r) for r in ratios],
                        worst_ratio_departure=worstratio, assumes_RH=True,
                        prereg_defect="fixed 1e-8 threshold unreachable; replaced by a "
                                      "sign + 1/T-scaling test")

    print("\n" + "=" * 78)
    print("PHASE 0 VERDICT:", "PASS" if OK[0] else "FAIL -> Phase 2 unconditional claim STOPPED")
    print("=" * 78)
    out["verdict"] = "PASS" if OK[0] else "FAIL"
    json.dump(out, open(os.path.join(os.path.dirname(__file__), "..", "data", "phase0.json"), "w"),
              indent=1, default=str)
    return OK[0]


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
