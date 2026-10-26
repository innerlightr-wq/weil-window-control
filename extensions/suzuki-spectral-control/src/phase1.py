r"""PHASE 1 -- reproduce Suzuki section 6.  Tests P1-T1..T5 of PREREG.md.

G_{a,N} is built ONCE per (a, dps) at N = 18 and SLICED for N = 14, 16: the index sets
{-14..14} and {-16..16} are sub-blocks of {-18..18} at the same a, so the submatrices are
exactly G_{a,14} and G_{a,16}.  Gram matrices are cached to data/gram_*.json.

Match criterion (PREREG): agreement with Suzuki's printed digits.  Any mismatch is
investigated and logged, never tuned away.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__))
import mpmath as mp
import suzuki as S

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
XS = (3, 5, 7, 11)
NS = (14, 16, 18)
NMAX = 18

# Suzuki Table 1 (10 d.p.), transcribed from arXiv:2610.32704v1
TABLE1 = {
    "gamma": ["14.1347251417", "21.0220396388", "25.0108575801",
              "30.4248761259", "32.9350615877", "37.5861781588"],
    3:  ["14.1347274566", "21.0226378897", "25.0192189912",
         "30.6272491287", "33.6051319793", "38.6158405572"],
    5:  ["14.1347251417", "21.0220396388", "25.0108575802",
         "30.4248761361", "32.9350616846", "37.5861810892"],
    7:  ["14.1347251417", "21.0220396388", "25.0108575801",
         "30.4248761259", "32.9350615877", "37.5861781588"],
    11: ["14.1347251417", "21.0220396388", "25.0108575801",
         "30.4248761259", "32.9350615877", "37.5861781588"],
}
# Suzuki Table 2: d6^{(14)}, d6^{(16)}, d6^{(18)}, kappa2(G_{a,18})
TABLE2 = {
    3:  ("1.04", "1.03", "1.03", "4.70e7"),
    5:  ("5.21e-6", "3.00e-6", "2.93e-6", "1.85e17"),
    7:  ("1.11e-9", "1.52e-11", "5.43e-13", "8.24e25"),
    11: ("1.28e-9", "3.70e-15", "1.03e-17", "2.33e34"),
}


def gram_cached(X, dps, progress=True):
    path = os.path.join(DATA, f"gram_X{X}_N{NMAX}_dps{dps}.json")
    if os.path.exists(path):
        d = json.load(open(path))
        mp.mp.dps = dps
        n = 2 * NMAX + 1
        G = mp.zeros(n, n)
        for i in range(n):
            for j in range(n):
                G[i, j] = mp.mpf(d["G"][i][j])
        return G
    mp.mp.dps = dps
    a = mp.log(mp.mpf(X)) / 2
    t0 = time.time()
    if progress:
        print(f"      building G_(a,{NMAX}) for X={X} at dps {dps} ...", flush=True)
    G = S.gram(a, NMAX, dps=dps)
    n = G.rows
    json.dump({"X": X, "N": NMAX, "dps": dps,
               "a": mp.nstr(a, dps),
               "G": [[mp.nstr(mp.re(G[i, j]), dps) for j in range(n)] for i in range(n)]},
              open(path, "w"))
    if progress:
        print(f"      done in {time.time()-t0:.1f}s", flush=True)
    return G


def slice_gram(G, N):
    """G_{a,N} as the central (2N+1) sub-block of G_{a,NMAX}."""
    off = NMAX - N
    n = 2 * N + 1
    return mp.matrix([[G[off + i, off + j] for j in range(n)] for i in range(n)])


def _bisect(f, lo, hi, dps):
    """Plain bisection to the working precision.  mpmath's findroot(solver='bisect')
    gives up at its default maxsteps well short of 140 digits; the function here is
    cheap (one Gram-solve dot product), so an explicit loop is both simpler and
    reliable.  Runs until the bracket is below 10^-(dps-20) or 8 dps iterations."""
    lo, hi = mp.mpf(lo), mp.mpf(hi)
    flo = f(lo)
    tol = mp.mpf(10) ** (-(dps - 20))
    for _ in range(8 * dps):
        mid = (lo + hi) / 2
        if hi - lo < tol:
            break
        fm = f(mid)
        if fm == 0:
            return mid
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


def roots_on(f, lo, hi, npts, dps):
    """Sign changes on a uniform grid, then full-precision bisection."""
    lo, hi = mp.mpf(lo), mp.mpf(hi)
    xs = [lo + (hi - lo) * k / (npts - 1) for k in range(npts)]
    vs = [f(x) for x in xs]
    out = []
    for i in range(npts - 1):
        if vs[i] == 0:
            out.append(xs[i]); continue
        if (vs[i] < 0) != (vs[i + 1] < 0):
            out.append(_bisect(f, xs[i], xs[i + 1], dps))
    return out


def main():
    OK = [True]
    out = {}
    def chk(name, cond, extra=""):
        OK[0] &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name} {extra}")

    print("=" * 78)
    print("PHASE 1 -- reproduce Suzuki Section 6")
    print("=" * 78)

    # ---------------------------------------------------------- P1-T5: the J convention
    print("\n-- P1-T5: the J convention in (6.4) --")
    print("""   PREREGISTRATION DEFECT, recorded not repaired silently.  PREREG asserted that the
   brief's restatement  J(w,z) = phi(z)^T conj(G^-1 phi(w))  and "the paper's (6.4)"
   J(w,z) = phi(z)^T G^-1 phi(w)  were different functions, to be separated empirically.
   THERE IS NO SUCH DISCREPANCY.  The paper's (6.4), as TYPESET, carries an overline over
   the whole group G^-1 phi(w) -- i.e. it IS the conjugated form.  Both pdftotext and
   pypdf silently drop that overline; rendering page 15 as an image shows it plainly.
   The brief transcribed the paper correctly and my text extraction did not.
   What follows is therefore not a choice between the paper and the brief, but a
   measurement of how far the overline-dropped MISREADING diverges.  See Retraction 25.""")
    dps = 140
    X = 5
    G18 = gram_cached(X, dps)
    mp.mp.dps = dps
    a = mp.log(mp.mpf(X)) / 2
    tab = {}
    for drop in (False, True):
        f = S.W_N_factory(G18, a, NMAX, drop_overline=drop)
        rr = [r for r in roots_on(f, mp.mpf('0.5'), 40, 4001, dps) if r > 1]
        tab[drop] = rr
    print(f"\n   {'j':>2} {'Suzuki Table 1 (X=5)':>22} {'(6.4) with overline':>22}"
          f" {'overline dropped':>22}")
    nmatch = {False: 0, True: 0}
    for j in range(6):
        cells = []
        for drop in (False, True):
            rr = tab[drop]
            v = mp.nstr(rr[j], 12) if j < len(rr) else "--"
            cells.append(v)
            if v != "--" and abs(rr[j] - mp.mpf(TABLE1[5][j])) <= mp.mpf('5e-10'):
                nmatch[drop] += 1
        print(f"   {j+1:>2} {TABLE1[5][j]:>22} {cells[0]:>22} {cells[1]:>22}")
    print(f"\n   of six Table-1 values reproduced to 10 d.p.:"
          f"  with overline {nmatch[False]}/6,  dropped {nmatch[True]}/6")
    chk("P1-T5: the typeset (6.4) (with overline) reproduces all six Table-1 values",
        nmatch[False] == 6)
    chk("P1-T5: the overline-dropped misreading does NOT", nmatch[True] < 6,
        f"it agrees on the first {nmatch[True]} and then diverges")
    DROP = False
    out["P1_T5"] = dict(with_overline=nmatch[False], dropped=nmatch[True],
                        adopted="typeset (6.4), overline present",
                        prereg_defect="PREREG asserted a paper/brief discrepancy that does "
                                      "not exist; it was a pdftotext+pypdf extraction "
                                      "artifact (dropped overline)")

    # ---------------------------------------------------------- P1-T1 / P1-T2 / P1-T3 / P1-T4
    res = {}
    for X in XS:
        print(f"\n{'='*78}\n   X = {X},  a = (1/2) log X\n{'='*78}")
        for dps in (140, 180):
            G18 = gram_cached(X, dps)
            mp.mp.dps = dps
            a = mp.log(mp.mpf(X)) / 2
            entry = res.setdefault((X, dps), {})
            # --- kappa2 and lambda_min at N = 18
            evs = mp.eigsy(G18, eigvals_only=True)
            lmin, lmax = min(evs), max(evs)
            entry["lmin"] = mp.nstr(lmin, 12); entry["lmax"] = mp.nstr(lmax, 12)
            entry["kappa2"] = mp.nstr(lmax / lmin, 6)
            # --- roots for each N
            for N in NS:
                GN = slice_gram(G18, N)
                f = S.W_N_factory(GN, a, N, drop_overline=DROP)
                rr = roots_on(f, mp.mpf('0.5'), 40, 4001, dps)
                entry[f"roots_N{N}"] = [mp.nstr(r, 14) for r in rr]
            if dps == 140:
                print(f"   dps 140: lambda_min = {mp.nstr(lmin,8)}, lambda_max = "
                      f"{mp.nstr(lmax,8)}, kappa2 = {mp.nstr(lmax/lmin,6)}"
                      f"   (Suzuki: {TABLE2[X][3]})")
                for N in NS:
                    print(f"      N={N}: {len(entry[f'roots_N{N}'])} roots in [0.5,40]")
        # --- precision doubling on the roots
        e140, e180 = res[(X, 140)], res[(X, 180)]
        worst = mp.mpf(0)
        for N in NS:
            r1, r2 = e140[f"roots_N{N}"], e180[f"roots_N{N}"]
            if len(r1) == len(r2):
                worst = max([worst] + [abs(mp.mpf(u) - mp.mpf(v)) for u, v in zip(r1, r2)])
            else:
                worst = mp.inf
        print(f"   precision doubling 140 -> 180, worst root shift: {mp.nstr(worst,4)}")
        chk(f"X={X}: root counts and values stable under precision doubling",
            worst < mp.mpf('1e-25'), f"worst shift {mp.nstr(worst,3)}")
        res[(X, 'dbl')] = mp.nstr(worst, 6)

    # ---------------------------------------------------------- Table 1
    print(f"\n{'='*78}\n   TABLE 1 -- first six positive roots of W_18\n{'='*78}")
    print(f"   {'j':>2} {'gamma_j':>16} " + " ".join(f"{'X='+str(X):>18}" for X in XS))
    t1ok = True
    for j in range(6):
        row = []
        for X in XS:
            rr = [mp.mpf(r) for r in res[(X, 140)]["roots_N18"]]
            rr = [r for r in rr if r > 1]
            row.append(mp.nstr(rr[j], 12) if j < len(rr) else "--")
        print(f"   {j+1:>2} {TABLE1['gamma'][j]:>16} " + " ".join(f"{v:>18}" for v in row))
        for X, v in zip(XS, row):
            if v == "--":
                t1ok = False; continue
            if abs(mp.mpf(v) - mp.mpf(TABLE1[X][j])) > mp.mpf('5e-10'):
                t1ok = False
    chk("P1-T1: Table 1 reproduced to all 10 printed decimal places", t1ok)
    out["P1_T1"] = dict(ok=t1ok, table={str(X): res[(X,140)]["roots_N18"][:8] for X in XS})

    # ---------------------------------------------------------- Table 2
    print(f"\n{'='*78}\n   TABLE 2 -- d6^(N) and kappa2(G_(a,18))\n{'='*78}")
    g6 = mp.mpf(TABLE1['gamma'][5])
    print(f"   {'X':>3} {'d6^(14)':>14} {'d6^(16)':>14} {'d6^(18)':>14} {'kappa2':>13}"
          f"   {'Suzuki d6/kappa2':>40}")
    t2rows = []
    for X in XS:
        e = res[(X, 140)]
        ds = []
        for N in NS:
            rr = [mp.mpf(r) for r in e[f"roots_N{N}"] if mp.mpf(r) > 1]
            ds.append(rr[5] - g6 if len(rr) > 5 else mp.nan)
        k2 = e["kappa2"]
        print(f"   {X:>3} " + " ".join(f"{mp.nstr(d,3):>14}" for d in ds)
              + f" {k2:>13}   {str(TABLE2[X]):>40}")
        t2rows.append(dict(X=X, d6=[mp.nstr(d, 6) for d in ds], kappa2=k2,
                           suzuki=TABLE2[X]))
    out["P1_T2"] = t2rows

    # ---------------------------------------------------------- P1-T3 forced lattice
    print(f"\n{'='*78}\n   P1-T3 -- the forced-zero lattice  W_N(k pi/a) = 0 for |k| > N\n{'='*78}")
    print("""   phi_n(k pi/a) = sqrt(2a) sincpi(k + n) = 0 whenever k + n != 0.  For |k| > N no
   n in [-N,N] has k + n = 0, so EVERY phi_n vanishes and W_N(k pi/a) = 0 identically.""")
    lat = []
    for X in XS:
        dps = 140
        G18 = gram_cached(X, dps); mp.mp.dps = dps
        a = mp.log(mp.mpf(X)) / 2
        for N in NS:
            GN = slice_gram(G18, N)
            f = S.W_N_factory(GN, a, N, drop_overline=DROP)
            scale = max(abs(f(mp.mpf(t))) for t in ('3.7', '11.3', '19.9', '28.1'))
            out_hi = [abs(f(k * mp.pi / a)) for k in (N + 1, N + 2, N + 5)]
            in_lo = [abs(f(k * mp.pi / a)) for k in (1, max(1, N // 2), N)]
            okl = all(v / scale < mp.mpf(10) ** -(dps - 25) for v in out_hi)
            oki = all(v / scale > mp.mpf(10) ** -8 for v in in_lo)
            lat.append(dict(X=X, N=N, scale=mp.nstr(scale, 6),
                            outside=[mp.nstr(v, 4) for v in out_hi],
                            inside=[mp.nstr(v, 4) for v in in_lo], ok=okl and oki))
            print(f"   X={X:>2} N={N}: |W_N| at k = {N+1},{N+2},{N+5} (|k|>N): "
                  + ", ".join(mp.nstr(v, 3) for v in out_hi)
                  + f"   |  at k = 1,{max(1,N//2)},{N} (|k|<=N): "
                  + ", ".join(mp.nstr(v, 3) for v in in_lo))
            chk(f"   X={X} N={N}: lattice vanishing / non-vanishing", okl and oki)
    out["P1_T3"] = lat

    # ---------------------------------------------------------- P1-T4 the 39.3044 artifact
    print(f"\n{'='*78}\n   P1-T4 -- Suzuki's extra zero of W_14 at 30 pi/log 11\n{'='*78}")
    dps = 140
    G18 = gram_cached(11, dps); mp.mp.dps = dps
    a = mp.log(mp.mpf(11)) / 2
    art = 30 * mp.pi / mp.log(mp.mpf(11))
    print(f"   30 pi/log 11 = {mp.nstr(art,14)}   (Suzuki prints 39.3044)")
    print(f"   15 pi/a      = {mp.nstr(15*mp.pi/a,14)}   -> the SAME point, lattice k = 15")
    chk("P1-T4: 30 pi/log 11 = 15 pi/a, i.e. lattice index k = 15 > N = 14",
        abs(art - 15 * mp.pi / a) < mp.mpf('1e-100'))
    f14 = S.W_N_factory(slice_gram(G18, 14), a, 14, drop_overline=DROP)
    f18 = S.W_N_factory(slice_gram(G18, 18), a, 18, drop_overline=DROP)
    sc14 = max(abs(f14(mp.mpf(t))) for t in ('3.7', '11.3', '19.9', '28.1'))
    print(f"   |W_14(39.3044)| = {mp.nstr(abs(f14(art)),4)}  (relative to scale "
          f"{mp.nstr(abs(f14(art))/sc14,4)})  -> a zero, as Suzuki reports")
    chk("P1-T4: W_14 vanishes there", abs(f14(art)) / sc14 < mp.mpf(10) ** -(dps - 25))
    k19 = 19 * mp.pi / a
    print(f"\n   Preregistered corollary: for N = 18 the first lattice point is k = 19,")
    print(f"   19 pi/a = 38 pi/log 11 = {mp.nstr(k19,12)}, which is OUTSIDE [0,40].")
    chk("P1-T4 corollary: the first N=18 lattice point lies above 40", k19 > 40,
        f"= {mp.nstr(k19,10)} > 40, so no such artifact can appear in Suzuki's N=18 table")
    out["P1_T4"] = dict(artifact=mp.nstr(art, 14), k=15,
                        W14_rel=mp.nstr(abs(f14(art)) / sc14, 6),
                        first_N18_lattice=mp.nstr(k19, 12), outside_range=bool(k19 > 40))

    print("\n" + "=" * 78)
    print("PHASE 1 VERDICT:", "PASS" if OK[0] else "see REPORT -- mismatches logged")
    print("=" * 78)
    out["verdict"] = "PASS" if OK[0] else "MISMATCH"
    out["res"] = {f"{k[0]}_{k[1]}": v for k, v in res.items()}
    json.dump(out, open(os.path.join(DATA, "phase1.json"), "w"), indent=1, default=str)
    return OK[0]


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
