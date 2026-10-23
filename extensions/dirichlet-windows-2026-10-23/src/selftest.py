#!/usr/bin/env python3
r"""Phase 1 self-test.  Everything here must pass before Gate A is run.

Six independent checks, each with an externally known answer:

  (1) character axioms: completely multiplicative, periodic mod q, chi(-1) = (-1)^a,
      chi(n) = 0 iff gcd(n,q) > 1, and PRIMITIVE (no smaller period);
  (2) Gauss sum tau(chi) = sqrt(q) for a = 0 and i sqrt(q) for a = 1 -- hence eps(chi) = 1;
  (3) L(1, chi_d) against the class-number formula (external, closed form);
  (4) the Hurwitz-zeta assembly of L(s,chi) against the Dirichlet series at s = 3;
  (5) arch_closed_a(., ., a=0) == the parent repository's arch_closed(., .) -- the
      regression that makes the archimedean generalisation a strict extension;
  (6) Im Z_chi(t) at the precision floor: a numerical check of the functional equation
      and of eps(chi) = 1.

Exits non-zero if any check fails.  The gates assert the inequalities themselves.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp                                                      # noqa: E402
import dirichlet_window as dw                                            # noqa: E402
from zeta_window import arch_closed, Csym_terms, EVEN_D, EVEN, ODD       # noqa: E402

FAIL = []


def check(name, ok, detail=''):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ''))
    if not ok:
        FAIL.append(name)


def main():
    mp.mp.dps = 40
    print("=" * 78)
    print("(1) character axioms")
    print("=" * 78)
    for d in dw.DISCRIMINANTS:
        q, a = dw.CHARS[d]['q'], dw.CHARS[d]['a']
        vals = dw.char_values(d)
        mult = all(dw.chi(d, m * n) == dw.chi(d, m) * dw.chi(d, n)
                   for m in range(1, 3 * q) for n in range(1, 3 * q))
        per = all(dw.chi(d, n) == dw.chi(d, n + q) for n in range(1, 4 * q))
        par = dw.chi(d, q - 1) == (-1) ** a
        from math import gcd
        zer = all((dw.chi(d, n) == 0) == (gcd(n, q) > 1) for n in range(1, 4 * q))
        # primitivity: no proper divisor of q is a period on the units
        prim = True
        for m in range(1, q):
            if q % m == 0 and m < q:
                if all(dw.chi(d, n) == dw.chi(d, n + m)
                       for n in range(1, 3 * q) if gcd(n, q) == 1):
                    prim = False
        check(f"d={d:>3d} q={q} a={a} values={vals}",
              mult and per and par and zer and prim,
              f"mult={mult} per={per} parity={par} support={zer} primitive={prim}")

    print("=" * 78)
    print("(2) Gauss sum tau(chi) = i^a sqrt(q), hence eps(chi) = 1")
    print("=" * 78)
    for d in dw.DISCRIMINANTS:
        q, a = dw.CHARS[d]['q'], dw.CHARS[d]['a']
        tau = sum(dw.chi(d, r) * mp.e ** (2j * mp.pi * r / q) for r in range(q))
        want = (1j ** a) * mp.sqrt(q)
        eps = tau / ((1j ** a) * mp.sqrt(q))
        check(f"d={d:>3d} tau={mp.nstr(tau, 12)}", abs(tau - want) < mp.mpf(10) ** -30,
              f"eps = {mp.nstr(eps, 12)}")

    print("=" * 78)
    print("(3) L(1, chi_d) vs the class-number formula")
    print("=" * 78)
    # h = 1 for all five discriminants; w = 6, 4, 2, 2, 2
    closed = {
        -3: mp.pi / (3 * mp.sqrt(3)),
        -4: mp.pi / 4,
        -7: mp.pi / mp.sqrt(7),
         5: 2 * mp.log((1 + mp.sqrt(5)) / 2) / mp.sqrt(5),
         8: 2 * mp.log(1 + mp.sqrt(2)) / mp.sqrt(8),
    }
    for d in dw.DISCRIMINANTS:
        got = dw.L_chi(1, d)
        want = closed[d]
        rel = abs(got - want) / abs(want)
        check(f"d={d:>3d} L(1,chi) = {mp.nstr(got, 20)}", rel < mp.mpf(10) ** -30,
              f"closed form {mp.nstr(want, 20)}  rel.err {mp.nstr(rel, 4)}")

    print("=" * 78)
    print("(4) Hurwitz assembly vs EXTERNAL closed forms and a bounded-tail series")
    print("=" * 78)
    # (4a) Dirichlet beta: L(3, chi_-4) = pi^3/32, exactly.
    got, want = dw.L_chi(3, -4), mp.pi ** 3 / 32
    check("L(3, chi_-4) = pi^3/32", abs(got - want) / want < mp.mpf(10) ** -32,
          f"{mp.nstr(got, 24)} vs {mp.nstr(want, 24)}")
    # (4b) L(0, chi_d) = 2h(d)/w for d < 0 (RATIONAL), and 0 for d > 0.
    L0 = {-3: mp.mpf(1) / 3, -4: mp.mpf(1) / 2, -7: mp.mpf(1),
          5: mp.mpf(0), 8: mp.mpf(0)}
    for d in dw.DISCRIMINANTS:
        got = dw.L_chi(0, d)
        want = L0[d]
        ok = (abs(got - want) < mp.mpf(10) ** -30 if want == 0
              else abs(got - want) / want < mp.mpf(10) ** -30)
        check(f"d={d:>3d} L(0,chi) = 2h/w (d<0) or 0 (d>0)", ok,
              f"{mp.nstr(got, 20)} vs exact {mp.nstr(want, 20)}")
    # (4c) truncated Dirichlet series at s = 6 with a proved tail bound.
    #      |sum_{n>N} chi(n) n^{-6}| <= q * (N^{-6} + 6 int_N^inf x^{-7} dx) = 2 q N^{-6},
    #      by partial summation and |sum_{n=M+1}^{M+k} chi(n)| <= q.
    N = 4000
    for d in dw.DISCRIMINANTS:
        q = dw.CHARS[d]['q']
        ser = mp.fsum(dw.chi(d, n) / mp.mpf(n) ** 6 for n in range(1, N + 1))
        got = dw.L_chi(6, d)
        tail = 2 * q * mp.mpf(N) ** -6
        check(f"d={d:>3d} L(6,chi) vs series(N={N}), tail bound {mp.nstr(tail, 3)}",
              abs(got - ser) < 10 * tail,
              f"|diff| = {mp.nstr(abs(got - ser), 4)}")
    # (4d) functional equation Lambda(s) = Lambda(1-s) off the critical line.
    #      This tests the conductor q, the parity a and the root number simultaneously.
    for d in dw.DISCRIMINANTS:
        worst = mp.mpf(0)
        for z in (mp.mpc('0.3', '0.7'), mp.mpc('2.1', '-1.4'), mp.mpc('-0.6', '3.2')):
            A, B = dw.Lambda_chi(z, d), dw.Lambda_chi(1 - z, d)
            worst = max(worst, abs(A - B) / max(abs(A), mp.mpf(1)))
        check(f"d={d:>3d} Lambda(s) = Lambda(1-s) at 3 off-line points",
              worst < mp.mpf(10) ** -28, f"max rel.err {mp.nstr(worst, 6)}")

    print("=" * 78)
    print("(5) arch_closed_a(a=0) == parent arch_closed   [strict-extension regression]")
    print("=" * 78)
    worst = mp.mpf(0)
    for sector in (EVEN_D, EVEN, ODD):
        for Lx in ('0.5', '0.8', '1.3'):
            L = mp.mpf(Lx)
            for j in range(4):
                for k in range(j, 4):
                    tm = Csym_terms(j, k, L, sector)
                    v0 = arch_closed(tm, L)
                    v1 = dw.arch_closed_a(tm, L, 0)
                    worst = max(worst, abs(v0 - v1))
    check("max |arch_closed_a(a=0) - arch_closed| over 3 sectors x 3 L x 10 pairs",
          worst < mp.mpf(10) ** -32, f"max abs diff = {mp.nstr(worst, 6)}")
    # and a=1 must DIFFER (otherwise the parity never entered)
    tm = Csym_terms(0, 0, mp.mpf('0.8'), EVEN_D)
    d0 = dw.arch_closed_a(tm, mp.mpf('0.8'), 0)
    d1 = dw.arch_closed_a(tm, mp.mpf('0.8'), 1)
    check("a=1 differs from a=0 (parity actually enters)", abs(d1 - d0) > mp.mpf('1e-6'),
          f"a=0: {mp.nstr(d0, 12)}   a=1: {mp.nstr(d1, 12)}   diff {mp.nstr(d1 - d0, 6)}")

    print("=" * 78)
    print("(6) Im Z_chi(t) at the precision floor  [functional equation / root number]")
    print("=" * 78)
    for d in dw.DISCRIMINANTS:
        worst = mp.mpf(0)
        for t in ('1.0', '5.5', '12.3', '30.7'):
            re, im = dw.Z_chi(mp.mpf(t), d)
            worst = max(worst, abs(im) / max(abs(re), mp.mpf(1)))
        check(f"d={d:>3d} max |Im Z| / max(|Re Z|,1)", worst < mp.mpf(10) ** -30,
              f"{mp.nstr(worst, 6)}")

    print()
    if FAIL:
        print(f"FAILED: {len(FAIL)} check(s): {FAIL}")
        return 1
    print("ALL PHASE-1 SELF-TESTS PASSED")
    return 0


if __name__ == '__main__':
    sys.exit(main())
