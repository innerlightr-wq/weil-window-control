#!/usr/bin/env python3
r"""Driver: validation, the two witnesses, and the exact-quartet detection scan.

Usage:  python3 run_certify.py            (writes results.json)
"""
import json, math, sys, time
sys.path.insert(0, '.')
import certify_witness as cw

L = 0.8
M = 2000000                     # archimedean truncation; tail bound ~ 4e-6 per entry
DELTA_MAX = 0.5                 # a zero at 1/2 + delta is in the critical strip iff |delta| <= 1/2
ORDS = [('gamma_exact_14', 14.0, 'EXACT chosen ordinate, gamma = 14 exactly'),
        ('gamma_1', float(cw.GAMMA1_STR),
         'documented planted ordinate gamma_1; NOT rigorously enclosed in this run')]
NS = [4, 8]
NS_TREND = [16]


def witness_vectors(L, g, N):
    av = [cw.Fk(k, L, g) for k in range(N)]
    bv = [-cw.dFk(k, L, g) for k in range(N)]
    na2 = sum(t*t for t in av)
    nb2 = sum(t*t for t in bv)
    ab = sum(av[i]*bv[i] for i in range(N))
    A_vec = bv[:]
    B_vec = [bv[i] - (ab/na2)*av[i] for i in range(N)] if na2 > 0 else bv[:]
    return av, bv, na2, nb2, ab, A_vec, B_vec


def scan(x, L, g, Qz, Qe, dmax, n=200000):
    """Exact-quartet scan on (0, dmax].  Returns the first-negative crossing enclosure."""
    prev_d, prev_v = 0.0, Qz + Qe
    first = None
    worst = (0.0, prev_v)
    for i in range(1, n+1):
        d = dmax*i/n
        v = Qz + Qe + cw.exact_response(x, L, g, d)['delta_Q']
        if v < worst[1]:
            worst = (d, v)
        if first is None and v < 0:
            first = (prev_d, d, prev_v, v)
        prev_d, prev_v = d, v
    return first, worst, prev_v


def controls(x, L, g):
    r0 = cw.exact_response(x, L, g, 0.0)
    rp = cw.exact_response(x, L, g, 0.137)
    rm = cw.exact_response(x, L, g, -0.137)
    Cc = r0['c_inner']
    return dict(delta0_deltaQ=r0['delta_Q'], delta0_expected=-2*Cc*Cc + 4*Cc*Cc,
                delta0_residual=abs(r0['delta_Q'] - 2*Cc*Cc),
                even_in_delta_dev=abs(rp['delta_Q'] - rm['delta_Q']),
                c_inner=Cc)


out = dict(meta=dict(L=L, archimedean_truncation_M=M, delta_max=DELTA_MAX,
                     basis='even_d  phi_k = cos((2k+1)pi x/2L)/sqrt(L), Gram = I',
                     route='exact finite-displacement quartet (no Taylor remainder)',
                     generated=time.strftime('%Y-%m-%dT%H:%M:%S')),
           witnesses=[])

for ordkey, g, ordnote in ORDS:
    K = cw.K_exact(L, g)
    for N in NS + NS_TREND:
        A, E, prim = cw.weil_matrix_fast(L, N, M)
        av, bv, na2, nb2, ab, A_vec, B_vec = witness_vectors(L, g, N)
        Keff = nb2 - ab*ab/na2
        for name, vec, rule in (
                ('A_projected_r', A_vec, 'P_E r normalised; full F F\'\' retained'),
                ('B_nodal', B_vec, 'P_E r with the F(gamma)=0 constraint projected out')):
            n2 = sum(t*t for t in vec)
            x = [t/math.sqrt(n2) for t in vec]
            Fg = sum(x[k]*av[k] for k in range(N))
            Fp = -sum(x[k]*bv[k] for k in range(N))
            Fpp = sum(x[k]*cw.ddFk(k, L, g) for k in range(N))
            a_f = 4*(Fp*Fp + Fg*Fpp)
            Qz = cw.quad_form(A, x)
            Qe = cw.quad_form_err(E, x)
            b0 = Qz + 2*Fg*Fg
            first, worst, vend = scan(x, L, g, Qz, Qe, DELTA_MAX)
            # Taylor route, for comparison only
            tay = None
            if a_f > 0:
                d2 = (b0 + Qe)/a_f
                if d2 > 0:
                    d = math.sqrt(d2)
                    R = (16/3)*L**5*d**4*math.exp(2*L*d)
                    tay = dict(delta_naive=d, R_upper_there=R,
                               signal_there=a_f*d2, R_over_signal=R/(a_f*d2),
                               in_strip=bool(d <= 0.5))
            rec = dict(
                ordinate_key=ordkey, ordinate=g, ordinate_note=ordnote,
                ordinate_is_exact=(ordkey == 'gamma_exact_14'),
                trend_only=(N in NS_TREND),
                witness=name, selection_rule=rule, basis='even_d', dim=N,
                coeffs=x,
                K=K, norm_PEr2=nb2, K_eff=Keff,
                ratio_Keff_over_K=Keff/K, ratio_Keff_over_PEr2=Keff/nb2,
                ratio_PEr2_over_K=nb2/K,
                F_gamma=Fg, Fp_gamma=Fp, Fpp_gamma=Fpp,
                a_f=a_f, a_lower=a_f, four_Keff=4*Keff,
                Q_zeta=Qz, Q_zeta_tail_bound=Qe, b_0=b0, b_upper=b0 + Qe,
                exact_first_negative=(None if first is None else
                                      dict(delta_lo=first[0], delta_hi=first[1],
                                           Q_at_lo=first[2], Q_at_hi=first[3])),
                exact_min_on_strip=dict(delta=worst[0], Q_upper=worst[1]),
                exact_Q_at_delta_max=vend,
                taylor_route=tay,
                controls=controls(x, L, g))
            out['witnesses'].append(rec)
            tag = 'TREND' if N in NS_TREND else 'FROZEN'
            fn = rec['exact_first_negative']
            print(f"{ordkey:16s} N={N:2d} {name:14s} [{tag}] "
                  f"K_eff/K={Keff/K:.4f} Q_zeta={Qz:.6f}+-{Qe:.1e} a_f={a_f:.6f}")
            if fn:
                print(f"{'':20s}EXACT first negative at delta in "
                      f"[{fn['delta_lo']:.6f}, {fn['delta_hi']:.6f}]   "
                      f"min on strip {worst[1]:+.6e} at delta={worst[0]:.4f}")
            else:
                print(f"{'':20s}NO negative value on (0, {DELTA_MAX}]; "
                      f"min {worst[1]:+.6e} at delta={worst[0]:.4f}")

with open('results.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print("\nwrote results.json")
