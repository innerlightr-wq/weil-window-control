r"""N1: torus measure.

Phases theta_p uniform and independent per PRIME; prime powers of p share theta_p, so the
contribution of p is the block

    block_p(theta) = sum_k a_{p^k} cos(k theta),    a_{p^k} = 2 Lambda(p^k)/sqrt(p^k).

P_L(theta) = sum_p block_p(theta_p).  The law of each block is the pushforward of the
uniform measure on the circle (computed exactly by fine sampling + histogram); the law of
P_L is their convolution, done by FFT.  pi_L(x) = Pr[P_L > x].
"""
import numpy as np


def prime_blocks(sym):
    """{p: (ks, coeffs)} grouping the active prime powers by their base prime."""
    out = {}
    for p, k, w in zip(sym.prime, sym.expo, sym.weights):
        out.setdefault(int(p), []).append((int(k), float(w)))
    return out


def block_values(ks_ws, n_theta=1 << 16):
    th = np.linspace(0.0, 2 * np.pi, n_theta, endpoint=False)
    v = np.zeros(n_theta)
    for k, w in ks_ws:
        v += w * np.cos(k * th)
    return v


def tail_pdf(sym, dx=0.002, n_theta=1 << 16):
    """Returns (grid x, pdf of P_L) by FFT-convolving the per-prime block laws."""
    blocks = prime_blocks(sym)
    A = sym.A
    # common grid, wide enough for the full support [-A, A]
    nx = int(2 * np.ceil(A / dx)) + 8
    nx = 1 << int(np.ceil(np.log2(nx * 2)))        # pad for linear convolution
    x = (np.arange(nx) - nx // 2) * dx
    acc = None
    for p, ks_ws in blocks.items():
        v = block_values(ks_ws, n_theta)
        h, _ = np.histogram(v, bins=np.append(x - dx / 2, x[-1] + dx / 2))
        h = h.astype(float) / h.sum()
        if acc is None:
            acc = h
        else:
            acc = np.real(np.fft.ifft(np.fft.fft(acc) * np.fft.fft(h)))
            acc = np.roll(acc, -nx // 2)
            acc = np.maximum(acc, 0.0)
            acc /= acc.sum()
    return x, acc


def tail_prob(x_grid, pdf):
    """pi(x) = Pr[P > x] as a callable on the same grid."""
    surv = np.concatenate([np.cumsum(pdf[::-1])[::-1][1:], [0.0]])
    def pi(xq):
        xq = np.asarray(xq, dtype=float)
        return np.interp(xq, x_grid, surv, left=1.0, right=0.0)
    return pi


def monte_carlo(sym, n=10_000_000, seed=0):
    rng = np.random.default_rng(seed)
    blocks = prime_blocks(sym)
    tot = np.zeros(n)
    for p, ks_ws in blocks.items():
        th = rng.uniform(0, 2 * np.pi, n)
        b = np.zeros(n)
        for k, w in ks_ws:
            b += w * np.cos(k * th)
        tot += b
    return tot


def predict(sym, beta, pi, t_hi, n_t=200_000):
    """Expected dip measure density pi_L(H(t)-beta); returns (t grid, density, cum tail)."""
    t = np.geomspace(2 * np.pi, t_hi, n_t)
    H = sym.H(t)
    dens = pi(H - beta)
    # expected measure beyond each t  (integrate the density from t to t_hi)
    dt = np.diff(t, prepend=t[0])
    cum_from_end = np.cumsum((dens * dt)[::-1])[::-1]
    return t, dens, cum_from_end


def predicted_t0(sym, beta, pi, t_hi, median_dip_len, n_t=200_000):
    t, dens, tail = predict(sym, beta, pi, t_hi, n_t)
    idx = np.flatnonzero(tail < median_dip_len)
    return float(t[idx[0]]) if idx.size else float(t[-1])
