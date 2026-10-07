r"""Small-height near-relations ("commas") among the active primes.

A comma is a rational a/b, gcd(a,b)=1, with a,b built only from the active primes and
max(a,b) <= height, whose DETUNING |log(a/b)| is small. These are the near-
commensurabilities that make the real prime comb different from a random-frequency one:
whenever  sum_p e_p log p ~ 0,  the corresponding phases on the torus are not independent
over long stretches of t.

Two scales matter and they are NOT the same:

  detuning           d = |sum_p e_p log p|          -- how close the relation is;
  preservation scale s* = d / sqrt(sum_p e_p^2 / 3) -- the surrogate spread below which a
                          U(-s,s) perturbation of each log p typically does NOT dominate d
                          (the perturbation of sum e_p log p has SD s sqrt(sum e_p^2/3)).

A surrogate at spread s destroys the relation when s >> s*, preserves it when s << s*.
So the dose-response markers belong at s*, while the raw detunings d are what the brief
asks to mark; we report both.
"""
import numpy as np
from math import gcd


def smooth_numbers(primes, height):
    vals = {1: {}}
    frontier = [1]
    while frontier:
        nxt = []
        for v in frontier:
            for p in primes:
                w = v * p
                if w <= height and w not in vals:
                    e = dict(vals[v]); e[p] = e.get(p, 0) + 1
                    vals[w] = e
                    nxt.append(w)
        frontier = nxt
    return vals


def commas(primes, height=200, max_detune=0.30):
    """[(a, b, detuning, sum e^2, s*)] sorted by detuning."""
    vals = smooth_numbers(sorted(primes), height)
    ns = sorted(vals)
    out = []
    for i, a in enumerate(ns):
        for b in ns:
            if b >= a or gcd(a, b) != 1:
                continue
            d = abs(np.log(a / b))
            if 0 < d <= max_detune:
                e = dict(vals[a])
                for p, k in vals[b].items():
                    e[p] = e.get(p, 0) - k
                s2 = sum(k * k for k in e.values())
                out.append((a, b, float(d), int(s2), float(d / np.sqrt(s2 / 3.0))))
    out.sort(key=lambda r: r[2])
    return out
