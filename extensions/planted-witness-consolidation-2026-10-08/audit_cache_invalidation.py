#!/usr/bin/env python3
r"""Does the cache key actually invalidate?

A stale cache is the classic way a 'reverification' silently reproduces the number it was
supposed to check.  For each field of the key, this test mutates ONE field of a COPY of the
cache (the original is never touched -- the test runs in a scratch directory) and requires
that verify_certificate_consolidated.py REFUSES the cached value and recomputes.

Detection: Q_zeta is replaced by a sentinel-raising stub.  Sentinel raised  = cache rejected,
recompute attempted (what we want for a mutated key).  No sentinel = cache accepted.
"""
import importlib.util, json, os, shutil, sys, tempfile
from fractions import Fraction as Q

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'verify_certificate_consolidated.py')
CACHE = os.path.join(HERE, 'regenerated_qzeta_cache.json')


class Recomputed(Exception):
    pass


def run(mutate, label):
    """Returns True if the cache was ACCEPTED, False if it was rejected (recompute attempted)."""
    d = tempfile.mkdtemp(prefix='cacheaudit_')
    try:
        shutil.copy(SRC, d)
        cj = json.load(open(CACHE))
        mutate(cj)
        json.dump(cj, open(os.path.join(d, 'regenerated_qzeta_cache.json'), 'w'))
        cwd = os.getcwd()
        os.chdir(d)
        try:
            spec = importlib.util.spec_from_file_location(
                'vc_' + label, os.path.join(d, 'verify_certificate_consolidated.py'))
            m = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = m
            spec.loader.exec_module(m)

            def stub():
                raise Recomputed()
            m.Q_zeta = stub
            devnull = open(os.devnull, 'w')
            so, se = sys.stdout, sys.stderr
            sys.stdout = sys.stderr = devnull
            try:
                m.main()
                return True            # ran to the end on the cached value
            except Recomputed:
                return False           # key mismatch -> recompute
            finally:
                sys.stdout, sys.stderr = so, se
                devnull.close()
        finally:
            os.chdir(cwd)
    finally:
        shutil.rmtree(d, ignore_errors=True)


fails = []
print('control: the unmutated cache must be ACCEPTED')
acc = run(lambda cj: None, 'control')
print(f'  {"ok   " if acc else "FAIL "} unmutated cache accepted = {acc}')
if not acc:
    fails.append('control')

MUTATIONS = [
    ('coefficient x_0 changed by 1 in the last digit',
     lambda cj: cj['key'].__setitem__('x', [cj['key']['x'][0] + 1] + cj['key']['x'][1:])),
    ('coefficient denominator changed',
     lambda cj: cj['key'].__setitem__('xd', cj['key']['xd'] * 10)),
    ('window half-length L changed',
     lambda cj: cj['key'].__setitem__('L', '7/9')),
    ('detection height gamma changed',
     lambda cj: cj['key'].__setitem__('gamma', '70')),
    ('archimedean truncation M changed',
     lambda cj: cj['key'].__setitem__('M', cj['key']['M'] - 1)),
    ('rounding grid changed',
     lambda cj: cj['key'].__setitem__('grid', cj['key']['grid'] - 1)),
    ('basis dimension changed',
     lambda cj: cj['key'].__setitem__('ndim', 8)),
    ('basis source tag changed',
     lambda cj: cj['key'].__setitem__('src', 'odd_d sin(k pi x/L)/sqrt(L)')),
    ('key field deleted entirely',
     lambda cj: cj['key'].pop('gamma')),
    ('key block absent',
     lambda cj: cj.pop('key')),
]
print('\nmutations: each must be REJECTED (recompute attempted)')
for i, (label, mut) in enumerate(MUTATIONS):
    acc = run(mut, f'm{i}')
    ok = not acc
    print(f'  {"ok   " if ok else "FAIL "} {label}: {"rejected" if ok else "ACCEPTED STALE VALUE"}')
    if not ok:
        fails.append(label)

print("\ncontrol: a cache whose KEY matches but whose NUMBERS are corrupted")
print("         (Q_zeta forced to [-1,-1/2]).  The key is an input fingerprint, not a")
print("         checksum of the output, so the cache CHECK cannot catch this; the question")
print("         is whether anything downstream does.")


def corrupt(cj):
    cj['Q']['lo'] = str(Q(-1))
    cj['Q']['hi'] = str(Q(-1, 2))


try:
    acc = run(corrupt, 'corrupt')
    print(f'  FAIL  corrupted cache ran to completion (accepted = {acc})')
    fails.append('corrupted cache not caught')
except Exception as e:
    cls = type(e).__name__
    if 'CertificationFailure' in cls:
        print(f'  ok    corrupted cache was caught downstream by the GATES: {e}')
    else:
        print(f'  FAIL  corrupted cache raised an unexpected {cls}: {e}')
        fails.append('corrupted cache wrong failure mode')

print()
if fails:
    print(f'CACHE AUDIT FAILED: {fails}')
    sys.exit(1)
print('CACHE AUDIT PASSED: every input change invalidates the cache.')
