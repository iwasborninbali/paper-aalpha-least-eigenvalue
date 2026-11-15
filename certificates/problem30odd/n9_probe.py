#!/usr/bin/env python3
"""Deterministic exact n=9 stress test; NOT exhaustive graph enumeration.

Checks every near-balanced K_{4,5} relabeling, each single-edge deletion,
each single within-part edge addition, and 4000 seeded random masks.
All non-K cases must pass the exact rational PD test for all 14 alphas.
The n=9 *exhaustive* certification is the proof in REPORT.md, not this probe.
"""
import math
import random
from fractions import Fraction
from exact_sweep import ALPHAS, SCALE, bipartite_masks, matrix, positive_definite_bareiss, pairs

n = 9
N = n*(n-1)//2
universe = (1 << N) - 1
K = bipartite_masks(n)
S = set(K)
for k in K:
    # Every one-edge deletion from each labelled complete near-balanced K.
    for j in range(N):
        bit = 1 << j
        if k & bit:
            S.add(k ^ bit)
        else:
            S.add(k | bit)
# 4000 seeded randomly chosen labelled graphs, independent of atlas data.
seed = 20261115
rng = random.Random(seed)
S.update(rng.getrandbits(N) for _ in range(4000))
sorted_masks = tuple(sorted(S))
print(f'n=9 labelled_candidates={len(sorted_masks)} balanced_K_in_set={sum(m in K for m in S)} '
      f'seed={seed} sample=K_relabelings + one_edge_edits + random_4000', flush=True)
for a_text in ALPHAS:
    f = Fraction(a_text)
    a, b = f.numerator, f.denominator
    R = a*a*n*n+(n*n-1)*b*(b-2*a)
    t = math.isqrt(R*SCALE*SCALE)
    if not all(positive_definite_bareiss(matrix(n, mask, a, b, t)) for mask in sorted_masks if mask not in K):
        raise AssertionError(f'n=9 alpha={a}/{b} candidate not proven PD')
    print(f'alpha={a}/{b} radicand={R} '
          f'nonK_PD={len(sorted_masks)-len(K)} PASS', flush=True)
print('PASS n9 stress test; NOTE not full 274668-class graph6 enumeration', flush=True)
