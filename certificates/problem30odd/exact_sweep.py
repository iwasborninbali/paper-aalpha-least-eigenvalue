#!/usr/bin/env python3
"""Exact finite sweep, complementing (not replacing) the all-n proof.

n=5: every one of 2^10 labelled graphs, constructed from edge masks.
n=7: every one of the 1044 isomorphism representatives in NetworkX atlas.
For every non-extremizer, an INTEGER Sylvester/Bareiss certificate proves
A_alpha(G) - qI positive definite, where q is a rigorously computed rational
upper bound on the proposed radical minimum (hence q >= that minimum).
For the extremizer the exact two-cell characteristic equation is checked by
certificate_checker.py. n=9 is certified for all graphs by REPORT.md's proof,
not by this enumerator; use --graph6 path to audit external n=9 graph6 data.
"""
import argparse
import hashlib
import math
from fractions import Fraction

ALPHAS = ('1/1000', '1/100', '1/20', '1/10', '1/8', '1/6', '1/5',
          '1/4', '3/10', '1/3', '2/5', '9/20', '49/100', '499/1000')
SCALE = 10**12

def pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

def bipartite_masks(n):
    p, q = n // 2, (n + 1) // 2
    pp = pairs(n)
    out = set()
    for subset in range(1 << n):
        if subset.bit_count() != p:
            continue
        out.add(sum(1 << k for k, (i, j) in enumerate(pp)
                    if bool(subset & (1 << i)) != bool(subset & (1 << j))))
    return out

def matrix(n, mask, a, b, t):
    es = pairs(n)
    # 2*b*SCALE * (A_{a/b} - q I), with
    # q = (a*n*SCALE - t)/(2*b*SCALE), t=floor(SCALE*sqrt(R)).
    M = [[0] * n for _ in range(n)]
    d = [0] * n
    for k, (u, v) in enumerate(es):
        if mask & (1 << k):
            d[u] += 1
            d[v] += 1
            M[u][v] = M[v][u] = 2 * SCALE * (b - a)
    for i in range(n):
        M[i][i] = 2 * SCALE * a * d[i] - (a * n * SCALE - t)
    return M

def positive_definite_bareiss(M):
    """Exact integer Sylvester test; each pivot is a leading minor.

    Symmetric fraction-free elimination is verified to divide exactly at
    every stage. No floating point is involved.
    """
    n = len(M)
    M = [row[:] for row in M]
    prev = 1
    for k in range(n - 1):
        pivot = M[k][k]
        if pivot <= 0:
            return False
        for i in range(k + 1, n):
            for j in range(i, n):
                numerator = pivot * M[i][j] - M[i][k] * M[k][j]
                val, rem = divmod(numerator, prev)
                if rem:
                    raise AssertionError('Bareiss division not exact')
                M[i][j] = M[j][i] = val
        prev = pivot
    return M[n - 1][n - 1] > 0

def encode_mask(g, n):
    return sum(1 << k for k, (u, v) in enumerate(pairs(n)) if g.has_edge(u, v))

def masks(n, external_graph6=None):
    if external_graph6:
        import networkx as nx
        with open(external_graph6, 'rb') as f:
            for line in f:
                if not line.strip() or line.startswith(b'>>'):
                    continue
                g = nx.from_graph6_bytes(line.strip())
                if len(g) != n:
                    raise AssertionError('Wrong graph order')
                yield encode_mask(g, n)
    elif n == 5:
        yield from range(1 << (n * (n - 1) // 2))
    elif n == 7:
        import networkx as nx
        for g in nx.graph_atlas_g():
            if len(g) == n:
                yield encode_mask(g, n)
    else:
        raise ValueError('n=9 requires --graph6 path; theorem covers n=9 otherwise')

def run(n, external_graph6=None):
    all_masks = tuple(masks(n, external_graph6))
    target_masks = bipartite_masks(n)
    count_target = sum(m in target_masks for m in all_masks)
    count_expected = {5: 1024, 7: 1044, 9: 274668}[n]
    assert len(all_masks) == count_expected, (len(all_masks), count_expected)
    assert len(set(all_masks)) == len(all_masks), 'Repeated labelled mask'
    assert count_target == (math.comb(n, n // 2) if n == 5 else 1), count_target
    fingerprint = hashlib.sha256(''.join(f'{m}\n' for m in all_masks).encode()).hexdigest()
    print(f'n={n} entries={len(all_masks)} balanced_K={count_target} '
          f'raw_mask_sha256={fingerprint}', flush=True)
    for a_text in ALPHAS:
        f = Fraction(a_text)
        a, b = f.numerator, f.denominator
        assert 0 < 2 * a < b
        R = a * a * n * n + (n * n - 1) * b * (b - 2 * a)
        root_floor = math.isqrt(R * SCALE * SCALE)
        # For non-square R, q > true minimum, hence PD at q is stronger.
        assert root_floor**2 <= R * SCALE**2 < (root_floor+1)**2
        failures = []
        for mask in all_masks:
            if mask in target_masks:
                continue
            if not positive_definite_bareiss(matrix(n, mask, a, b, root_floor)):
                failures.append(mask)
                if len(failures) >= 5:
                    break
        if failures:
            raise AssertionError(f'n={n} alpha={f}, nonK fails PD at q: {failures}')
        q = Fraction(a*n*SCALE - root_floor, 2*b*SCALE)
        print(f'alpha={f} radicand={R} exact=( {a*n} - sqrt({R}) )/{2*b} '
              f'nonK_PD={len(all_masks)-count_target} '
              f'upper_bound={q} upper_decimal={float(q):.12f} PASS', flush=True)
    print(f'PASS n={n} every listed graph and every listed rational alpha', flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, choices=(5, 7, 9), required=True)
    parser.add_argument('--graph6', default=None, help='external complete nonisomorphic graph6 file')
    args = parser.parse_args()
    run(args.n, args.graph6)
