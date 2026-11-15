#!/usr/bin/env python3
"""Exact rational separating-cut certificates for the runner-up, n=3..7.

No floating point enters any certificate comparison. NumPy is used only to
choose a rational t between candidate and next-nearest values; the ensuing
integer Sylvester tests and negative determinants establish the result.
"""
from fractions import Fraction
from math import floor
import networkx as nx
import numpy as np

ALPHAS = [Fraction(1, 10), Fraction(1, 4), Fraction(2, 5), Fraction(49, 100)]

def matrices(g, a, t):
    n = len(g)
    p, q = a.numerator, a.denominator
    scale = q * t.denominator
    B = [[-scale * t.numerator // t.denominator if i == j else 0 for j in range(n)] for i in range(n)]
    # scale=t.denominator*q: A_alpha - t I, all entries integral
    for i in range(n):
        B[i][i] += p * t.denominator * g.degree(i)
    for i, j in g.edges():
        B[i][j] = B[j][i] = (q - p) * t.denominator
    return B

def det_bareiss(B):
    B = [r[:] for r in B]
    n = len(B)
    sign, prev = 1, 1
    for k in range(n - 1):
        if B[k][k] == 0:
            m = next((i for i in range(k+1, n) if B[i][k]), None)
            if m is None:
                return 0
            B[k], B[m] = B[m], B[k]
            sign = -sign
        pivot = B[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                num = B[i][j] * pivot - B[i][k] * B[k][j]
                assert num % prev == 0
                B[i][j] = num // prev
            B[i][k] = 0
        prev = pivot
    return sign * B[-1][-1]

def positive_definite(B):
    return all(det_bareiss([row[:k] for row in B[:k]]) > 0 for k in range(1, len(B)+1))

def candidate(n):
    p, q = n//2, n-n//2
    if n == 3:
        g = nx.Graph();g.add_nodes_from(range(3));g.add_edge(0,1);return g
    if n % 2 == 0:
        return nx.complete_bipartite_graph(p-1, q+1)
    g = nx.complete_bipartite_graph(p,q)
    g.remove_edge(0,p)
    return g

def least(g, a):
    A = nx.to_numpy_array(g)
    return float(np.linalg.eigvalsh(float(a)*np.diag(A.sum(1))+(1-float(a))*A)[0])

def main():
    atlas = nx.graph_atlas_g()
    print('2026-11-15 atlas exact runner-up certificate; NetworkX atlas n<=7; arithmetic=integers')
    for n in range(3,8):
        gs = [nx.convert_node_labels_to_integers(g) for g in atlas if len(g)==n]
        K = nx.complete_bipartite_graph(n//2,n-n//2)
        T = candidate(n)
        assert sum(nx.is_isomorphic(g,K) for g in gs)==1
        assert sum(nx.is_isomorphic(g,T) for g in gs)==1
        for a in ALPHAS:
            other = [g for g in gs if not nx.is_isomorphic(g,K) and not nx.is_isomorphic(g,T)]
            guess = min(least(g,a) for g in other)
            runner = least(T,a)
            assert runner < guess - 1e-7, (n,a,runner,guess)
            t = Fraction(floor(((runner+guess)/2)*1000000),1000000)
            assert det_bareiss(matrices(T,a,t))<0, (n,a,'runner not below cutoff')
            assert det_bareiss(matrices(K,a,t))<0, (n,a,'minimum not below cutoff')
            fail = sum(not positive_definite(matrices(g,a,t)) for g in other)
            assert fail == 0
            print(f'n={n} alpha={a} atlas={len(gs)} competitor=1 cut={t} other_PD={len(other)} failures={fail} signs=negative_determinants')
    print('PASS: 20 exact rational separating cuts; 0 uncertified competitors')

if __name__=='__main__': main()
