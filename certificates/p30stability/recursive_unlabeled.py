#!/usr/bin/env python3
"""Independent recursive unlabeled generator (no NetworkX atlas), numerical rank check.

Every n-vertex graph has a vertex whose removal gives an (n-1)-vertex graph;
all vertex-neighborhood subsets are tried, then duplicate isomorphism classes
are removed with NetworkX VF2. Independent catalogue construction from atlas.
"""
import networkx as nx
import numpy as np
from fractions import Fraction

ALPHAS = (0.1, 0.25, 0.4, 0.49)

def signature(g):
    return (g.number_of_edges(), tuple(sorted(d for _,d in g.degree())),
            sum(1 for _ in nx.enumerate_all_cliques(g) if len(_)==3))

def generate():
    levels={1:[nx.empty_graph(1)]}
    print('2026-11-15 recursive non-atlas isomorphism-class generator, VF2 + LAPACK')
    for n in range(2,8):
        buckets={}
        for h in levels[n-1]:
            for mask in range(1<<(n-1)):
                g=h.copy();g.add_node(n-1)
                g.add_edges_from((i,n-1) for i in range(n-1) if mask>>i&1)
                key=signature(g)
                bucket=buckets.setdefault(key,[])
                if not any(nx.is_isomorphic(g,old) for old in bucket):
                    bucket.append(g)
        levels[n]=[g for bucket in buckets.values() for g in bucket]
        assert len(levels[n])=={2:2,3:4,4:11,5:34,6:156,7:1044}[n], (n,len(levels[n]))
        print(f'n={n} generated_unlabeled={len(levels[n])} expected_class_count=PASS')
    return levels

def runner(n):
    p,q=n//2,n-n//2
    if n==3:
        g=nx.empty_graph(3);g.add_edge(0,1);return g
    if n%2==0:
        return nx.complete_bipartite_graph(p-1,q+1)
    g=nx.complete_bipartite_graph(p,q);g.remove_edge(0,p)
    return g

def main():
    levels=generate()
    for n in range(3,8):
        K=nx.complete_bipartite_graph(n//2,n-n//2)
        T=runner(n)
        for alpha in ALPHAS:
            vals=[]
            for g in levels[n]:
                A=nx.to_numpy_array(g)
                val=float(np.linalg.eigvalsh((1-alpha)*A+alpha*np.diag(A.sum(axis=1)))[0])
                vals.append((val,g))
            k=[v for v,g in vals if nx.is_isomorphic(g,K)]
            t=[v for v,g in vals if nx.is_isomorphic(g,T)]
            other=[v for v,g in vals if not nx.is_isomorphic(g,K) and not nx.is_isomorphic(g,T)]
            assert len(k)==len(t)==1
            assert k[0]<t[0]-1e-8 and t[0]<min(other)-1e-8
            print(f'n={n} alpha={alpha:.2f} runner_unique=PASS next_gap={min(other)-t[0]:.8f}')
    print('PASS: independent generation through n=7; 20 numerical second-place checks')
if __name__=='__main__':main()
