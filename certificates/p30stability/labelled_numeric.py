#!/usr/bin/env python3
"""Independent labelled-bitmask generator and LAPACK eigenvalue crosscheck."""
import itertools
from math import comb
import numpy as np

ALPHAS = [0.1,0.25,0.4,0.49]

def masks_for(g_edges,n):
    edges=tuple(itertools.combinations(range(n),2))
    lookup={e:i for i,e in enumerate(edges)}
    result=set()
    for perm in itertools.permutations(range(n)):
        mask=0
        for u,v in g_edges:
            a,b=sorted((perm[u],perm[v]))
            mask |= 1<<lookup[a,b]
        result.add(mask)
    return result

def target(n):
    p,q=n//2,n-n//2
    K=[(i,j) for i in range(p) for j in range(p,n)]
    if n==3: C=[(0,1)]
    elif n%2==0: C=[(i,j) for i in range(p-1) for j in range(p-1,n)]
    else: C=K[1:]
    return K,C

def main():
    print('2026-11-15 independent labelled exhaustive numeric crosscheck (n=3..6)')
    for n in range(3,7):
        edges=tuple(itertools.combinations(range(n),2))
        K,C=target(n)
        mink=masks_for(K,n)
        cand=masks_for(C,n)
        assert not mink & cand
        for alpha in ALPHAS:
            win=1e100; second=1e100; third=1e100; errors=0
            winners=runner_count=0
            for mask in range(1<<len(edges)):
                A=np.zeros((n,n),dtype=float)
                for i,(u,v) in enumerate(edges):
                    if mask>>i&1:A[u,v]=A[v,u]=1
                M=(1-alpha)*A+alpha*np.diag(A.sum(axis=1))
                value=float(np.linalg.eigvalsh(M)[0])
                if mask in mink:
                    winners+=1
                    win=min(win,value)
                elif mask in cand:
                    runner_count+=1
                    second=min(second,value)
                else:third=min(third,value)
            assert win < second-1e-7 and second < third-1e-7, (n,alpha,win,second,third)
            print(f'n={n} alpha={alpha:.2f} labelled={1<<len(edges)} K_copies={winners} runner_copies={runner_count} gap_to_third={third-second:.8f} ordering=PASS')
    print('PASS: 16 complete labelled-order/parameter checks; independent graph generation')
if __name__=='__main__':main()
