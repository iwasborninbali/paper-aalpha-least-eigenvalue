#!/usr/bin/env python3
"""Independent, exact Fraction check of the ground-state identity.

Uses rational Perron eigenpairs of K_(r,s), including odd unequal parts;
all 2^6 graphs n=4 and 2^10 graphs n=5, and deterministic samples n=6,7.
No NumPy, NetworkX, or floating-point calculations.
"""
from fractions import Fraction as F
from itertools import combinations

# (r,s, alpha, mu, z_on_left, z_on_right) exact quotient eigenpairs
CASES = [(2,2,F(1,4),F(1),1,1),
         (2,3,F(1,11),F(2),6,5),
         (3,3,F(2,5),F(3,5),1,1),
         (3,4,F(1,15),F(3),8,7)]

def test(r,s,alpha,mu,zp,zq,mask):
    n=r+s;beta=1-alpha
    z=[F(zp) if i<r else F(zq) for i in range(n)]
    y=[F((i*7+mask%11)%13+1) for i in range(n)]
    pairs=list(combinations(range(n),2))
    edges=[e for k,e in enumerate(pairs) if mask>>k&1]
    cross=[e for e in edges if (e[0]<r)!=(e[1]<r)]
    internal=[e for e in edges if (e[0]<r)==(e[1]<r)]
    missing=[(i,j) for i in range(r) for j in range(r,n) if (i,j) not in cross]
    L=mu*sum(v*v for v in y)+sum(alpha*(y[i]**2+y[j]**2)-2*beta*y[i]*y[j] for i,j in cross)
    J=sum(alpha*(y[i]**2+y[j]**2)+2*beta*y[i]*y[j] for i,j in internal)
    E=sum(beta*z[i]*z[j]*(y[i]/z[i]-y[j]/z[j])**2 for i,j in cross)
    T=sum(mu*(y[i]**2/(s if i<r else r)+y[j]**2/(r if j>=r else s)) for i,j in missing)
    assert L+J == E+T+J, (r,s,mask,L,E,T)

def main():
    print('2026-11-15 rational ground-state identity checker (no third-party libraries)')
    for r,s,a,mu,zp,zq in CASES:
        assert (mu+a*s)*zp==(1-a)*s*zq
        assert (mu+a*r)*zq==(1-a)*r*zp
        n=r+s
        count=1<<(n*(n-1)//2) if n<=5 else 2500
        for i in range(count):
            mask=i if n<=5 else ((i*0x9E3779B97F4A7C15+0xA5A5A5A5) & ((1<<(n*(n-1)//2))-1))
            test(r,s,a,mu,zp,zq,mask)
        print(f'n={n} r={r} s={s} alpha={a} mu={mu} masks={count} rational_identity=PASS')
    print('PASS: exact identity checks with both equal and unequal bipartitions')
if __name__=='__main__':main()
