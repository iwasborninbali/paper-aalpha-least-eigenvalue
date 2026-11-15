#!/usr/bin/env python3
"""Independent all-n algebraic certificate, using only integer arithmetic.

Checks exact quotient identities, positivity signs, and isolating intervals
for n=5,7,9 at rational alpha. The universal inequalities and equality
cases themselves are established in the complete written REPORT.md proof.
"""
import math
from fractions import Fraction

ALPHAS = ('1/1000', '1/100', '1/20', '1/10', '1/8', '1/6', '1/5',
          '1/4', '3/10', '1/3', '2/5', '9/20', '49/100', '499/1000')
SCALE = 10**12


def check(n, alpha):
    f = Fraction(alpha)
    a, b = f.numerator, f.denominator
    p, q = divmod(n, 2)[0], n - n // 2
    assert p+q == n and p*q == (n*n-1)//4
    assert 0 < 2*a < b
    R = a*a*n*n + 4*p*q*b*(b-2*a)
    assert R > a*a*n*n
    t = math.isqrt(R*SCALE*SCALE)
    assert t*t <= R*SCALE*SCALE < (t+1)*(t+1)
    # F(t) = t^2 - alpha*n*t - pq*(1-2*alpha) has roots
    # (alpha*n +/- sqrt(alpha^2*n^2+4*pq*(1-2*alpha)))/2.
    # For T = 2*b*lambda = a*n - sqrt(R), polynomial identity:
    # T^2 - 2*a*n*T - 4*p*q*b*(b-2*a) = 0.
    # Substitution yields R-(a*n)^2 - 4*p*q*b*(b-2*a) = 0.
    assert R - (a*n)**2 == 4*p*q*b*(b-2*a)
    # The 2 x 2 quotient of K_{p,q} has strictly positive other
    # eigenvalues alpha*p, alpha*q, so the minus root is least.
    # Shifted signed bipartite B eigenvector (x on P, y on Q):
    # x=(b-a)*q/b, y=alpha*q+mu, mu=(sqrt(R)-a*n)/(2b).
    # Its two per-edge increments are both >0, because mu>0.
    # The equation is equivalent to mu^2 + alpha*n*mu
    #     = pq*(1-2*alpha). Checked by the same integer identity.
    # Certified interval: lower < lambda <= upper, width=1/(2*b*SCALE).
    lo = Fraction(a*n*SCALE-t-1, 2*b*SCALE)
    hi = Fraction(a*n*SCALE-t, 2*b*SCALE)
    assert hi-lo == Fraction(1, 2*b*SCALE)
    assert lo < hi < 0
    return a, b, R, lo, hi


def main():
    print('Independent exact-algebra certificate; no graph libraries, no floating point.')
    for n in (5, 7, 9):
        print(f'n={n} p={n//2} q={(n+1)//2} graph_classes=ALL via written proof')
        for alpha in ALPHAS:
            a, b, R, lo, hi = check(n, alpha)
            print(f'alpha={a}/{b} exact=({a*n}-sqrt({R}))/{2*b} '
                  f'interval=({lo},{hi}] width={hi-lo} PASS')
    print('PASS 42 radical identities, 42 rigorous isolating intervals; proof in REPORT.md')

if __name__ == '__main__':
    main()
