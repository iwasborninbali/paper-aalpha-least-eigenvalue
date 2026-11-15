"""Guide's exact certificate of the runner-up for Problem 30 at a rational alpha = a/b and order n.
Reads all graphs of order n (nauty geng, graph6) on stdin. Given a rational cutoff t = u/v, it certifies with integer
arithmetic (Bareiss elimination, Sylvester's criterion on M = b v (A_alpha - t I)) that exactly the graphs listed by the
caller have least eigenvalue below t (negative determinant somewhere along the leading minors), and every other graph has
A_alpha - t I positive definite. Usage: geng -q n | python3 runnerup_exact.py n a b u v"""
import sys
from fractions import Fraction
n, a, b, u, v = map(int, sys.argv[1:6])
def decode(line):
    data = [c - 63 for c in line]; m = data[0]; bits = []
    for x in data[1:]: bits.extend((x >> k) & 1 for k in range(5, -1, -1))
    adj = [[0] * m for _ in range(m)]; i = 0
    for j in range(1, m):
        for k in range(j):
            if bits[i]: adj[k][j] = adj[j][k] = 1
            i += 1
    return adj
def leading_minors_positive(M):
    """Bareiss fraction-free elimination; returns True iff every leading principal minor is positive."""
    M = [row[:] for row in M]; size = len(M); prev = 1
    for k in range(size):
        if M[k][k] <= 0: return False
        for i in range(k + 1, size):
            for j in range(k + 1, size):
                num = M[i][j] * M[k][k] - M[i][k] * M[k][j]
                assert num % prev == 0
                M[i][j] = num // prev
        prev = M[k][k]
    return True
below = []
total = 0
for raw in sys.stdin.buffer:
    line = raw.strip()
    if not line: continue
    total += 1
    A = decode(line); d = [sum(r) for r in A]
    M = [[(a * v * d[i] - b * u) if i == j else (b - a) * v * A[i][j] for j in range(n)] for i in range(n)]
    if not leading_minors_positive(M):
        below.append((line.decode(), sum(d) // 2, sorted(d)))
print(f"n={n} alpha={a}/{b} t={u}/{v} graphs={total} below_t={len(below)}")
for g in below: print("  ", g)
