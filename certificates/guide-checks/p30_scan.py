"""Guide's finite corroboration of the Problem 30 theorem: for every graph of order n (nauty geng, all graphs), the least
eigenvalue of A_alpha = alpha D + (1-alpha) A is at least F_n(alpha), with equality only for K_{floor(n/2),ceil(n/2)}.
Floating point (numpy eigvalsh) with a margin; it corroborates the written proof and is not itself a proof."""
import sys
import numpy as np
n = int(sys.argv[1]); alphas = [0.01, 0.1, 0.2, 0.25, 1/3, 0.4, 0.45, 0.49]
p, q = n // 2, n - n // 2
F = {a: (a * n - np.sqrt(a * a * n * n + 4 * p * q * (1 - 2 * a))) / 2 for a in alphas}
def decode(line):
    data = [c - 63 for c in line]; m = data[0]; bits = []
    for x in data[1:]: bits.extend((x >> k) & 1 for k in range(5, -1, -1))
    A = np.zeros((m, m)); i = 0
    for v in range(1, m):
        for u in range(v):
            if bits[i]: A[u, v] = A[v, u] = 1
            i += 1
    return A
mats = [decode(l.strip()) for l in sys.stdin.buffer if l.strip()]
A = np.array(mats); D = np.zeros_like(A); idx = np.arange(n); D[:, idx, idx] = A.sum(axis=2)
deg = A.sum(axis=2); edges = A.sum(axis=(1, 2)) / 2
kpq = (edges == p * q) & np.all(np.sort(deg, axis=1) == np.array(sorted([q] * p + [p] * q)), axis=1)  # K_{p,q} is the only graph with these degrees and pq edges? check below
for a in alphas:
    lmin = np.linalg.eigvalsh(a * D + (1 - a) * A)[:, 0]
    best = lmin.min(); arg = np.where(lmin < F[a] + 1e-9)[0]
    second = np.sort(lmin)[1]
    print(f"n={n} graphs={len(A)} alpha={a:.4f} F={F[a]:.9f} min={best:.9f} attaining={len(arg)} kpq_flag={kpq[arg].tolist()} gap_to_second={second - best:.6f} violations={(lmin < F[a] - 1e-9).sum()}")
