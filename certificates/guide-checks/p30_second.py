"""Guide's check of the runner-up for Problem 30: second smallest distinct least eigenvalue of A_alpha over all graphs of
order n (nauty geng), and which graphs attain it (numpy, tolerance 1e-9; corroboration, not proof)."""
import sys
import numpy as np
n = int(sys.argv[1]); alphas = [0.1, 0.25, 0.4, 0.49]
def decode(line):
    data = [c - 63 for c in line]; m = data[0]; bits = []
    for x in data[1:]: bits.extend((x >> k) & 1 for k in range(5, -1, -1))
    A = np.zeros((m, m)); i = 0
    for v in range(1, m):
        for u in range(v):
            if bits[i]: A[u, v] = A[v, u] = 1
            i += 1
    return A
lines = [l.strip() for l in sys.stdin.buffer if l.strip()]
A = np.array([decode(l) for l in lines]); D = np.zeros_like(A); idx = np.arange(n); D[:, idx, idx] = A.sum(axis=2)
deg = A.sum(axis=2); edges = (A.sum(axis=(1, 2)) / 2).astype(int)
for a in alphas:
    lmin = np.linalg.eigvalsh(a * D + (1 - a) * A)[:, 0]
    order = np.argsort(lmin); first = lmin[order[0]]
    rest = lmin[lmin > first + 1e-9]; second = rest.min()
    at = np.where(np.abs(lmin - second) < 1e-9)[0]
    desc = [f"{lines[i].decode()} e={edges[i]} deg={sorted(deg[i].astype(int).tolist())}" for i in at]
    print(f"n={n} alpha={a} second={second:.9f} count={len(at)} {desc}")
