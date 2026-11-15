# Guide's independent checks (Claude, 2026-11-15)

All graphs of a given order come from nauty `geng` 2.9.3, a generator independent of the programs in the other folders.

| File | What | Status |
|---|---|---|
| `p30_scan.py`, `p30-n4..9.log` | least eigenvalue minimum over all graphs of order 4-9 at 8 values of alpha: no violations, unique minimiser K_{p,q} | floating-point corroboration of Theorem 2.2 |
| `p30_second.py`, `p30-second.log` | runner-up at alpha in {0.1, 0.25, 0.4, 0.49}, orders 4-9: unique; K_{1,3}, K_{2,3}-e, K_{2,4}, K_{3,4}-e, K_{3,5}, K_{4,5}-e | floating-point corroboration |
| `p30_second_a0.py`, `p30-second-alpha0.log` | runner-up at alpha = 0 (adjacency), orders 4-9: for odd n it is K_{p,q} plus one edge inside the larger part, not K_{p,q}-e | floating-point corroboration |
| `p30_second_small.py`, `p30-second-small-alpha.log` | runner-up at alpha in {0.005, ..., 0.09}, orders 5-9: odd n switches from K_{p,q}+e to K_{p,q}-e | floating-point corroboration |
| `p30-threshold.log` | the switching value alpha_n for odd n = 5..25 (bisection between the two candidates): 0.0731, 0.0601, 0.0492, ...; n alpha_n increases towards about 1/2 | numerical |
| `runnerup_exact.py`, `runnerup-exact.log` | exact integer certificate (Bareiss, Sylvester) over **all** graphs of order 5, 7 (alpha = 1/20) and 9 (alpha = 1/40): exactly K_{p,q} and K_{p,q}+e (edge inside the larger part) have least eigenvalue below the rational cutoff t | **exact**; refutes the runner-up conjecture in the form "K_{p,q}-e for every alpha in (0,1/2)" for odd n |

Replay: `geng -q N | python3 SCRIPT N [args]` as in the first lines of each log; `runnerup_exact.py n a b u v` takes alpha = a/b and t = u/v.
