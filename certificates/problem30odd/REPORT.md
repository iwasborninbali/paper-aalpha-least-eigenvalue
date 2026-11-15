# Problem 30, odd order: exact value and uniqueness (2026-11-15)

Athanor lab (Aleksei Kudriashov). Base repository: `iwasborninbali/athanor-aalpha-spectral` (`main`); this contribution resides entirely in `problem30odd/`.

## Question

For finite simple undirected graphs `G` of order `n`, define
`A_alpha(G) = alpha D(G) + (1-alpha) A(G)` and write `lmin(G,alpha)` for its least eigenvalue. Nikiforov's **Problem 30** [N17, Section 6.1, p.102; arXiv v1 p.21] asks for

`min_{|V(G)|=n} lmin(G,alpha)` for `0 < alpha < 1/2`.

The prior [P] priority audit explains the even-order result from published adjacency extremality and the known complete-bipartite spectrum. The focus here is **odd n**. Unlike a finite-grid conjecture, the following argument works for **all n >= 2 and all real alpha in (0,1/2)**.

## Results with evidence status

### Theorem (full written proof; odd and even n)

For all `n >= 2` and `0 < alpha < 1/2`, put `p=floor(n/2)`, `q=ceil(n/2)`. Then

```
F_n(alpha) = [alpha*n - sqrt(alpha^2*n^2 + 4*p*q*(1-2*alpha))]/2
min_{|V(G)|=n} lmin(G,alpha) = F_n(alpha).
```

**The unique minimizer, up to isomorphism, is `K_{p,q}`.** For odd `n`, `4pq=n^2-1`. This is a complete mathematical solution of the optimization question, *conditional only on the standard Rayleigh minimum principle and the elementary nonnegative-matrix Perron theorem*, which are explicitly invoked below. **Novelty is not established** by our bounded priority search.

#### Proof

Write `beta=1-alpha > alpha`.

**1. A minimizer can be taken bipartite.** Let `G` be any n-vertex graph and let `x != 0` be a real least-eigenvalue eigenvector of `A_alpha(G)`. Recall the edge-by-edge quadratic form

```
x^T A_alpha(G) x = sum_{uv in E(G)}
                       [alpha*(x_u^2+x_v^2)+2*beta*x_u*x_v].        (1)
```

Delete every edge `uv` for which `x_u*x_v >= 0`, obtaining `H`. For every deleted edge its contribution to (1) is nonnegative, because `alpha>0` and `beta>0`. Thus

```
lmin(H,alpha) <= (x^T A_alpha(H) x)/(x^T x)
              <= (x^T A_alpha(G) x)/(x^T x) = lmin(G,alpha).       (2)
```

All remaining edges join a strictly positive coordinate of `x` to a strictly negative coordinate. Isolated zero-coordinate vertices may be assigned to either side, so `H` is bipartite on the original n vertices. In particular, a global minimizer exists among bipartite graphs.

**2. For fixed bipartition sizes, the complete bipartite graph is uniquely best.** Let `H` be any bipartite spanning graph with sides `P,Q` of sizes `r,s >= 1` and `r+s=n`. Put `K=K_{r,s}` on those same vertex sets. Conjugation by the diagonal sign matrix `S` (`+1` on `P`, `-1` on `Q`) gives

```
S A_alpha(H) S = alpha*D(H) - beta*A(H) = -B(H),
B(H) := beta*A(H) - alpha*D(H),
lmin(H,alpha) = -lambda_max(B(H)).                              (3)
```

The constant-on-parts quotient of `B(K)` is

```
           [ -alpha*s   beta*s ]
           [  beta*r  -alpha*r ].
```

Its positive largest eigenvalue is

```
mu = [sqrt(alpha^2*n^2+4*r*s*(1-2*alpha))-alpha*n]/2 > 0.        (4)
```

`B(K)` has a strictly positive eigenvector `z`, with value `u>0` on `P` and `v>0` on `Q`, because the quotient has positive off-diagonal entries and `mu>0`. Its eigen-equations give

```
beta*v - alpha*u = mu*u/s > 0,
beta*u - alpha*v = mu*v/r > 0.                               (5)
```

For any subgraph `H` of `K`, coordinatewise

```
(B(H) z)_i = d_H(i)*(beta*z_opposite - alpha*z_i)
           <= d_K(i)*(beta*z_opposite - alpha*z_i) = mu*z_i.      (6)
```

To justify the spectral comparison even when `H` is disconnected, add the scalar matrix `alpha*(n-1)I` to `B(H)`; the result is symmetric and entrywise nonnegative. The Perron theorem supplies a nonnegative largest-eigenvalue eigenvector `y != 0` of `B(H)`. Multiplying (6) by `y^T` gives

```
[mu-lambda_max(B(H))]*(y^T z) >= 0.
```

Since `y^T z>0`, `lambda_max(B(H))<=mu`. The inequality is **strict if H is a proper subgraph of K**: take a component on which `y` is positive and which supports `lambda_max(B(H))`. Unless `H=K`, some vertex in that component lacks an edge to the opposite full side, making its inequality in (6) strict by (5). Indeed, if every vertex in a component had all its opposite-side edges, that component would contain all of `P` and `Q`, forcing `H=K`. Thus, by (3),

```
lmin(H,alpha) >= lmin(K_{r,s},alpha),
```

with equality exactly for `H=K_{r,s}`. If one side is empty, `H` is edgeless with least eigenvalue 0; it cannot minimize because (4) is positive for nonempty balanced parts.

**3. Balance is optimal and uniqueness lifts to all graphs.** For `K_{r,s}` with `r+s=n`, (3)-(4) give

```
lmin(K_{r,s},alpha) = [alpha*n-sqrt(alpha^2*n^2+4*r*s*(1-2*alpha))]/2.
```

Because `1-2*alpha>0`, this expression is *strictly decreasing* in `rs`. Among nonnegative integers `r+s=n`, `rs` is uniquely maximized (up to exchanging sides) by `r=p,s=q`. Combining this with (2) proves the claimed minimum.

For uniqueness, suppose `G` attains that minimum, and form `H` using a least eigenvector `x` as in step 1. Every inequality in (2) and the fixed-partition comparison must be equality. Hence `H` is precisely `K_{p,q}`. Its negative least eigenvalue is **simple** (the negative root of its two-cell quotient, with all remaining eigenvalues `alpha*p` or `alpha*q` nonnegative), and its least eigenvector has nonzero coordinates at every vertex, strictly opposite in sign on the two classes. Equality in (2) therefore forces the same `x` to be a least eigenvector of `H` and to have no zeros. Every edge deleted from `G` would then have both endpoints of the same nonzero sign and would contribute `alpha*(x_u^2+x_v^2)+2*beta*x_u*x_v > 0` to (1), contradicting equality in (2). No edge was deleted, so `G=H=K_{p,q}`. QED.

**Sideways corollary (proof above).** The fixed-partition inequality is strict for every proper bipartite subgraph of `K_{r,s}`, for `0<alpha<1/2`. More locally, the bipartition sizes obtained by the sign pattern of *any* least eigenvector give the corresponding `K_{r,s}` lower bound to `lmin(G,alpha)` (zero coordinates can be allocated arbitrarily). The underlying strict comparison applies to `B=beta*A-alpha*D`, not to the ordinary adjacency Perron root, so it supplies a reusable extremal tool.

### Exact rational checkpoints (independent checker plus finite sweeps)

For odd `n`, write `alpha=a/b` reduced. Define the integer radicand

```
R(n,a,b) = a^2*n^2 + (n^2-1)*b*(b-2*a).
F_n(a/b) = (a*n - sqrt(R(n,a,b)))/(2*b).
```

The independent `certificate_checker.py` constructs exact algebraic identities and a certified **rational isolating interval** of width `1/(2*b*10^12)` using only integer `isqrt` (no floating point). It checks all `3 x 14 = 42` cases, with `n=5,7,9` and

`alpha in {1/1000, 1/100, 1/20, 1/10, 1/8, 1/6, 1/5, 1/4, 3/10, 1/3, 2/5, 9/20, 49/100, 499/1000}`.

At `alpha=1/4`, the **exact minima** and unique graph types are:

| n | exact minimum | unique minimizer | certificate status |
|--:|:--|:--|:--|
| 5 | `(5-sqrt(217))/8` | `K_{2,3}` | Full proof; 1024 labelled graphs checked with integer PD test |
| 7 | `(7-sqrt(433))/8` | `K_{3,4}` | Full proof; 1044 atlas isomorphism classes checked with integer PD test |
| 9 | `(9-sqrt(721))/8` | `K_{4,5}` | Full proof certifies all graphs; 8662 *sampled labelled masks* checked with integer PD test, **not a graph6 enumeration** |

For **every** one of the 14 rational parameters (not merely `1/4`), `exact_sweep.py` tests `A_alpha(G)-tI` positive definite for each non-extremizer on `n=5` and `n=7`, using exact integer Bareiss elimination/Sylvester's criterion; `t` is a rational **upper** bound on the radical minimum, hence this is stronger than the target inequality. `n5_exact.log` proves positive definiteness for 1014 of 1024 labelled graphs at each parameter, with 10 labelled `K_{2,3}` copies omitted because their exact least eigenvalue is the target. `n7_exact.log` proves positive definiteness for 1043 of 1044 atlas graphs at each parameter; the omitted graph is `K_{3,4}`. Graph-list coverage for `n=7` relies on NetworkX's published atlas generation (count 1044), whereas the analytic proof is independent of graph-list completeness. The program compares exact integers, not approximate eigenvalues.

`n9_probe.py` tests the 126 labelled copies of `K_{4,5}`, their one-edge edits, and 4000 deterministic random masks, giving 8662 distinct labelled graphs (8536 non-extremizers) for each alpha. These are genuine **exact** positive-definiteness checks, but explicitly **not** an exhaustive computational enumeration at `n=9`. The **mathematical exhaustive n=9 minimum** follows instead from the full real-alpha proof, covering all `2^36` labelled graphs without iterating through them. The optional `--graph6` mode of `exact_sweep.py` permits a separate replay of all 274668 nonisomorphic n=9 graphs against a free, externally supplied `graph9.g6` file, but **that optional scan has not been run**.

**Evidence status summary:** complete analytical proof (all n, all real `0<alpha<1/2`); exact integer spectral certificates (n=5 labelled, n=7 atlas); independent integer radical checker (n=5,7,9, 14 alpha values); n=9 sampled exact stress tests (not computational exhaustion). No floating-point eigenvalues and no paid services were used.

### Replay and integrity

Run from repository root with Python 3 and `networkx` for n=7 (no other third-party Python libraries are needed):

```sh
python problem30odd/certificate_checker.py > /tmp/certificate_checker.replay.log
python problem30odd/exact_sweep.py --n 5 > /tmp/n5_exact.replay.log
python problem30odd/exact_sweep.py --n 7 > /tmp/n7_exact.replay.log
python problem30odd/n9_probe.py > /tmp/n9_probe.replay.log
for name in certificate_checker n5_exact n7_exact n9_probe; do
    diff -u problem30odd/$name.log /tmp/$name.replay.log
done
(cd problem30odd && sha256sum -c SHA256SUMS)
```

`exact_sweep.py --n 9 --graph6 /path/to/graph9.g6` additionally requires the full nonisomorphic data file, e.g. the free Brendan McKay graph collection [G9] (274668 simple graphs of order 9). It is **optional**, not a prerequisite to the all-n theorem, and is not represented as run here. If supplied, the script checks the count and masks, and validates every non-`K_{4,5}` graph at all 14 parameters. `SHA256SUMS` includes every submitted script, log, and this report, but not itself.

## Priority check

**Verdict: the exact formula for `K_{r,s}` is ALREADY PUBLISHED; the universal minimization/unique odd-n conclusion was NOT FOUND in the bounded search, but we do not claim that it is new.** The theorem is presented as an independently proved deduction, not as an established first publication.

Sources examined and scope:

- [N17] V. Nikiforov, *Merging the A- and Q-spectral theories*, Applicable Analysis and Discrete Mathematics **11** (2017), 81–107, DOI [10.2298/AADM1701081N](https://doi.org/10.2298/AADM1701081N), [journal PDF](https://doiserbia.nb.rs/img/doi/1452-8630/2017/1452-86301701081N.pdf), [arXiv:1607.03015](https://arxiv.org/pdf/1607.03015). Problem 30 is p.102 in the journal (p.21 in preprint). Proposition 38 (journal p.105; preprint Proposition 37, p.23) already gives the `K_{r,s}` eigenvalue formula. Theorem 25 gives a related but different extremal claim: **largest** `A_alpha` eigenvalue among **bipartite** graphs for `alpha<1/2`; it does not itself state the least-eigenvalue global optimization.
- [P] Repo [`priority/REPORT.md`](https://github.com/iwasborninbali/athanor-aalpha-spectral/blob/main/priority/REPORT.md), sections 3 and "Search record": even-order result from published ingredients, explicit caution about unreviewed citation graph. Repo [`verify/REPORT.md`](https://github.com/iwasborninbali/athanor-aalpha-spectral/blob/main/verify/REPORT.md) tested small orders numerically; [`status/REPORT.md`](https://github.com/iwasborninbali/athanor-aalpha-spectral/blob/main/status/REPORT.md) correctly marked odd Problem 30 unsolved at that stage. [`attack/REPORT.md`](https://github.com/iwasborninbali/athanor-aalpha-spectral/blob/main/attack/REPORT.md) also contains no odd-order proof. These are prior lab work, not claimed as current discoveries.
- [LXS] H. Lin, J. Xue, J. Shu, *On the A_alpha-spectra of graphs*, Linear Algebra and its Applications **556** (2018), 210–219, DOI [10.1016/j.laa.2018.07.003](https://doi.org/10.1016/j.laa.2018.07.003), [arXiv:1709.00182, Section 3](https://arxiv.org/pdf/1709.00182). The arXiv text *explicitly recalls Problem 30 as the small-alpha question* (p.8 of PDF, lines around the start of Section 3) and then treats **alpha>1/2**; Theorem 3.3 goes in the opposite direction (`lambda_min(G) <= lambda_min(K_{p,q})`) for bipartite graphs when `alpha>1/2`. This is a relevant source, **not** publication of the present theorem.
- [LDSSh] S. Liu, K. C. Das, S. Sun, J. Shu, *On the least eigenvalue of A_alpha-matrix of graphs*, Linear Algebra and its Applications **586** (2020), 347–376, DOI [10.1016/j.laa.2019.10.025](https://doi.org/10.1016/j.laa.2019.10.025). Abstract explicitly concerns an upper bound for least eigenvalues when `1/2 <= alpha <= 1` and equality questions; the accessible abstract does **not** settle Problem 30 in our range. Full theorem-by-theorem text was **not** available for inspection here.
- [WWT] S. Wang, D. Wong, F. Tian, *Bounds for the largest and the smallest A_alpha eigenvalues of a graph in terms of vertex degrees*, Linear Algebra and its Applications **590** (2020), 210–223, DOI [10.1016/j.laa.2019.12.039](https://doi.org/10.1016/j.laa.2019.12.039). The accessible full-text/abstract gives Theorem 1.3: a degree-based lower bound for least eigenvalues of **bipartite** graphs, equality when connected regular/semiregular. It does not display the universal order-n odd minimum in the checked sections; it is closer prior art than a generic spectral paper.
- [LDS] S. Liu, K. C. Das, J. Shu, *On the eigenvalues of A_alpha-matrix of graphs*, Discrete Mathematics **343** (2020), 111917, DOI [10.1016/j.disc.2020.111917](https://doi.org/10.1016/j.disc.2020.111917). Lemma 2.5 records the sharp least-adjacency-eigenvalue bound and extremal `K_{p,q}`, as reported by [P].
- [G9] B. McKay, [Combinatorial Data: Graphs](https://users.cecs.anu.edu.au/~bdm/data/graphs.html) and free graph6 distribution. Source lists **34**, **1044**, **274668** nonisomorphic simple graphs on 5, 7, 9 vertices respectively; the present execution uses NetworkX atlas for n=7, not McKay's n=9 file.

Additional public queries on 2026-11-15: `"Nikiforov Problem 30 minimum least eigenvalue"`, `"smallest eigenvalue" "A_alpha" "complete bipartite" minimum odd`, `"least A_alpha eigenvalue" "balanced"`, `"On the least eigenvalue of Aα-matrix of graphs"`, and searches targeting theorem names/DOIs, `Turán` and `K_{floor(n/2),ceil(n/2)}`. Relevant 2021 and 2024–2026 works on least eigenvalue or A_alpha-index were checked by accessible summaries but no direct duplicate of the precise theorem was verified. **This is not an exhaustive cited-by harvest**; see [P] for earlier trace and the same limitation. No priority novelty claim should be made until the entire citation graph and potentially inaccessible full texts are audited.

## Not covered

1. **No actual computational enumeration of all 274668 unlabeled n=9 graphs or 2^36 labelled n=9 graphs was run.** The **general proof does mathematically cover them**, so n=9's minimum is exact; the `n9_probe.log` is correctly described as a stress test. Replaying the optional `--graph6` scan would add independent finite corroboration.
2. The `n=7` set of isomorphism representatives is taken from the installed NetworkX atlas, not generated by an independently implemented isomorphism-class enumerator; finite coverage depends on that known data source.
3. Irrational alpha is covered by the written proof, **not** by finite program sweeps or a formal Lean/Coq mechanization. The scripts certify 14 rational alphas and exact algebraic roots; none of their finite output alone proves the continuous theorem.
4. The literature review is bounded, not a complete paper-by-paper review of every citation of Nikiforov (2017) or of all 2026 material. The full text of [LDSSh] was not accessible. Novelty/original priority is therefore **not certified**.
5. No explicit lower bound on the **gap to the runner-up graph**, no exact second-minimizer classification, and no generalization to unequal edge weights or other matrix pencils is proved.

## Next steps

1. **Priority first:** search the complete citation list and inspect all accessible work on least `A_alpha` eigenvalues (`alpha<1/2`) for an earlier statement of the global extremum and its equality case. If prior, cite it prominently rather than calling this discovery new.
2. **Deepen sideways:** determine the *second-smallest* possible least eigenvalue and a quantitative stability theorem (gap in terms of deleted cross edges, imbalance, or edges inside eigenvector-sign classes). The strict comparison in (5)-(6) identifies where a sharp gap might come from.
3. Independently run the optional full n=9 graph6 PD scan at the same 14 rational parameters, verifying that all 274667 nonextremal classes pass; record the input graph6 hash and output log. Use only the freely available data; do not represent this as completed until actually executed.
4. Independently formalize the two decisive steps: sign-edge deletion (Rayleigh) and strict subgraph Perron domination. The full written proof above is complete, but a separate formal checker would strengthen publication readiness.
