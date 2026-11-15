# p30stability: runner-up and quantitative stability for the least A-alpha eigenvalue

Athanor lab (Aleksei Kudriashov). Investigation and literature search: **2026-11-15**. New work is confined to this directory. Base: [athanor-aalpha-spectral](https://github.com/iwasborninbali/athanor-aalpha-spectral), whose version 1.0 is archived at [Zenodo 10.5281/zenodo.24239159](https://doi.org/10.5281/zenodo.24239159). All programs use free tools.

## Question

For a simple graph G of order n, set A_alpha(G)=alpha D(G)+(1-alpha)A(G), with 0<alpha<1/2. Put beta=1-alpha, c=1-2alpha, p=floor(n/2), q=ceil(n/2), and

F_n(alpha)=[alpha*n-sqrt(alpha^2*n^2+4*c*p*q)]/2.

The earlier [problem30odd proof](../problem30odd/REPORT.md) and [guide's independent check](../guide/2026-11-15-independent-checks/README.md) show that F_n is the unique global minimum and K_{p,q} is the unique minimizer. The present questions are the **second distinct minimum and its graphs**, a quantitative edit-distance stability theorem, and fixed-size/fixed-bipartition analogues. Evidence below is labeled separately as proof, exact finite certificate, numerical crosscheck or unresolved.

## Results with evidence status

### R1: Quantitative stability -- PROVED for every n and every real 0<alpha<1/2

Let d_edit(G,K_{p,q}) be the least number of additions/deletions needed to make G isomorphic to the balanced complete bipartite graph. Define explicit positive constants

- mu_* = 2*c*(n-1)/n;
- C_* = 1+beta*sqrt(n-1)*(n-1)/mu_*;
- L_* = n^2*(n-1)^2*[1/alpha + (n-1)*(1+2*C_*)/mu_*].

**Theorem (full proof below).** Suppose lambda_min(A_alpha(G))<=F_n+epsilon with 0<=epsilon< -F_n. Then

**d_edit(G,K_{p,q}) <= n*sqrt(beta*n*epsilon/c) + L_* * epsilon.**  (S)

If epsilon>= -F_n, use the unconditional bound d_edit<=n(n-1)/2. In particular (S) forces G=K_{p,q} up to isomorphism when its right-hand side is less than 1. The constants are very crude; no claim of optimal scaling or novelty of ground-state transformations is made.

**Proof, part A: an exact, nonnegative deficit identity.** Fix a unit least-eigenvalue eigenvector x of G. Its eigenvalue is negative by the hypothesis. Thus x has both positive and negative coordinates: if x had just one sign, each edge contribution alpha*(x_i^2+x_j^2)+2*beta*x_i*x_j would be nonnegative. Partition the vertices into two nonempty classes P,Q by the signs of x, assigning any zero coordinates arbitrarily. Set r=|P|, s=|Q| and y_i=|x_i|, so ||y||=1. Define H to contain **only** those G-edges for which x_i*x_j<0. It is bipartite across P,Q. In particular edges incident to zeros are removed from H, regardless of where the zero vertices were assigned.

Let K=K_{r,s} on these two classes and B(J)=beta A(J)-alpha D(J) for a bipartite graph J across them. Let mu=lambda_max(B(K))=-lambda_min(A_alpha(K))>0 and let z>0 be its unit Perron eigenvector, constant within each class. Write w_i=y_i/z_i and d_i equal to the size of the opposite class in K. The complete-bipartite eigenvector equations say that for every opposite-class j,

beta*z_j-alpha*z_i=mu*z_i/d_i.  (A)

The edgewise Rayleigh contribution of H to x^T A_alpha(G)x is -y^T B(H)y. Every deleted edge of G has x_i*x_j>=0 and contributes alpha*(y_i^2+y_j^2)+2*beta*y_i*y_j>=0. Call the total of these deleted-edge contributions J.

Expanding both sides of the symmetric ground-state identity using (A) gives the exact formula

Delta := lambda_min(A_alpha(G)) - lambda_min(A_alpha(K))
= beta*SUM_{ij in E(H)} z_i*z_j*(w_i-w_j)^2
+ mu*SUM_i [miss_H(i)/d_i]*y_i^2
+ J.                                                             (1)

Here miss_H(i) is the number of K-cross-edges incident to i which are absent from H. Crucially, a cross-edge of G touching a zero x-coordinate may be counted in miss_H; this is harmless and necessary. All three terms are nonnegative. In particular Delta>=0 and Delta<=epsilon, because the complete-bipartite formula gives lambda_min(A_alpha(K))>=F_n.

For a separately replayable algebraic check of (1), identity_check.py uses exact Fraction arithmetic on rational Perron eigenpairs of balanced and unbalanced complete bipartite graphs. This check is supporting evidence, not a substitute for the displayed general proof.

**Proof, part B: missing and internal edges.** Write M=|E(K)\E(H)|, I=the number of G-edges with both endpoints in P or both endpoints in Q, h=<y,z>, and v=y-hz. Let z_min=min_i z_i. For r<=s, the ratio t between the two constant Perron coordinates satisfies

beta*r*t^2+alpha*(s-r)*t-beta*s=0.

The left-hand side is nonpositive at t=1 and nonnegative at t=sqrt(s/r). Hence 1<=t<=sqrt(s/r)<=sqrt(n-1). Therefore, with S=sqrt(n-1),

z_min^2>=1/[n(n-1)], and h>=z_min,

where h>=z_min follows from y>=0, ||y||_2=1 and ||y||_1>=1.

B(K) has Perron eigenvalue mu>0, while **all** its other eigenvalues are <=0: they are -alpha*r, -alpha*s (when the relevant multiplicities are nonzero), and the negative root of the two-part quotient. The complete-bipartite ground-state energy

E_K=beta*SUM_{ij in E(K)} z_i*z_j*(w_i-w_j)^2

is therefore at least mu*||v||^2 by orthogonal spectral decomposition. The portion of E_K along H is at most Delta by (1). For a missing edge ij, since y_i,y_j>=0 and z_i/z_j,z_j/z_i<=S,

beta*z_i*z_j*(w_i-w_j)^2
<= beta*S*(y_i^2+y_j^2).

Let Dmax=max(r,s)<=n-1. The middle term of (1) is at least (mu/Dmax)*SUM_i miss_H(i)*y_i^2. Consequently

E_K <= [1+beta*S*Dmax/mu]*Delta.

Also the exact positive-root formula

mu=2*c*r*s/[sqrt(alpha^2*n^2+4*c*r*s)+alpha*n]

and r*s>=n-1 show mu>=2*c*(n-1)/n=mu_* (the denominator is at most n). It follows that

||v||^2 <= C_* * Delta/mu_*.                            (2)

The elementary inequality (a+b)^2>=a^2/2-b^2, with a=h z_i and b=v_i, implies
y_i^2>=h^2*z_i^2/2-v_i^2. Since SUM_i miss_H(i)=2M and miss_H(i)<=Dmax,

SUM_i miss_H(i)*y_i^2>=h^2*M*z_min^2-Dmax*||v||^2.

Combining this with (1) and (2) gives

M <= z_min^(-4)*[Dmax*(1+C_*)/mu_*]*Delta.                (3)

Every within-class edge is among those deleted from H. Its nonnegative J-contribution is at least alpha*(y_i^2+y_j^2). With internal vertex degrees at most n-1, the same squared-coordinate inequality gives

I <= z_min^(-4)*[1/alpha+(n-1)*C_*/mu_*]*Delta.            (4)

Equations (3)-(4), z_min^(-4)<=n^2*(n-1)^2, and Dmax<=n-1 imply

M+I <= L_* * Delta <= L_* * epsilon.                       (5)

The *actual* G-cross-edges missing from K form a subset of E(K)\E(H), so M+I is an upper bound for the edits transforming G into K_{r,s} on its current sign classes. Note that none of this claims arbitrary bipartite edge addition is spectrally monotone; that shortcut would be invalid in general because A_alpha changes its degree diagonal too.

**Proof, part C: move to balanced sides.** Let f(r,s) be the known K_{r,s} least eigenvalue. Direct rationalization of the radical gives

f(r,s)-F_n
=2*c*(p*q-r*s)/[sqrt(alpha^2*n^2+4*c*p*q)+sqrt(alpha^2*n^2+4*c*r*s)]
>= c*(p*q-r*s)/(beta*n).

Let k=min(|r-p|,|r-q|). For even n, p*q-r*s=k^2; for odd n it equals k(k+1), so k^2<=p*q-r*s<=beta*n*epsilon/c. Reassigning k vertices between the two classes transforms K_{r,s} to a balanced K_{p,q} by at most n*k edge edits. Together with (5), this is (S). QED.

**Additional fixed-partition stability.** For bipartite H with *designated* sides r,s, the unique minimum is f(r,s) at K_{r,s}, as already proved in problem30odd. For its spectral deficit Delta=lambda_min(A_alpha(H))-f(r,s), (3) immediately gives the explicit edit bound

|E(K_{r,s})\E(H)| <= n^2*(n-1)^2 * [(n-1)*(1+C_*)/mu_*] * Delta.

This includes isolated vertices and disconnected graphs; a least eigenvector may have zeros and is handled in (1). The fixed-partition **minimum itself is earlier Athanor work**, not a new claim here.

### R2: Exact second minimum for n=3 -- PROVED for all 0<alpha<1/2

The four unlabeled graphs are K_{1,2}, K_2 union K_1, K_3, and the empty graph. The first is uniquely minimal by the baseline theorem. Their remaining least eigenvalues are, respectively, **2alpha-1**, 3alpha-1, and 0. Since 2alpha-1<3alpha-1 and 2alpha-1<0, the second value is exactly **2alpha-1** and its unique extremal graph is **K_2 union K_1**. This is a complete elementary proof, not a numerical observation.

### R3: Second-place finite classification for n=3..7 -- EXACT FOUR-POINT GRID ONLY

For alpha=1/10, 1/4, 2/5 and 49/100, the second distinct least eigenvalue and *unique isomorphism class* are as follows:

| n | second graph | value |
|---|---|---|
| 3 | K_2 union K_1 | 2alpha-1 |
| 4 | K_{1,3} | f(1,3) |
| 5 | K_{2,3} minus any single cross edge | lambda_min of this graph |
| 6 | K_{2,4} | f(2,4) |
| 7 | K_{3,4} minus any single cross edge | lambda_min of this graph |

Here f(r,s)=[alpha*n-sqrt(alpha^2*n^2+4*c*r*s)]/2. For the one-edge deletion, the value is unambiguously the least eigenvalue of its finite explicit adjacency/degree matrix; no unsupported simplified radical is claimed. The third suggested candidate, K_{p,q} plus an internal edge, does not win at these checkpoints.

**Exact method.** atlas_exact.py uses the standard NetworkX atlas of every isomorphism representative through order seven (4, 11, 34, 156, 1044 at n=3..7). At each of 20 (n,alpha) pairs a numerical exploratory search *only proposes* a rational threshold t between the candidate and its next-nearest competitor. An independent integer Bareiss determinant routine proves that A_alpha(K_{p,q})-tI and A_alpha(candidate)-tI each have negative determinant (hence a negative direction), whereas **every other atlas graph** has all positive leading principal minors of A_alpha(G)-tI (hence is positive definite by Sylvester). Exact separating t values and zero failures are in atlas_exact.log. All inequality certificates use integers; the exploratory floating-point values are not relied upon for certification. This proves the finite classification conditional on the standard atlas inventory.

**Second method.** labelled_numeric.py independently constructs **all** labeled graph masks for n=3..6 (8, 64, 1024, 32768 graphs), builds graph matrices without NetworkX, identifies the candidate's labeled orbit through permutations and evaluates least eigenvalues by NumPy eigvalsh. All 16 order/parameter cases agree (labelled_numeric.log). This genuinely independent graph generation and spectral method is **floating-point corroboration**, not a second exact proof. **Additional independently generated n=7 method.** recursive_unlabeled.py constructs all isomorphism classes from scratch, starting at one vertex, appending a new vertex with every possible neighborhood, and deduplicating each level by VF2 graph isomorphism. It does **not** read NetworkX's atlas. The generated class counts are 2,4,11,34,156,1044 at orders 2 through 7, and its NumPy/LAPACK second-place check agrees for **all 20** (n,alpha) pairs, including every n=7 class; see recursive_unlabeled.log. It is independently generated but numerical; the atlas scan remains the exact arithmetic method. Completeness follows inductively because every n-vertex graph has an (n-1)-vertex induced subgraph and each possible new-neighbor subset is visited.

**Unresolved beyond these cases.** Results suggest K_{n/2-1,n/2+1} at every even n>=4, and K_{(n-1)/2,(n+1)/2}-e at every odd n>=5. **This is a conjecture, not a proved all-n second-minimum result.** No claim about all real alpha is inferred from the four-point grid. The n=3 exception is proved separately.

### R4: Fixed edge count among connected graphs -- NOT SOLVED

Sign-deletion replaces G with a possibly disconnected bipartite H having fewer edges, so the baseline argument does **not** preserve either m or connectivity. When m>floor(n^2/4), there is no bipartite candidate with those parameters. For example the unique connected (n,m)=(3,3) graph is K_3, with least A_alpha eigenvalue 3alpha-1. For m below the bipartite edge ceiling, a general fixed-(n,m) minimization is still open here. The adjacency-alpha=0 fixed-size literature of Bell--Cvetkovic--Rowlinson--Simic (2008) is an essential precedent.

## Priority check (dated 2026-11-15)

**No novelty certification.** The order-n minimum, its uniqueness and fixed-bipartition minimizing formula are inherited, while ground-state or Doob representations are a familiar mathematical technique. An exact match for the particular quantitative bound (S) or a continuous-parameter runner-up theorem was *not found* in the bounded search; absence from results does not establish first publication.

Sources inspected and their specific scopes:

1. [Nikiforov 2017](https://doi.org/10.2298/AADM1701081N), *Merging the A- and Q-spectral theories*, Applicable Analysis and Discrete Mathematics **11** (2017), 81-107: Problem 30 and Proposition 38 already give the problem and K_{r,s} spectrum.
2. [Wang--Wong--Tian 2020](https://doi.org/10.1016/j.laa.2019.12.039), *Bounds for the largest and the smallest A_alpha eigenvalues of a graph in terms of vertex degrees*, Linear Algebra and its Applications **590** (2020), 210-223, Theorem 1.3: relevant bipartite least-eigenvalue degree inequality; the even-order result is obtainable from published ingredients, as the Athanor priority report already records.
3. [Liu--Das--Shu 2020](https://doi.org/10.1016/j.disc.2020.111917), *On the eigenvalues of A_alpha-matrix of graphs*, Discrete Mathematics **343** (2020), 111917, Lemma 2.5: already quotes the sharp adjacency extremum, associated with the earlier **Constantine 1985** theorem and equality case.
4. [Bell--Cvetkovic--Rowlinson--Simic 2008 I](https://doi.org/10.1016/j.laa.2008.02.032), Linear Algebra and its Applications **429** (2008), 234-241: abstract explicitly concerns minimum adjacency least eigenvalue at fixed connected order and size and sign-vector structure, including nonbipartite possibilities.
5. [Bell--Cvetkovic--Rowlinson--Simic 2008 II](https://doi.org/10.1016/j.laa.2008.06.018), Linear Algebra and its Applications **429** (2008), 2168-2179: fixed-size bipartite extremal structure and size transitions. [Stirling archive](https://www.stir.ac.uk/research/hub/publication/799844) also lists errata, which have **not** been exhaustively checked.
6. Prior internal reports, read **before** calculation: [verify](../verify/REPORT.md), [attack](../attack/REPORT.md), [status](../status/REPORT.md), [priority](../priority/REPORT.md), [classify22](../classify22/REPORT.md), [problem30odd](../problem30odd/REPORT.md), [oddlambda2](../oddlambda2/REPORT.md), [forkmatrices](../forkmatrices/REPORT.md) and the guide's [independent checks](../guide/2026-11-15-independent-checks/README.md). The prior problem30odd report identifies stability and second place as not covered.

**Public searches conducted on 2026-11-15:** "second least eigenvalue" graphs complete bipartite Constantine; "least eigenvalue second smallest among n vertex graphs"; "A_alpha least eigenvalue second minimum stability"; "Hong least adjacency eigenvalue second extremal bipartite Bell Cvetkovic Rowlinson Simic"; "second least adjacency eigenvalue" complete bipartite; "Constantine 1985 lower bound smallest eigenvalue graph second minimum". Searches included DOI pages, Elsevier abstracts, the Stirling archive and indexed introductions. The full Bell et al. papers/errata and every citing article were **not** checked theorem by theorem; this is not a complete cited-by audit. In particular, the guide's roughly 304-work OpenAlex review was directed at the *first* minimum, not second-place priority.

## Replay and SHA-256

From the repository root, using free Python 3 and free NetworkX/NumPy only for the first two programs:

    OPENBLAS_NUM_THREADS=1 python3 p30stability/atlas_exact.py > /tmp/p30_atlas.log
    OPENBLAS_NUM_THREADS=1 python3 p30stability/labelled_numeric.py > /tmp/p30_labelled.log
    python3 p30stability/identity_check.py > /tmp/p30_identity.log
    OPENBLAS_NUM_THREADS=1 python3 p30stability/recursive_unlabeled.py > /tmp/p30_recursive.log
    diff -u p30stability/atlas_exact.log /tmp/p30_atlas.log
    diff -u p30stability/labelled_numeric.log /tmp/p30_labelled.log
    diff -u p30stability/identity_check.log /tmp/p30_identity.log
    diff -u p30stability/recursive_unlabeled.log /tmp/p30_recursive.log
    (cd p30stability && sha256sum -c SHA256SUMS)

SHA256SUMS includes this report and every program/log but does not include itself. atlas_exact is a rational certificate only at the four rational alpha checkpoints and depends on the atlas inventory; labelled_numeric independently enumerates all labeled n<=6; recursive_unlabeled independently generates all unlabeled classes through n=7 and verifies the four-point ranking numerically. identity_check verifies (1) for rational Perron eigenpairs (balanced and unbalanced) with exact arithmetic, but does **not** mechanically certify the full theorem (S).

## Not covered

- **No all-n/all-alpha proof of the second distinct least eigenvalue or classification of ties.** Parity candidates are finite evidence only. No special transition in alpha has been excluded for larger n.
- No n>=8 runner-up enumeration and no exact irrational-alpha finite enumeration. At n=7, two independent catalogue constructions (atlas and recursive VF2) agree, but only the atlas spectral comparison uses exact arithmetic. The all-n result is the *stability bound*, not the runner-up conjecture.
- No asymptotically sharp constant in (S); L_* is huge and degenerates near alpha=0 and 1/2. No optimal epsilon scaling or inverse theorem.
- No fixed-(n,m) minimizer classification, even for all feasible connected bipartite edge counts.
- No comprehensive literature/citation priority certification, especially for adjacency runner-up results and the BCRS errata.

## Next steps

**Highest value:** close the all-n runner-up gap by proving strict comparison of (i) the nearly balanced complete bipartite candidates, (ii) balanced K_{p,q} minus one cross edge, (iii) balanced K_{p,q} plus one internal edge, **and all other proper bipartite subgraphs**. One must not assume that adding an edge always decreases lambda_min(A_alpha): adding it also increases the diagonal by alpha. The exact nonnegative identity (1) is the best starting point for extracting sharper one-edge gap estimates. Independently audit the complete BCRS adjacency papers and errata, and then tackle fixed-(n,m) optimization preserving connectivity.
