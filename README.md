# The minimum least eigenvalue of A_α(G) for 0 < α < 1/2

Athanor lab (Aleksei Kudriashov, ORCID 0009-0006-8885-5338). Article: [`paper/main.pdf`](paper/main.pdf) (LaTeX source
`paper/main.tex`, `paper/refs.bib`); arXiv package and metadata in [`arxiv/`](arxiv/); certificates for every result in
[`certificates/`](certificates/); literature check in [`PRIOR-ART.md`](PRIOR-ART.md) (two independent audits in `priority/`).

## Results

1. **Nikiforov's Problem 30, all n (proof).** For 0 < α < 1/2 the least eigenvalue of A_α(G) = αD(G) + (1−α)A(G) over all graphs
   of order n ≥ 2 is at least [αn − √(α²n² + 4pq(1−2α))]/2, p = ⌊n/2⌋, q = ⌈n/2⌉, with equality only for K_{p,q}
   (`certificates/problem30odd/`). For even n this follows from published results; for odd n it was not found published in two
   independent bounded audits and the author's search (15 Nov 2026); this is not a novelty certificate.
2. **Stability (proof).** A graph whose least eigenvalue is close to the minimum is close to K_{p,q} in edit distance, with an
   explicit bound (`certificates/p30stability/`).
3. **Runner-up (exact finite certificates, conjecture).** At α ∈ {1/10, 1/4, 2/5, 49/100} and n ≤ 7 the second smallest value is
   attained only by K_{1,3}, K_{2,3}−e, K_{2,4}, K_{3,4}−e; at small α (1/20 for n = 5, 7; 1/40 for n = 9) the odd-order
   runner-up is K_{p,q} plus one edge inside the larger part (exact integer certificates over all graphs,
   `certificates/guide-checks/runnerup-exact.log`). For odd n the runner-up switches at α_n (≈ 0.0731, 0.0601, 0.0492 for
   n = 5, 7, 9; n·α_n → about 1/2); an algebraic crossing equation from equitable quotients and a two-branch conjecture are in
   the article. The adjacency case α = 0 is credited to Yu–Fan (2012).

## Provenance

Found and checked in the audit of V. Nikiforov, *Merging the A- and Q-spectral theories*, Appl. Anal. Discrete Math. 11 (2017)
81–107 (arXiv:1607.03015): repository [athanor-aalpha-spectral](https://github.com/iwasborninbali/athanor-aalpha-spectral),
Zenodo [10.5281/zenodo.24239158](https://doi.org/10.5281/zenodo.24239158). ChatGPT research agents (OpenAI) wrote the
proofs, programs and the article draft in separate threads; Anthropic Claude, the lab's guide, read the proofs line by line,
re-derived the finite results with independent programs (nauty `geng`, exact integer arithmetic,
`certificates/guide-checks/`), found the small-α runner-up switch, and checked every reference. The author of record takes
responsibility for the content. Text and data CC BY 4.0, code MIT.

Archived on Zenodo: version 1.0, [doi:10.5281/zenodo.24241615](https://doi.org/10.5281/zenodo.24241615) (all versions: [doi:10.5281/zenodo.24241614](https://doi.org/10.5281/zenodo.24241614)).
