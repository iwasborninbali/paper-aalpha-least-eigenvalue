# The minimum least eigenvalue of A_α(G) for 0 < α < 1/2

Athanor lab (Aleksei Kudriashov, ORCID 0009-0006-8885-5338). Article in preparation; the LaTeX source and PDF will be in
`paper/`. Every result below has its certificate in `certificates/`.

## Results

1. **Nikiforov's Problem 30, all n (proof).** For 0 < α < 1/2 the least eigenvalue of A_α(G) = αD(G) + (1−α)A(G) over all graphs of
   order n ≥ 2 is at least [αn − √(α²n² + 4pq(1−2α))]/2, p = ⌊n/2⌋, q = ⌈n/2⌉, with equality only for K_{p,q}
   (`certificates/problem30odd/`). For even n this follows from published results; for odd n we did not find it published.
2. **Stability (proof).** A graph whose least eigenvalue is close to the minimum is close to K_{p,q} in edit distance, with an
   explicit bound (`certificates/p30stability/`).
3. **Runner-up (corroborated for n ≤ 9, conjectured for all n).** The second smallest least eigenvalue is attained only by
   K_{n/2−1,n/2+1} (n even) and K_{p,q} minus one edge (n odd) (`certificates/p30stability/`, `certificates/guide-checks/`).

## Provenance

Results were found and checked in the audit of V. Nikiforov, *Merging the A- and Q-spectral theories*, Appl. Anal. Discrete
Math. 11 (2017) 81–107 (arXiv:1607.03015): repository
[athanor-aalpha-spectral](https://github.com/iwasborninbali/athanor-aalpha-spectral), Zenodo
[10.5281/zenodo.24239158](https://doi.org/10.5281/zenodo.24239158). ChatGPT research agents (OpenAI) wrote the proofs and
programs in separate threads; Anthropic Claude, the lab's guide, read the proofs line by line and re-derived the finite
results with independent programs (nauty `geng` and exact arithmetic, `certificates/guide-checks/`). The author of record
takes responsibility for the content. Text and data CC BY 4.0, code MIT.
