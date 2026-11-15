# arXiv submission metadata (prepared by the guide; the owner submits)

- **Archive to upload:** `arxiv/paper-aalpha-least-eigenvalue.tar.gz` (`main.tex` + `main.bbl`; no figures). Checked: compiles
  with `pdflatex` twice from a clean folder (TeX Live 2025, pdfTeX 1.40.27), 14 pages, no warnings.
- **Title:** The minimum least eigenvalue of $A_\alpha(G)$ for $0<\alpha<1/2$
- **Authors:** Aleksei Kudriashov
- **Abstract** (1331 characters, limit 1920):

  For a simple graph $G$ of order $n$, let $A_\alpha(G)=\alpha D(G)+(1-\alpha)A(G)$, with $0<\alpha<1/2$. We determine the minimum of its least eigenvalue over all graphs of order $n\ge2$: writing $p=\lfloor n/2\rfloor$ and $q=\lceil n/2\rceil$, the minimum is $(\alpha n-\sqrt{\alpha^2n^2+4(1-2\alpha)pq})/2$, attained precisely by $K_{p,q}$ up to isomorphism. This addresses Problem~30 of Nikiforov (2017). For even $n$, the statement is a direct consequence of published adjacency extremality and the complete-bipartite spectrum; for odd $n$, that deduction is not sharp. We give a uniform proof by sign-pattern deletion and strict Perron comparison. A nonnegative ground-state deficit identity yields a quantitative edit-distance stability estimate. We give exact rational certificates for the second distinct minimum at four parameter values through order seven and at smaller parameters through order nine. The latter show that the odd-order runner-up changes with the parameter. Equitable quotients give an algebraic equation for the candidate crossing; we state a revised two-branch conjecture and distinguish exact certificates from numerical scans. We distinguish mathematical proofs, exact computations, and floating-point evidence, and make no claim of established publication priority beyond a bounded literature search.

- **Comments:** 14 pages, 3 tables. Programs, logs and SHA-256 certificates: https://github.com/iwasborninbali/paper-aalpha-least-eigenvalue ; Zenodo DOI 10.5281/zenodo.24241614.
- **Primary category:** math.CO (Combinatorics). **Cross-list:** math.SP (Spectral Theory) is optional.
- **MSC 2020:** 05C50 (primary), 15A18, 05C35.
- **ACM / report number / journal ref:** none.
- **Licence:** CC BY 4.0 (matches the repository).
- **AI assistance:** stated in the acknowledgement of the paper (ChatGPT research agents by OpenAI; Anthropic's Claude as
  orchestrator and verifier; the author takes responsibility). arXiv policy: AI tools cannot be listed as authors; the
  submitting author is responsible for the content.
- **Endorsement:** a first submission to math.CO from an account without earlier math submissions usually needs an
  endorsement by an established math.CO author; arXiv shows this after the category is chosen. Asking is the owner's decision.
- **ORCID:** 0009-0006-8885-5338 (link it in the arXiv account settings).
