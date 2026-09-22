# Borwein-type sign patterns: the cubic Borwein conjecture and the Third Borwein conjecture

Paper and computer-assisted certificates for two results, obtained by one method.

Write `P_{p,n}(q) = prod_{k <= pn, p ∤ k} (1 - q^k)`.

- **Theorem A (cubic Borwein conjecture, remaining class).** For all `n >= 1` and all `m ≡ 2 (mod 3)` with
  `m <= 9n^2/2`, the coefficient `[q^m] P_{3,n}(q)^3` is `0` if `m < 3n` and `< 0` if `m >= 3n`. Together with
  Wang–Krattenthaler (classes 0 and 1) this proves the cubic Borwein conjecture.
- **Theorem B (Third Borwein conjecture).** For all `n >= 1` and all `m`, `[q^m] P_{5,n}(q)` is `>= 0` if `5 | m` and
  `<= 0` otherwise.

The method: factor `P_{p,n}^δ = G_p^δ E_n` (`G_p = (q;q)_∞/(q^p;q^p)_∞`, `E_n` a tail with non-negative
coefficients), use an explicit circle-method expansion of `G_p^δ` with a proved error constant, and cover the
coefficient range by an exact range, a band, a Laplace (saddle-point) region and a centre. Every analytic inequality is
verified in Arb ball arithmetic, uniformly in `n`. The cubic theorem is **not** used in the proof of Theorem B; only the
method is shared.

## Layout

```
paper/          main.tex, refs.bib, sections/*.tex (one file per section), main.pdf
cubic/          certificates, logs and tables for Theorem A (copy of the cubic-borwein-conjecture package)
  certificates/   programs (Arb via python-flint)
  logs/           shipped logs (the certificate is the required line in each log)
  tables/         centre tables
  run_all.sh      regenerates every cubic log
  paper_original/ the earlier stand-alone manuscript ms.tex / ms.pdf
  cubic-borwein-class2/  older package, for reference only (not used)
third/          certificates and logs for Theorem B
  certificates/   the programs cited in Section 5 of the paper, the modules they import
                  (offpeak_common.py, cls012_lib.py, lap34_lib.py) and the exact table g5_coeffs.pkl
  logs/           the shipped logs cited in the paper
  verification/   referee and audit reports (REFEREE_*.md, AUDIT*.md, GAP_*.md, ARB_F_RUN.md) and the source
                  write-ups THIRD_BORWEIN_PROOF.md, SEC3_PROOF.md
  run_all.sh      regenerates every Theorem-B log
Makefile        targets pdf, verify-cubic, verify-third, smoke-third
requirements.txt
LICENSE         (placeholder; no licence chosen yet)
```

## Building the PDF

Requires a TeX distribution with `pdflatex`, `bibtex` and `latexmk`.

    make pdf          # = cd paper && latexmk -pdf main.tex

## Rerunning the certificates

Python 3.12 and `pip install -r requirements.txt` (python-flint 0.9.0, mpmath 1.3.0, numpy 2.4.4, scipy 1.17.1,
sympy 1.14.0). All scripts must be run from their `certificates/` directory (they read `g5_coeffs.pkl`,
`cls012_lib.py`, `lap34_lib.py` by relative path); the `run_all.sh` scripts and the Makefile do this.

    make smoke-third     # < 2 min: saddle endpoint, saddle coverage, off-peak piece (N), Laplace tail cells, band table
    make verify-cubic    # all Theorem-A logs, ~3 h on one core (JOBS=16 parallelises the Laplace grid)
    make verify-third    # all Theorem-B logs; see below

`third/run_all.sh` accepts `SKIP_EXACT=1` (skip the exact range `n <= 979`) and `SKIP_F=1` (skip the Arb run of the
off-peak piece (F)), and `JOBS=k` for piece (F). Piece (F) reruns write `rerun_arb_F_*.log` so that the shipped
`arb_F_*.log` files (which were produced in several part-runs) are not overwritten.

A single certificate can be rerun by hand, e.g.

    cd third/certificates && python3 offpeak_tb.py 821 0.25     # must print "N0=821: piece (N) CERTIFIED for all n >= N0"

The table of certificates, commands, logs and required log lines is Section 5 of the paper (and
`cubic/README.md` for Theorem A).

## Runtimes and memory

Measured on a 2-core machine (from the logs or the audits).

| step | time | memory |
|---|---|---|
| cubic: exact range `n <= 243` | ~7 min | < 1 GB |
| cubic: band lemma (exact to i = 36000) | ~5 min | small |
| cubic: uniform Laplace grid (36 chunks) | the bulk of the 3 h; seconds to minutes per chunk, parallel over chunks | small |
| cubic: everything (`run_all.sh`) | ~3 h single core | |
| third: exact range `n <= 920` | ~80 min total | grows with n |
| third: exact segment 920–1000 (valid up to n = 979) | ~63 min to n = 979 | ~6 GB needed for NB = 1000 (killed after n = 979 on the original machine) |
| third: band exact run `291 <= n <= 2000` | 50–192 s | moderate |
| third: Laplace classes 0–2 / 3–4 | 8 s / 76 s | small |
| third: near peak + combination (two halves) | ~11 min per earlier single run | small |
| third: off-peak (E), (E0), (N) | 57 s, 6 s, 2 s | small |
| third: off-peak (F) in Arb | ~1.5 h on 2 cores | small |

Smoke test performed from the new location (`third/certificates/`): `centre3_saddle_end.py`, `cls012_cov.py`,
`offpeak_tb.py 821 0.25`, `cls012_lap.py 2.0 12.0 0.02`, `cls012_lap.py tail 11.99`, `lap34_cert.py tail 11.99`,
`offpeak_em.py 821 6.0 32`, `offpeak_em0.py 821 20 32`, `bandfix_coeffs.py`, `bandfix_analytic.py 10001`,
`cls012_band.py` all reproduce the required lines of the shipped logs.

## Trust base

1. The written arguments of the paper (Sections 2–4) and the standard facts of Appendix B.
2. Correctness of Arb ball arithmetic and FLINT exact integer arithmetic (through python-flint), and of the Python
   interpreter running the certificate programs as written.
3. IEEE-754 floating point is used only in the superseded float-interval version of off-peak piece (F); mpmath only in
   consistency checks and the independent recomputation of constants.

## Licence

TODO: no licence has been chosen. See `LICENSE`.

## Validation

Every certificate of both proofs was re-run in full from a clean copy and matched its shipped log; see
`validation/VALIDATION.md`. (`validation/run_remaining.sh` re-runs the two longest jobs, ~1 h, ~8 GB RAM.)

## Before submission (open items)

- **Wang–Krattenthaler scope (resolved by symmetry; one verbatim check left).** Berkovich–Dhar describe their cubic
  result as "all of κ_n, one half of δ_n, one half of θ_n". Since P³ is palindromic of degree 9n² ≡ 0 (mod 3), the
  mirror swaps classes 1 and 2, so this is class 0 plus one mirror pair of halves; the open half is the thin one
  (class 2, m ≤ 9n²/2, where Jacobi's identity cancels the main term), which WK §11 says resists their method and which
  Theorem A proves. Still worth reading WK Theorem 10.x verbatim before submission.
- **Related unrefereed work.** Two anonymous GitHub repositories (user `aasojsdoagaosgj`, AI-drafted, unreviewed)
  claim partial results: p = 5 powers 2 and 3 for large n, and the cubic class m ≡ 2 (mod 3) for n ≥ 80000.
  Neither covers the Third conjecture (power 1) or all n. Decide whether to mention them.
- Remaining bibliography TODOs (edition/theorem numbers, a journal version of Wang–Krattenthaler, arXiv:2609.22324
  author list): `grep -rn TODO paper/`.
- Licence and archival DOI for the certificate repository.
- Library versions: the cubic logs were produced with mpmath 1.3.0, the p = 5 logs with mpmath 1.4.1; both
  re-validate under the versions in `requirements.txt`.
