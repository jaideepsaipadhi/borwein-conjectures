# Cubic Borwein conjecture, residue class m ≡ 2 (mod 3): certificates

Computer-assisted proof accompanying `paper/ms.tex`. Every analytic inequality the proof uses is checked in ball
arithmetic (Arb via python-flint). This directory contains only the programs used by the current proof.

## Requirements
Python 3.12 and `pip install -r requirements.txt` (python-flint 0.9.0, mpmath 1.3.0, numpy 2.4.4, scipy 1.17.1, sympy 1.14.0).

## Reproduce everything
    ./run_all.sh              # ~3 h on one core
    JOBS=16 ./run_all.sh      # parallelises the uniform Laplace grid

## Certificates, logs and the line that must appear

| Paper | Program | Log | Required output |
|---|---|---|---|
| Cert. 3.1 finite range n ≤ 243 | exact_small_n.py 1 200; 201 225; 226 243 | logs/01_exact_small_n_{a,b,c}.log | `EXACT n in [..]: violations 0` in all three |
| Prop. 3.2 band | band_lemma.py 36000, band_tail_rigorous.py | logs/02_band_lemma.log, 02b_band_tail.log | `(a) ... True`, `(b) ... True`, `BAND TAIL: g(i) < 0 for all i >= 36000` |
| Thm 2.3 / App. A expansion | sz_expansion_audit.py, expansion_check.py 36000 | logs/03_sz_audit.log, 03b_expansion_check.log | `E_con=50.150108`; `max |rest|/bound = 0.8660` |
| Cert. 2.5 Poisson | poisson_check.py | logs/04a_poisson.log | `max |eps| ... <= 1.420` |
| App. B checks (not part of proof) | hand_checks.py, koksma_cauchy_check.py | logs/04b, 04c | Peano residual 0; ratios <= 1 |
| Prop. 2.6 end-to-end (not part of proof) | identity_check.py | logs/04d_identity.log | rel.err < 1e-8 on every line |
| Cert. 5.2 centre tables | centre_mid_uniform{,2,3}.py, centre_g3_uniform{,2,3}.py, centre_clow_uniform.py | certificates/*_table.txt | no Python assertion errors |
| Cert. 4.3 uniform Laplace grid | run_uniform.py (36 chunks) | logs/uni_*.log | every line `UNRESOLVED 0` |
| Prop. 4.4 closing lemma | ray_asymptotic2.py | logs/07_closing_lemma.log | `CLOSING LEMMA (v >= 3000, all N): CERTIFIED` |
| Prop. 6.1 coverage (+ centre err, check E) | coverage.py | logs/08_coverage.log | `COVERAGE COMPLETE` |

## Notes
- Ball arithmetic: a printed bound is an upper (or lower) end of a certified enclosure; assertions abort the run on failure.
- `expansion_check.py` uses mpmath floating point with adaptive precision; it is corroboration of Appendix A, not part of the proof.
  The shipped log (i <= 108000) is from the development run.
