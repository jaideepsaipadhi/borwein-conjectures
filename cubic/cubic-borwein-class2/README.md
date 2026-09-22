# The remaining residue class of the cubic Borwein conjecture — code and certificates

Paper: `paper/cubic_borwein_class2.tex` (compiled PDF alongside).

**Claim.** For all n ≥ 1 and m ≡ 2 (mod 3), m ≤ 9n²/2: [q^m] P_n(q)³ = 0 if m < 3n and < 0 otherwise,
where P_n(q) = ∏_{k≤3n, 3∤k} (1−q^k). With Krattenthaler–Wang (arXiv:2201.12415) this proves the cubic Borwein conjecture.

**This is a computer-assisted proof.** Every certificate below is either exact integer arithmetic (GMP) or ball arithmetic (Arb via
python-flint) on an explicit covering of a compact parameter domain. Nothing uses floating point as a proof step.

## Reproduce
```
sudo apt-get install libgmp-dev          # or brew install gmp
pip install -r requirements.txt
./run_all.sh                             # ~2 h on one core;  JOBS=8 ./run_all.sh to parallelise the Laplace grid
```
Everything is regenerated from scratch into `logs/`; `logs/SUMMARY.txt` collects the verdicts. The expected summary is in
`logs_reference/SUMMARY.txt`.

## What each script certifies
| step | script | statement (paper) | expected |
|---|---|---|---|
| 1 | `exact/cubic_check.c` | Certificate 3.1: n ≤ 243 | 243 lines, all `violations=0 zeros_in_3n_to_D/2=0 nonzero_below_3n=0` |
| 2 | `certificates/band_certificate.py` | Certificate 4.2: g(i) < 0 for all i | `BAND LEMMA CERTIFICATE: PASS` |
| 3 | `certificates/centre_{mid,g3,clow}_uniform*.py` | centre tables, uniform in N ≥ 733 | 110 rows each, no assertion failures |
| 4 | `certificates/run_uniform.py` → `B_uniform.py`, `uniform_lib.py` | Certificate 5.1: Laplace grid v ∈ [2.95, 3000] | 36 lines `UNRESOLVED 0`, worst B ≈ 0.117 |
| 5 | `certificates/ray_asymptotic2.py` | Proposition 5.2: v ≥ 3000 | `CLOSING LEMMA ... CERTIFIED`, B ≤ 0.0082 |
| 6 | `certificates/coverage.py` | Certificate 6.1, Lemma 6.3, tiling of all (n, m) | checks A–G PASS, `COVERAGE COMPLETE` |
| 7 | `certificates/poisson_check.py` | Lemma 2.4 | max \|η\| ≤ 1.42 |

Independent consistency checks (not part of the proof): `identity_check.py`, `expansion_check.py`, `sz_expansion_audit.py`,
`hand_checks.py`, `koksma_cauchy_check.py`, `B_exact.py`/`run_exact.py`, `exact_accept.py`, `centre_accept.py`.

## Software used for the reference run
Python 3.12.3, python-flint 0.9.0, mpmath 1.3.0, numpy 2.4.4, scipy 1.17.1, sympy 1.14.0, GMP 6.3.0, gcc 13.

## Notes for readers checking the proof
- The effective Laplace parameter is μ̃ = (m − N − ¼)/N², not m/N² (paper, Remark after Prop. 2.5).
- Uniformity in N comes from treating ε = 1/N as an interval containing 0 plus asserted monotonicity (paper §5.2, §6.2).
- The sz_joint expansion (Appendix A) is proved in the paper for the case used (G₃, δ = 3).
