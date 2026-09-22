# End-to-end validation, 19–22 September 2026

**Every certificate of both proofs was re-run in full from a clean copy; every result matches the shipped logs.**

Every certificate was re-run from a **clean copy** of the sources, and its output compared with the log the paper
cites. The logs of these re-runs are in `validation/third/` and `validation/cubic/`.

## Third Borwein conjecture

| Certificate | Re-run result | Matches shipped log |
|---|---|---|
| G_5 coefficients to 2·10⁵ (`g5_coef.py`) | rebuilt in 4 s | yes |
| Band, one-part block (`bandfix_coeffs.py`) | g(t) < 0 for t ≤ 39999 except t = 1, 5; G(T) ≤ 0 | yes |
| Band, exact n = 291–2000, m < 4N, all classes (`bandfix_exact.py`) | 0 sign failures; max ratio 1.1302e-4 | yes |
| Band, analytic n ≥ 2000 (`bandfix_analytic.py 10001`) | Ψ < 0.1977 for all N ≥ 10001 | yes |
| Band, classes 0–2, n > 2000 (`cls012_band.py`) | Q ≤ 5.1390e-20 | yes |
| Laplace, classes 0–2 (`cls012_lap.py`, tiles + tail from 11.99) | 500 tiles, 0 failures; margins 3.0700 / 1.6879 / 0.8339; tail 3.0950 / 1.7130 / 0.8589 | yes |
| Saddle coverage (`cls012_cov.py`) | root V ≈ 2.0816 | yes |
| Laplace, classes 3–4 (`lap34_cert.py`, tiles + tail from 11.99) | 500 tiles, 0 failures; 0.8363 / 0.8361; tail 0.8146 | yes |
| Centre saddle end (`centre3_saddle_end.py`) | μ(2.1) < 0.105: True | yes |
| Centre near-peak (`centre3_peak.py 821 210 64`) | worst err/margin 0.7977 / 0.8469 (classes 3 / 4) | yes |
| Centre combination (`centre3_comb.py`, two halves) | COMBINATION ASSERTED on all 210 boxes | yes |
| Off-peak (N) (`offpeak_tb.py`) | certified, worst −1.8 | yes |
| Off-peak (E0) (`offpeak_em0.py`) | certified, worst −0.23 | yes |
| Off-peak (E) (`offpeak_em.py`) | certified, worst −2.18 | yes |
| Off-peak (F), Arb (`offpeak_chunk_arb.py`, all 8 pieces) | CERTIFIED, all 163 rects OK, min margin 47.91 | yes |
| Exact, every coefficient, n = 1–820 (`centre3_exact.py`; 580–820 in `rerun_exact_580_820.log`, 1927 s) | 0 violations | yes |

The analytic certificates start at n = 821, so exact coverage is needed exactly up to n = 820, and all of it was re-run.

## Cubic Borwein conjecture (class m ≡ 2 mod 3)

| Step | Re-run result | Matches shipped log |
|---|---|---|
| 1 Exact n ≤ 243 (three segments) | 0 violations in each | yes |
| 2 Band lemma (`band_lemma.py 36000`) and band tail | g(i) < 0 for all i ≤ 36000; tail certified for i ≥ 36000 | yes |
| 3 Expansion audit (`sz_expansion_audit.py`) and corroboration (`expansion_check.py 36000`) | reproduced; max ratio 0.8660 | yes |
| 4 Hand lemmas (Poisson, hand checks, Koksma–Cauchy, identity) | reproduced (hand_checks needs `sympy`, now in requirements) | yes |
| 5 Centre tables (mid, g3, clow) | **byte-identical** to `tables/` | yes |
| 6 Uniform Laplace grid (36 chunks) | all 36: UNRESOLVED 0; worst B = 0.1169; every line identical to shipped logs | yes |
| 7 Closing lemma (`ray_asymptotic2.py`) | CERTIFIED | yes |
| 8 Coverage (`coverage.py`) | COVERAGE COMPLETE | yes |


## Discrepancy found and fixed during validation
Off-peak piece (E): the proof document quoted worst value −12.03; the current certificate gives −2.18. The −12.03
predates widening the lowest L-box to start at 0.0999 (done to cover the exact L = 0.1). The certificate only
requires the value to be negative; the documents now quote −2.18.
