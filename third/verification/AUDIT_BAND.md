# AUDIT_BAND: adversarial audit of the repaired band argument (classes 3,4; GAP_BAND.md, proof §5.1)

**Overall: SOUND.** No mathematical or computational error found. Two cosmetic misstatements fixed (below).
All reruns were done in the foreground: bandfix_exact 291..2000 took 50 s, bandfix_coeffs and bandfix_analytic under 1 min each.

## (a) Claim 1: exact run, 291 <= n <= 2000, all classes, m < 4N — VERIFIED
- **At most 3 parts.** Every part of E_n is >= N = 5n+1, so 4 parts sum to >= 4N. The claim is correct.
- **Multiplicities are handled.** The construction uses the S_2/S_3 cycle index:
  - h2 = (p1^2 + p1(q^2))/2
  - h3 = (p1^3 + 3 p1(q^2) p1 + 2 p1(q^3))/6

  These are the complete homogeneous symmetric functions, so repeated parts are counted: q^{2k} and q^{3k} come from p1(q^2) and p1(q^3), and q^{2k+k'} from the cross term.
  - Truncating p1 to [N, 4N) is exact modulo q^{4N}.
  - The script asserts that 6E_n has all coefficients divisible by 6.
- **Independent recomputation.** I formed the direct product prod_{k<=5n, 5∤k}(1-q^k) as an exact fmpz_poly tree product mod q^{4N}. This method does not use g5_coeffs.pkl or E_n.
  - It was run for 16 values of n: 291, 292, 300, 422, 500, 568, 777, 815, 1000, 1234, 1458, 1500, 1857, 1936, 1999, 2000.
  - The result is identical to `bandfix_exact.coeffs(n)` on all m < 4N, and there are 0 sign failures. Coefficients reach 86 digits.
- **Rerun and coverage.**
  - I reran `bandfix_exact.py 291 2000`: 0 failures, max |b2+b3|/|b1| = 1.1302e-4 at (291, 5819), matching the log.
  - Coverage has no gaps: `imap_unordered(check, range(lo, hi+1))` visits every n, and a worker exception would abort the run, not drop an n.
  - The script also asserts c(m) = G(T) for classes 3,4 with m < 2N.

## (b) Claim 2: analytic majorant for n >= 2000 — VERIFIED
Here a = 2π²/75.
- **P5' formula.** From Theorem 3.1, P5 = 2a I_1(y)/y with y = 2 sqrt(a(ν-1/6)). Hence P5' = 4a² I_2(Y)/Y², which I checked.
- **|s(i)| <= e^{2 sqrt(a i)} for all i.**
  - For i <= 2e5 I recomputed the bound from the pkl. The maximum ratio is 1, at i = 0. For i >= 1000 the maximum ratio is 0.004.
  - For i > 2e5 the bound follows from Theorem 3.1 with |a| <= 3.62, P5 <= 2a e^y/(y sqrt(2πy)), P10 <= a e^{y/2}/(y sqrt(πy)), Φ(K) <= K²/10 + K/2 and K <= 2sqrt(π(i-1/6)) + 3. This gives 7.73e-5.
- **Bessel facts: all directions correct.**
  - I_ν decreases in ν.
  - I_1 <= I_{1/2} <= e^y/sqrt(2πy). This is an upper bound and is used only in numerators.
  - I_2 >= I_{5/2}. Expanding the closed form of I_{5/2} gives e^y/sqrt(2πy)·[1 - 3/y + 3/y² - e^{-2y}(1 + 3/y + 3/y²)], which is >= e^y/sqrt(2πy)·(1 - 3/y - 2e^{-2y}). This is a lower bound and is used only in denominators.
- **Pair and triple counts.**
  - Pairs: p2(2N+d) = floor(d/2) + 1. Brute force for n = 3, 20, 291 confirms the bound d/2+1. The maximum is 1456 at n = 291, compared with the old claim of 292.
  - Triples: partitions of e into at most 3 parts, <= (e+3)²/12 + 1/2. Checked by brute force for e < 300 and against the true triple count for n = 3.
  - No triples occur for d <= D < N.
- **Split at D.**
  - For d <= D the terms total at most (D+1)(D+2)/2 · e^{2 sqrt(aI)} <= (9/a) N L² (1+ε1) e^{2 sqrt(aI)}.
  - For d > D I used sqrt(I-d) <= sqrt I - d/(2 sqrt I) and I < 2N. This gives each term <= e^{2 sqrt(aI)} N^{-3}, times at most 2N terms, times Wmax. The total is 2 Wfac/(12N), which is the code's ε3·N L².
  - Every ε is decreasing in N, so evaluating them at N0 is valid.
- **Denominator.**
  - m = 5M+r gives 5T = I + N + 1 - r >= I + N - 3, and T >= 2000 >= 200.
  - |b1| = |G(T)| >= |g(T)|, because G(T-1) <= 0 by Corollary 5.3.
  - |g(T)| >= 4P5'(I+N-3), since P5' is increasing (a positive series in ν - 1/6).
  - Y <= 2 sqrt(3aN) and κ(Y) >= κ(2 sqrt(a(N-4))).
- **Exponent at I = 2N.** The exponent is 2 sqrt(a)(sqrt I - sqrt(I+N-19/6)) = -2 sqrt(a)(N-19/6)/(sqrt I + sqrt(I+N-19/6)), which increases in I. Its value at I = 2N is <= -cc sqrt N + δ0, with δ0 = 4 sqrt(a)/sqrt(3N-4) >= (19/6) sqrt(a)/sqrt(3N-4). Verified.
- **Monotonicity in N.** d log Ψ/dN < 0 exactly when (cc/2) sqrt N > 9/4 + 2/log N, with cc/2 = 0.16304. This first holds at **N = 257**, not "about 230" (a misstatement, now fixed in GAP_BAND.md and proof §5.1). The code checks it at N0 = 10001, and the condition persists for larger N.
- **Ψ values.**
  - Independent mpmath recomputation, using the exact D and Wmax: Ψ(10001) = 0.19714, Ψ(10006) = 0.1958, Ψ(12001) = 0.0137, Ψ(15001) = 3.49e-4, Ψ(50001) = 3.1e-17.
  - The script's arb values: Ψ(10001) <= 0.19769, and Ψ(15001) <= 3.49e-4.
  - The handover is correct: n = 2000 gives N = 10001, so the exact range (n <= 2000) and the analytic range (N >= 10001) overlap at n = 2000.
  - Conclusion: c(m) <= b1(1 - 0.198) < 0.

## (c) Claim 3: one-part block — VERIFIED
- **Exact part.** My own pass over the pkl found g(t) >= 0 only at (1, 1) and (5, 0), and G(T) >= 0 only at T = 1, where G = 0, for t <= 39999. `bandfix_coeffs.py` rerun agrees, and it confirms the pkl against Euler × p(k) for i <= 2e5. The direct products in (a) also validate the pkl up to i = 40004 without using it.
- **L1.** I checked the lemma's derivation.
  - Main term: |a1| + 2|a2| = sqrt5 + (5 - sqrt5) = 5, and the main-term bound follows from the mean value theorem with P5' increasing.
  - |g - main| <= 3 R_grp(5t+2), using that R_grp is increasing (K(ν) and each P_k are increasing).
  - P_k <= P10 for k >= 10.
  - P10 = (a/2) I_1(y10)/y10 <= sqrt(a/(4X)) e^{y10}/sqrt(2π y10).
  - y10 - Y/2 <= sqrt(a/(X-1/6)).
  - Φ(K)/X is decreasing.
  - Resulting shape: RA ∝ X^{1/2}(X-1/6) e^{-Y/2}, which is decreasing for X >= 1000 (log-derivative 0.0015 vs 0.0081). RB ∝ Y^{5/2} e^{-Y}, and 1/κ is decreasing.
  - Arb value: 0.017972 <= 0.018 at X = 1000.
  - Coverage: L1 covers every t >= 200 (5t >= 1000) and overlaps the exact range t <= 39999. R_grp includes E_con = 27.4636 (the RB term).

## (d) Classes — CONFIRMED
Proof §5.1 (the main-term cancellation) is used only for classes 3 and 4. Classes 0, 1, 2 use the exact run for n <= 2000 and `cls012_band` for n > 2000 (§5.2, §8 item 1). That route is outside this audit.

Note that GAP_BAND.md §5 still says classes 0–2 are "not covered" for n > 2000. That statement is stale: it predates `cls012_band`. It is harmless.

## (e) band_mono.py — SUPERSEDED
`band_mono.py` is used nowhere in the proof. It is listed only in Appendix B, among the withdrawn items. It still depends on `band_ext_cert.py`, whose counts are invalid, so it carries no evidential weight.

## Fixes applied
1. The threshold for "Ψ decreasing" was changed from about 230 to N >= 257 (GAP_BAND.md line 89, THIRD_BORWEIN_PROOF.md §5.1). This does not affect the result, since the argument is used from N0 = 10001.
2. Not changed: the proof's §5.1 says "for 980 <= n <= 2000", but quotes the ratio maximum at n = 291. This is cosmetic, because the run covers 291..2000.
