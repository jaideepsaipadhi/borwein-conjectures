# GAP_BAND: repair of the band region (NOTES section 5, PROOF section 4; audit items 19-21)

Region: n >= 291, 0 <= m < 4N, N = 5n+1 (n <= 290 is handled by the separate exhaustive computation).
Notation as in NOTES: s(i) = [q^i] G_5, a = a_5 = 2 pi^2/75, E_n = prod_{k>5n, 5 not| k} (1-q^k)^{-1},
c_n(m) = s(m) + b1 + b2 + b3, where b_j sums over partitions of m - i into j parts, each part >= N. Below 4N
there are at most 3 parts. For m = 3, 4 (mod 5), s(m) = 0 and b1 = G(T) = sum_{t<=T} g(t), with
g(t) = s(5t)+s(5t+1)+s(5t+2) and T = floor(m/5) - n.

**Verdict.** The band is now covered for every n >= 291 and every m < 4N in classes 3 and 4 (mod 5). The
proof has three parts: exact integer computation for 291 <= n <= 2000, an analytic majorant for n >= 2000
(N >= 10001), and an exact-plus-analytic one-part-block lemma. One input is cited rather than proved here:
the G_5 expansion bound with E_con = 27.4636, which is being checked separately. **Classes 0, 1, 2 are
covered only for n <= 2000** (see section 5).

## 1. Exact certification, 291 <= n <= 2000 (`bandfix_exact.py`, log `bandfix_exact.log`)

For m < 4N only partitions into at most 3 parts >= N contribute. Let p1 = sum_{N<=k<4N, 5 not| k} q^k. By the
cycle index of S_2 and S_3, working mod q^{4N}:

    6 E_n = 6 + 6 p1 + 3(p1^2 + p1(q^2)) + (p1^3 + 3 p1(q^2) p1 + 2 p1(q^3)).

Then c_n = G_5 E_n mod q^{4N} is computed exactly with FLINT `fmpz_poly.mul_low` (integers, no rounding).

* Validation: for n = 1..14 the result equals the directly expanded product S_n on all m < 4N.
* Input: `g5_coeffs.pkl` agrees with an independent recomputation, (Euler pentagonal series) times p(k) at
  q^5, for every i <= 200000 (`bandfix_coeffs.py`).
* **Result: zero sign failures in all five classes, for every n in [291, 2000] and every m < 4N.**
  Runtime is 3 minutes on 2 cores.
* For classes 3 and 4, 2N <= m < 4N, the largest exact ratio |b2+b3|/|b1| is **1.13e-4**, at n = 291,
  m = 5819 (t = 3.997). The old certificate's bound there was 0.4576.

This replaces the n = 291..~480 certification (audit 19). No monotonicity or sampling is used.

## 2. One-part block, and m < 2N (audit (d))

* **Exact** (`bandfix_coeffs.py`, t <= 39999, exact integers): g(t) < 0 for every t except t = 1 (g = 1) and
  t = 5 (g = 0). Every partial sum satisfies G(T) <= 0, with G(T) = 0 only at T = 1. The NOTES claim, stated
  for t <= 3999, is confirmed and extended tenfold. (NOTES says "g < 0 except t = 1, 5". At t = 5, g is 0,
  not positive, which is harmless.)
* **Lemma L1** (`bandfix_analytic.py`): 3 Rgrp(X+2) <= 0.018 P5'(X) for all X >= 1000. This gives
  g(t) <= -5P5'(5t) + 3Rgrp(5t+2) < -4P5'(5t) < 0 for all t >= 200, which overlaps the exact range.
  Proof sketch:
  * P5'(X) = 4a^2 I_2(Y)/Y^2 with Y = 2 sqrt(a(X-1/6)).
  * I_2 >= I_{5/2} (I_nu decreases in nu), and I_{5/2} has a closed form. Together these give
    I_2(Y) >= e^Y (1 - 3/Y - 2e^{-2Y})/sqrt(2 pi Y).
  * I_1 <= I_{1/2} gives I_1(y) <= e^y/sqrt(2 pi y). Also P_k <= P_10 for k >= 10, and
    sum_{5|k<=K} phi(k) <= K^2/10 + K/2.
  * The ratio is bounded by const * X^{3/2} e^{-Y/2} + const * Y^{5/2} e^{-Y}. Both terms decrease for
    X >= 1000, which is checked by the sign of the log-derivative, so the bound at X = 1000 holds for all
    larger X.
* Hence **G(T) <= 0 for every T**, and G(T) <= g(T) < 0 for T >= 200. So for classes 3 and 4 and m < 2N,
  c(m) = G(T) <= 0 for every n.
* `band_tail5.py` is now unnecessary. It is also sound only in spirit:
  * Its default T0 is 200000, not 200. Run at 200 it gives 0.0071.
  * It uses kappa2 measured at y0 as a bound for all y >= y0. That needs y^{1/2} e^{-y} I_2(y) to be
    increasing, which is true but not stated.
  * Its claim "y10 <= y0/2 + 0.001" is false at T0 = 200 (the true gap is 0.016). A slack factor of 2
    absorbs this.
  * The claim "ratio decreasing for y >= 9" is asserted, not checked.

  L1 replaces it.

## 3. Analytic majorant for N >= 10001 (n >= 2000), classes 3 and 4, 2N <= m < 4N (`bandfix_analytic.py`)

Write I = m - 2N, so 0 <= I < 2N. From the class of m, 5T >= I + N - 3.

* **Denominator.** G(T-1) <= 0 and g(T) < 0, so |b1| >= |g(T)|. By L1, |g(T)| >= 4P5'(5T) >= 4P5'(I+N-3),
  since P5' is increasing (P5 is a positive power series in x). The I_{5/2} bound above then gives an
  explicit lower bound.
* **Numerator.**

      |b2| + |b3| <= sum_{d=0}^{I} w(2N+d) |s(I-d)|,

  where w counts pairs and triples of parts. Two facts are used:
  * p2(2N+d) <= d/2 + 1, and p3(3N+e) <= (e+3)^2/12 + 1/2 (partitions of e into at most 3 parts).
  * |s(i)| <= alpha e^{2 sqrt(a i)} for all i >= 0, with **alpha = 1**. For i <= 200000 this is exact from
    the coefficients (the maximum, 1, is at i = 0). For larger i the cited bound gives at most 7.7e-5.

  Split the sum at D = floor(3 sqrt(2N/a) log N), which is less than N:
  * Terms with d <= D contain no triples and contribute at most (D+1)(D/2+1) alpha e^{2 sqrt(aI)}.
  * Terms with d > D contribute at most 2N Wmax alpha e^{2 sqrt(aI)} e^{-(D+1) sqrt(a/2N)}, which is at most
    about 1/6 times alpha e^{2 sqrt(aI)}.
* **Ratio.** The exponent 2 sqrt(a)(sqrt(I) - sqrt(I+N-19/6)) increases in I, so it is largest at I = 2N.
  Bounding each factor at its worst gives

      (|b2|+|b3|)/|b1| <= Psi(N) = C * N^{9/4} (log N)^2 exp(-2 sqrt(a)(sqrt3 - sqrt2) sqrt N).

  The constant C is evaluated at N0 from quantities that decrease in N. Psi is decreasing whenever
  0.163 sqrt N > 9/4 + 2/log N, which holds for N >= 257 (checked in the script at N0; earlier text said "about 230", corrected by AUDIT_BAND).
* **Psi(10001) <= 0.198** (arb). For comparison: Psi(15001) <= 3.5e-4, and Psi(9001) = 0.81 is also below 1.
  So c(m) < 0 for all n >= 2000 in this range, with a proven (not sampled) monotonicity. Overlap with the
  exact range is at n = 2000.

This fixes audit items (a) and (b):
* The old majorant's I_2 envelope 0.99 e^y/sqrt(2 pi y) is replaced by the rigorous I_{5/2} bound.
* Monotonicity in n and in t is now a derivative argument. The t-dependence is handled by bounding the
  worst I = 2N, not by sampling.

## 4. Other findings

* **`band_ext_cert.py` counts are wrong, not only loose.** It uses #pairs <= I/10 + 1 and
  #triples <= ((m-3N)/15+1)^2.
  * At n = 291, t = 4, the true maximum pair count is 1160 against the claimed 292.
  * The true triple count at m = 4N - 2 is 92223 against the claimed 9591.

  So the "0.4576" is not a certified bound. It survives in practice only because |s(i)| is concentrated where
  the weights are about 1 (the exact ratio is 1.1e-4). The certificate is superseded and no longer used.
* **`band_mono.py`**: fixed. Out-of-scope cells (t < 2) are now skipped explicitly. Any None, nan, inf or
  >= 1 result is counted as a failure. Rerun: 0 failures, worst 0.4576, but that uses the invalid counts
  above. The script is labelled as sampled and superseded.
* `band_majorant4.py`: superseded by section 3 (its crossing near n ~ 480 is moot).

## 5. Scope: which classes

NOTES section 5 treats **only classes 3 and 4**. The main-term cancellation |a(0)| = |a(1)| + |a(2)| is
specific to them. NOTES claims nothing about classes 0, 1, 2 in the band; section 1 says only that they "keep
an O(1) main term".
* Section 1 here certifies all five classes, but only for n <= 2000.
* **For n > 2000, classes 0, 1, 2 with m < 4N are not covered by the band argument.**
  * They must be covered by another region of the proof, or by an easy analogue: s(m) dominates with
    exponent gap 2 sqrt(a)(2 - sqrt3) sqrt N. That gap is smaller, so it needs a larger N0 or a longer exact
    range.
  * I have not checked whether PROOF.md's other regions reach m < 4N for classes 0, 1, 2. This is an open
    coverage question for the architecture, not for this region.

## 6. What remains non-rigorous

* The cited expansion bound (E_con = 27.4636), used for i > 200000 and in L1. Taken as given, per the brief.
* Standard Bessel facts: monotonicity of I_nu in nu (Cochran 1967; Watson), and the closed forms of I_{1/2}
  and I_{5/2}.
* Nothing else. Every finite range is exact integer arithmetic. Every infinite range is an arb evaluation at
  its endpoint plus an explicit monotonicity argument.
