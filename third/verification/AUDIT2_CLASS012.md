# AUDIT2_CLASS012: second independent audit, classes 0, 1, 2 (band and Laplace)

*19 September 2026. Audit scripts are in the session scratchpad: `crt_engine.py`, `run_chunk.py`, `analyse.py` (exact c(m) by a new method), `amp.py` (an mpmath re-implementation of the saddle quantities that does not import `cls012_lib`), `lapX.py` + `mono.py` (the certificate evaluated at X > X_min). I wrote this review before reading GAP_AUDIT012.md or REFEREE_FULL.md; §8 compares the two.*

**Overall verdict: SOUND, for claims (1), (2) and (3) as now run.**
- I found no mathematical error, no term without a bound, and no coverage gap.
- Every number in the logs reproduces.
- 76 new exact points (69 inside the region), computed by a method independent of the first audit, all have the correct sign. Their error is 15–260 times smaller than the certified bound.
- What remains are two latent float issues in the code (neither is triggered) and some stale documentation (§7).

---

## 1. Re-derivation of the representation: agrees with the code

Derived from scratch:

- **Starting identity.** c(m) = Σ_j e(j)s(m−j) = e(m)·s(0) + Σ_{j<m} e(j)[a(m−j)P₅(m−j) + r(m−j)], with |r(i)| ≤ R(i) for i ≥ 1 (Theorem 3.1).
- **The amplitude.** a(i) = Σ_{h=1}^4 A_hζ^{hi}. I checked that A₁ζ^i + A₄ζ^{4i} = 2cos(2πi/5 + π/5), and that this agrees with w_h ζ^{−νh} in §3.1. There is no h = 0 term, and a(3) = a(4) = 0.
- **The Laplace step.** P₅(i) = Σ_{r≥1} a₅^r x^{r−1}/(r!(r−1)!) is L⁻¹[e^{a₅/u} − 1] at x = i − 1/6. Put r₁ = e^{a₅/u} − 1 − a₅/u. Its inverse is P₅ − a₅ for x > 0 and 0 for x < 0, so the terms with j ≥ m drop out. Then J_h = a₅H_h + K_h holds exactly, for any W > 0.
- **E_n on the contour.** log E_n(ζ^{−h}e^{−u}) = −Σ_c Σ_{t≥n} Φ_c(u(5t+c)). Lemma 6.1 with θ = c/5 gives −γ(σ)/u + Σ_c B₁(c/5)Φ_c(σ) + R_h. I re-derived the lemma's kernel, the value B₁(θ) = mean of k, and the constant 11/300 ≤ 1/12. The sign of the Ω term, ∫_σ^∞ Q = 5γ(σ), and the h-independence of Q all check.
- **Numerical check of the expansion.** Against a direct sum (n = 291 and 980; V = 2, 2.08, 5, 12; η = 0, 0.3, 1, 2; h = 1, 2), the maximum of |R|/R̄ is 0.092. There is no violation.
- **"Depends on V alone."** The exponent is a₅/u + (m−1/6)u − γ(σ)/u = X·Ξ_V(η), once (m−1/6)W² = α(V) is imposed. Ξ′_V(0) = 0 is exactly the definition α = a₅ − γ + Vγ′. n and m enter only through X = 1/W.
  - The remaining pieces are Ω_h (bounded; Re Ω_h(W) = 0 by the pairing c ↔ 5−c) and R_h (∝ W).
  - The normalisation (W/2π)·e^{XΞ₀}·√(2π/(XΞ₂)) equals Amp.
- **Where a(m) comes from.** Factoring out Amp gives the main term M = Σ_h A_hζ^{hm}e^{Ω_h(W)}. This is a(cls) with each root weighted by the unimodular factor e^{Ω_h}. It does not cancel, because |Ω| ≤ 0.8·Li₁(e^{−2}) = 0.117.
- **Code match.** `Mval`, `Xi_series`/`Xi_val`, `crude`, `window` and the mid loop implement exactly this derivation.

## 2. Term inventory: no holes

| term of the exact representation | bound | code |
|---|---|---|
| e(m) | e^{mW}E_n(e^{−W}) | `Te` |
| a₅A_hζ^{hm}H_h, 4 values of h | a₅e^{(m−1)W}E_n(e^{−W}) | `4*a5*Fm` |
| ρ (Theorem 3.1 remainder), Bessel part k = 10 and k ≥ 15, and E_con | Rankin at W; φ(5j) ≤ 4j; K ≤ √(4π(a₅X²+1)) + 2 | `Trho` |
| K_h: the −(1 + a₅/u) part of r₁ on \|η\| ≤ 2 | 2η₁(W + a₅)·Fm/2π | `Tcr`, first term |
| K_h, \|y\| ∈ [2W, 1] (full r₁) | e^{a₅X/5} + 1 + a₅X/2 | `Tcr`, second term |
| K_h, \|y\| > 1 | a₅²e^{a₅} | `Tcr`, third term |
| K_h, e^{a₅/u} part, window | the ψ bound; Riemann sum at right endpoints | `eps_win` |
| Gaussian mass outside the window | erfc(cc/√2) | `gtail` |
| K_h, e^{a₅/u} part, η₀ < \|η\| ≤ 2 | d_k-tiles × e^{2ω̄ + R̄} | `eps_mid` |
| E_n on the contour | EM remainder R̄ | `Rbar` |
| E_n(e^{−W}) in all crude terms | exp(0.8X·Li₂ + 4Li₁), from monotone blocks of 5 | `Fm` |

- I checked the constants myself: φ(5j) ≤ 4j; P_k ≤ W(e^{a_kX} − 1)e^{xW}; a₁₀ = a₅/4; a_k ≤ a₅/9 for k ≥ 15; and K(i) ≤ K for i ≤ m.
- The mid grid covers [0, 2.0] exactly, since the last float is exactly 2.0.
- There are **170** η-tiles, not the "169" written in §6.3 (documentation only).

## 3. Monotonicity in X: holds

`lapX.py` replaces X_min by a given X. I evaluated each term at X/X_min ∈ {1, 1+10⁻⁷, 10⁻⁵, 10⁻⁴, 10⁻³, 0.01, 0.03, 0.1, 0.5, 1, 4, 9, 99, 999} (i.e. X/X_min − 1 for the entries after the first) in five tiles: [2, 2.02], [4, 4.02], [7, 7.02], [9.78, 9.80] (the worst), and the snapped last tile [11.979999999999832, 12]. I did the same for the tail cell (V_T = 11.99).

- ε_win, ε_mid and E₂ are all non-increasing, and every margin is non-decreasing.
- Worst tile, class 2: margin 0.8339 at X_min, 0.8381 at 1.01X_min, 1.0458 at 2X_min.
- Tail cell, class 2: margin 0.8589 at X_min, rising to 1.3695 at 1000X_min.

My first run did flag a non-monotonicity. That was an artefact of my own test, which used the X_min rounded to 0.1 in the `info` output, and so evaluated slightly below X_min. With the exact X_min the flag disappears.

## 4. Independent exact validation (new method)

**Method.** c(m) is computed as S_n = G₅·E_n modulo 26–82 primes below 2³¹, followed by symmetric CRT.

- G₅ = pentagonal series × `inverse_series_trunc` of the pentagonal series in q⁵.
- E_n = [Σ_r q^{r(5n+1)}/(q;q)_r] · [Σ_r (−1)^r q^{5(rn + r(r+1)/2)}/(q⁵;q⁵)_r]. The first factor counts partitions into parts > 5n; the second is Euler's identity for ∏_{k>n}(1 − q^{5k}).
- Neither factor uses fmpz truncated products or the pickle.

**Validation of the method.**

- It matches brute-force products for every m at n = 3, 7, 12.
- G₅ mod p matches `g5_coeffs.pkl` for all i ≤ 2·10⁵.
- E_n has e(j) = 0 on (0, N) and e(j) = [5∤j] on [N, 2N).
- The number of primes is set from a rigorous-style bound log|c| ≤ −m log r + Σ log(1 + r^k), plus 70 bits. The result is unchanged when 4 primes are dropped.

**Points: 76 (m, n) in total (7 of them, at n = 291, have t < 4 and lie outside the region), all three classes each time, all with the correct sign.**

| n | V (points) | \|c/Amp − M\| | certified 4δ + E₂ at the actual X |
|---|---|---|---|
| 980 | 2.0814 (μ = 0.10500, m = 1008415–7), 2.09, 5, 9.78, 11.98, 12.00 (either side), 12.02, 15, 17.95 (t = 4.00) | 0.0003–0.0095 | 0.149–0.393 (0.296–0.548 at X_min) |
| 291 | 2.0815 (μ = 0.10500), 2.10, 5, 9.78 (t = 4.00), 11.99 and 12.01 (the last two have t = 2.66, which is outside the region) | 0.0009–0.0215 | 0.294–0.548 |
| 1500 | 9.78, 11.99, 12.00, 12.01, 15, 22.2 (t = 4.00) | 0.0013–0.0077 | 0.19–0.34 |
| 2500 | 11.99, 12.01, 20, 28.7 (t = 4.00) | 0.0010–0.0059 | 0.16–0.29 |

- The largest error is 0.0215, at n = 291 with V ≈ 12.
- The error scales like 1/X. At n ≥ 980 the error is 25–260 times below the bound.
- At V = 12.0001 (class 0, in the tail cell) the bound is 0.304. At V = 11.9998 (classes 1 and 2, in the last tile) it is 0.277. Both sides of 12 were exercised.

## 5. Band, claim (2): correct

- **Weights.**
  - For j < 2N, e(j) ≤ 1.
  - For j < 4N, e(j) ≤ 1 + (N+1) + ((N+3)²/12 + 1/2) ≤ N² ≤ m².
  - Checked exactly against the true e(j) for n = 5, 20, 57, 200, 400. For example, at n = 400 the largest e(j) with j < 4N is 175336, against the bound 336338.
  - GAP_CLASS012.md §A still writes the pair count as "N/2+1". That is false for j ∈ [3N, 4N); the correct count is N + 1. The proof text and the code comment have the right count, and N² absorbs either (documentation only).
- **N_min device.** Valid: m < 4N ⟺ N ≥ ⌊m/4⌋ + 1, and A is non-decreasing. It is uniform over every n > 2000.
- **Part A.** Exact integer comparisons (`abs(s[m])>corr`), with a worst ratio of 5.058·10⁻¹⁰. The only float in Part A is the check |s(i)| ≤ e^{c√i}. That check is harmless: the largest value of |s(i)|e^{−c√i} for i ≥ 1 is e^{−1.026} = 0.36, and i = 0 is exact.
- **Part B.**
  - I re-derived h₀: I₁ ≥ I_{3/2} ≥ e^Y(1 − 1/Y)/√(2πY), together with the √m − √(m−1/6) correction.
  - The P₁₀ factor uses √(πY/2) where √(πY) would do, which is conservative.
  - The φ-sum ≤ K²/10 + K/2 and the four exponents k₁–k₄ check.
  - The powers p = 1.25, 3.25, 1.75 and 0.75 satisfy √m₀ > 2p/k.
  - Q ≤ 5.139·10⁻²⁰ reproduces.
  - Wording: Q bounds (corrections + R)/(|a|_min·P₅), not corrections/|s(m)|.

## 6. Soundness of numeric comparisons (float, .mid, max/min, loop drift)

**V-tiling.**

- The float loop starts at exactly 2.0, the tiles are contiguous (`v=v2`), and the last tile is [11.979999999999832, 12.0].
- For all 500 tiles, the ball `arb((Va+Vb)/2,(Vb-Va)/2)` contains [Va, Vb], because the radius is rounded up. I checked this for every tile.
- The tail cell is run from `arb(11.99)` = 11.99 + 2·10⁻¹⁶, which lies inside the last tile.
- Coverage no longer depends on the snap.

**Coverage lemma (`cls012_cov.py`).** The loop variable is an arb (`v=v+w`), so there is no float drift. The tile bound α(v)/(v+w)² is taken over the balls, the start ball reaches below 0.001, and the loop ends only once v.lower ≥ 2. The comparison `>0.042` uses the double 0.042000000000000003, which errs in the safe direction.

**Band thresholds.** m₀ = 200000, N_min and 10001 are exact integers. n ≤ 2000 is covered exactly and n ≥ 2001 analytically.

**`.mid()` and `float(`.** `.mid()` is used only in the c₁ sanity assertion and in diagnostics. `amax`/`min` with a float key only select between valid bounds. After the fix to the pass test, a NaN fails safely: NaN > 0 is False, and `assert dk.lower()>0` raises.

**Latent issue L1 (new; not triggered).** `X4sup(b)` uses `ceil(b/0.01 - 1e-9)`. If η₀.upper fell in (0.01k, 0.01k + 10⁻¹¹], the sup|Ξ⁗| sub-balls would stop just short of η₀. This is the same kind of float endpoint as the V = 12 bug.

- I checked all 500 tiles and every mid-grid b ≤ 1: the sub-ball union always covers.
- The width-0.01 sub-balls have no gaps between them.
- Suggested fix: drop the `-1e-9`, or take one extra ball.

**Latent issue L2 (flagged by the first audit, still unfixed).** `X4.append(float(up(...).upper()))` rounds an upper bound to the nearest double, which can lower it by one ulp. This is immaterial. Suggested fix: keep it as an `arb`.

**Coverage of the V range.** The certified inequality needs V > 2 and V ∈ [2, 12] ∪ [11.99, ∞). Both hold.

## 7. Documentation defects (none affects a number)

- GAP_CLASS012.md is stale in four places:
  - it states the tail as "V ≥ 12 … 3.096/1.714/0.860" and gives the command `tail 12`;
  - it gives the crude rate as "κ − a₅/5" (the code now uses a₅/4);
  - it gives the pair count as "N/2+1";
  - it still points to `cls012_lap.log`, which holds the pre-snap run.
- THIRD_BORWEIN_PROOF.md §6.3 says "169 η-tiles"; the code has 170.

## 8. Comparison with the first audit (GAP_AUDIT012.md) and the referee

- **Agreement.** Both audits find the representation, the EM constant, the tail-cell Cauchy constants, the band and the coverage correct, and both reproduce every logged number. My exact points (76, via mod-p/CRT with the E_n Euler-identity construction) are disjoint in method from the first audit's 45 (fmpz truncated products), and the error sizes agree: 0.0003–0.0215 here against 0.0005–0.0185 there.
- **Their fixes are in place.**
  - The NaN-safe pass test is now `not all(>0)`.
  - The crude rate is now κ − a₅/4.
- **Their item still open.** The `float(upper)` for sup|Ξ⁗| is still in the code (L2).
- **What they missed.** The first audit (§7) stated "the tiles are contiguous on [2, 12] and the tail cell covers V ≥ 12". That was false as run: this is the referee's M1 sliver. I confirm M1 is now closed twice over, by the snap and by the tail from 11.99.
- **New in this audit.**
  - The latent `-1e-9` endpoint in `X4sup` (L1).
  - The false "N/2+1" pair count in GAP_CLASS012.
  - The stale tail and κ text in GAP_CLASS012.
  - The 169 vs 170 tile count.
  - A direct test that monotonicity in X holds at X > X_min.

## Verdicts

| claim | verdict | numbers |
|---|---|---|
| (1) Laplace, classes 0–2 | **SOUND** | 500 tiles, 0 failures; margins 3.0700 / 1.6879 / 0.8339; tail from 11.99 at 3.0950 / 1.7130 / 0.8589; monotone in X; 76 exact points correct, error ≤ 0.0215 against bound ≥ 0.149 |
| (2) Band, n > 2000 | **SOUND** | 5.058·10⁻¹⁰ (exact, m ≤ 2·10⁵); Q ≤ 5.139·10⁻²⁰; e(j) bound verified exactly |
| (3) Coverage V > 2 | **SOUND** | A ≥ 0.044242 on (0, 2], compared against 0.042; the exact roots at μ = 0.105 are V = 2.0814 (n = 291, 980) |
