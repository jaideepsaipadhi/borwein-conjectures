# GAP_AUDIT012: adversarial audit of the classes-0,1,2 Laplace and band certificates

*18 September 2026. The audit scripts are in the session scratchpad: `aud_exact.py` (exact c(m) against Amp·M and the certified bound) and `aud_nan.py` (all tiles rerun with a finiteness check).*

**Overall verdict: SOUND.** I found no mathematical error and no error that changes a number. There are two small defects in the code or its documentation, and each needs only a one-line fix (items 2 and 4). This certificate does **not** share the defects found in its siblings:

- the main term is centred on the true saddle, and it includes the E_n factor;
- no NaN value is replaced by a finite number, and no NaN occurs in any tile;
- the Euler–Maclaurin (EM) constant is correct.

## 1. Main term — correct; agrees with exact data at ~1%

- **The representation is exact.** J_h = a₅H_h + K_h holds by Fubini, and the inverse Laplace transform of r₁ vanishes for x < 0.
- **Parametrisation.** The contour is u = W(1+iη). The exponent is a₅/u + (m−1/6)u − γ(5nu)/u.
- **The saddle includes the E_n part.** The −γ term is the EM integral term of log E_n(ζ^{−h}e^{−u}). By the same h-independence as ∏_c(1−ζ^{−hc}w) = (1−w⁵)/(1−w), it is identical for all four h.
- **The saddle condition is right.** α = a₅ − γ − VQ/5 is exactly Ξ′(0) = 0.
- **What is left over is O(1), not O(X).** The only part of log E_n not in the saddle is Ω_h + R_h. Ω_h is bounded, with |Ω| ≤ 0.8·Li₁(ξ), and Re Ω_h(W) = 0 holds exactly: the terms c and 5−c pair off with opposite B₁ and conjugate roots.
- **Amp keeps every factor.** Amp = W^{3/2}e^{XΞ(0)}/√(2πΞ₂) comes from (1/2π)·W·√(2π/(XΞ₂)); I checked that no factor is lost.

**Exact validation.** I built exact fmpz truncated products of S_n. There are 45 points: n ∈ {291, 437, 500, 700, 1500, 2500}, the three classes at each point, V from 2.08 to 40, including t ≈ 3.99–4.06 and μ = 0.1038/0.1051. At every point:

- the sign is correct;
- |c/Amp − M| is between 0.0005 and 0.0185;
- the certified bound 4δ + E₂ for the containing tile is between 0.27 and 0.54 (0.52 in the tail cell).

The observed error scales like 1/X. The margin is 15–30× everywhere.

| point | V | t | μ | \|c/Amp − M\| (class 0) | certified bound |
|---|---|---|---|---|---|
| n=500 | 2.08 | 105 | 0.1051 | 0.0022 | 0.296 |
| n=291 | 9.70 | 4.06 | 0.0070 | 0.0178 | 0.545 |
| n=437 | 12.0 | 3.99 | 0.0046 | 0.0143 | 0.522 |
| n=2500 | 25.0 | 5.26 | 0.0011 | 0.0052 | 0.522 |

## 2. Enclosures, NaN handling, `.mid()` and `float(` — sound as run; one latent driver defect

**What is enclosed.** Every error term is an arb or acb enclosure over the whole V-tile ball:

- the window is an upper Riemann sum with monotone factors evaluated at the right endpoints;
- sup|Ξ⁗| is taken over η-balls;
- the lower bound for the mid-region D uses a Taylor bound or a direct evaluation on a ball.

**Uses of `.mid()` and `float(`.** Each use is harmless:

- `.mid()` appears only in a sanity assertion that c₁ = 0 (c₁ is identically 0 by the definition of α) and in the diagnostic `info` output;
- `float(...)` appears only in selection keys and in the final comparison of the lower bound with 0.

One of them is not directed rounding. Line 59 stores sup|Ξ⁗| as `float(upper)`, which can round down by one ulp. That is immaterial, but `arb(...)` should be used there.

**Latent defect in the driver.** The pass test is `if min(...) <= 0`. If a margin were NaN, this test is False, and `min(worst, nan)` keeps the old value, so **a NaN tile would pass silently**. The same thing happens if `window`'s `min(..., key=)` receives a NaN. I reran all 500 tiles and the tail cell with a finiteness check: **there were 0 non-finite values**. The minimum margins are 3.0700 / 1.6879 / 0.8339, and the maximum δ is 0.13703. These reproduce the log.

- **Fix:** change the test to `if not (x > 0)`.

## 3. Euler–Maclaurin remainder — correct

I derived it independently. With P(x) = (B̃₂(x−θ) − B₂(θ))/2:

- the constant part telescopes to (B₂(θ)/2)·f′(n), and |B₂(c/5)|/2 ≤ 0.0367;
- the periodic part is bounded by sup|B̃₂|/2 = 1/12.

So |R| ≤ (1/12)(|f′(n)| + ∫|f″|) holds. Summing over c = 1..4 with f(x) = Φ_c(5ux) gives (5W/3)|ζ|ξ[...].

**The trap is avoided.** If one uses P with P(0) = 0 without splitting off the constant part, sup|P| reaches 0.12 > 1/12 at θ = 2/5. The split used in the certificate sidesteps this.

**Where things are evaluated.** The derivatives come from the closed-form majorants |Φ′| ≤ ξ/(1−ξ) and ∫_ray|Φ″| ≤ |ζ|ξ/(1−ξ). The certificate uses (1−ξ)² in the second, which is conservative. Nothing is evaluated at a single point. In `cls012_em_check.py` the bound exceeds the true remainder by a factor of 8–20.

## 4. V-only structure with X monotone — true; the stated justification is incomplete

**Terms that do not depend on X.** M, Ξ and d_k depend only on ξ and V.

**Terms that decrease in X:**

- k₃ ∝ X^{−1/2}, k₁ ∝ X^{−1/2} and k₄ ∝ X^{−1};
- R̄ ∝ W;
- the mid-region erfc argument max(cc√(d/Ξ₂), a√(Xd)) is non-decreasing in X.

The sup|Ξ⁗| range [0, η₀(X_min)] covers every X ≥ X_min.

**How m and n enter.** m enters only through α, and n only through V. The K(m) in the R-term uses m ≤ a₅X² + 1.

**Defect.** `crude()` asserts that every crude term decays at a rate of at least c_min = κ − a₅/5. That is false for the k = 10 part of the R-term, which is W·e^{a₅X/4}·Fm and decays only at the rate κ − a₅/4.

The conclusion still holds. κ − a₅/4 ≥ 0.0572 (at V = 2), so this term decreases for X > 0.5/0.0572 = 8.7, and X_min ≥ 120 in every tile. E₂ is at most 2.3e-10 in any case.

- **Fix:** assert κ − a₅/4 > 0 and X_min > 4.5/(κ − a₅/4).

## 5. Tail cell V ≥ 12 — correct

**Cauchy estimates.** On the strip |Im η| ≤ 1/2:

- Re ζ ≥ 1/2, so |e^{−Vζ}| ≤ e^{−V/2} and |1/ζ| ≤ 2, which gives |γ(Vζ)/ζ| ≤ (12/25)·Li₂(e^{−V/2});
- the discs of radius 1/2 centred on real η lie inside the strip, which gives the constants 8, 8 and 384 for Ξ₂, c₃ and Ξ⁗.

**Other bounds.**

- The a₅/ζ contributions are 2a₅, a₅ and 24a₅.
- The mid-region bound a₅η²/(1+η²) − (12/25)·Li₂(ξ) is valid.
- M ≥ |a| − 4(e^{ω̄} − 1) holds.

**Monotonicity in V.** Every V-dependent quantity is monotone in the safe direction:

- Li₁, Li₂, ξ and Vξ/(1−ξ) are all decreasing;
- X_min = 4V/a₅ grows with V.

## 6. Band, n > 2000 — correct

**Weight counts.**

- 2 parts: at most N + 1.
- 3 parts: partitions of D < N into at most 3 parts, at most (D+3)²/12 + 1/2.
- Total: 1 + N + 1 + N²/12 + … ≤ N² ≤ m² for N ≥ 10001.

**Reduction to one inequality per m.** Using N_min(m) together with the monotonicity of A in n is valid.

**Part B.**

- **Integral comparison:** Σe^{c√i} ≤ e^{c√I}(1 + 2√I/c) holds.
- **Bessel bounds:** the lower bound for P₅ via I₁ ≥ I_{3/2} ≥ e^Y(1−1/Y)/√(2πY), and the upper bound for P₁₀ via I_{1/2}, are correct.
- **Exponents:** k₁ = c(1 − √3/2), k₂ = c(1 − 1/√2), c/2 and c are correct.
- **Monotonicity:** the powers p are 1.25, 3.25, 1.75 and 0.75, and the monotonicity checks hold for them.

**Reruns.** The log values reproduce: worst correction ratio 5.06e-10, Q ≤ 5.14e-20, and the only zero is s(7) = 0.

## 7. Coverage — correct

- **The root lies above 2.** μ′ = (m − 1/6)/(25n²) < 0.4μ < 0.042. Because α is increasing, A(V) ≥ α(v₁)/v₂² on each tile, and the certified minimum is A ≥ 0.04424 on (0, 2]. So every root has V > 2.
- **X_min is right.** The condition t ≥ 4 gives αX² ≥ 4VX + 23/6, so X ≥ 4V/α ≥ 4V/a₅. The condition n ≥ 291 gives X ≥ 1455/V.
- **Tiles and roots.** The tiles are contiguous on [2, 12] and the tail cell covers V ≥ 12. The representation is exact for any W, so any root may be used; uniqueness is not needed.
- **No dependence on the other v.** The certificate never uses the v = √(a₅N/t) of §1 and §9. It is self-contained in (V, X), and it needs only the region definitions t ≥ 4, μ < 0.105 and n ≥ 291.

## Residual dependence

The only external input is the §3 expansion bound for all i ≥ 1. It is cited by transplant and checked only up to i = 2e5. That dependence is unchanged by this audit.
