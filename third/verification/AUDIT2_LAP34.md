# AUDIT2_LAP34: second independent audit of the classes-3,4 Laplace certificate

*19 September 2026. Scripts are in the session scratchpad (`ex2.py`, `mymain.py`, `cmp.py`, `mono.py`, `em.py`, `eta0chk.py`). I wrote §§1–5 before reading AUDIT_LAP34.md.*

**Overall verdict: SOUND.** I found no mathematical hole and no unbounded term. Nothing in the certificate path is unsound. The only items are cosmetic ulp-level ones (§5) and a stale log (§6).

## 1. Re-derivation from scratch — matches the paper and the code

1. **Starting point.** By Lemma 2.1 and Thm 3.1, c(m) = e(m) + Σ_h A_h ζ^{hm} J_h + ρ, with J_h = Σ_{j<m} e(j) ζ^{−hj} P₅(m−j).
   - The A_h ↔ w_h re-indexing h → 5−h checks out: w₄ = e^{iπ/5} = A₁.
2. **Laplace form.** P₅ − a₅ is the inverse Laplace transform of r₁. The terms with j ≥ m have x = m − j − 1/6 < 0 and give 0, so the sum over j closes into E_n. This gives J_h = a₅H_h + K_h exactly, for every W > 0.
3. **Euler–Maclaurin.** Write k = 5t + c and f(t) = −Φ_c(5ut). Then:
   - ∫f = −(1/5u)∫_σ^∞ Q = −γ(σ)/u, using γ′ = −Q/5, which I checked;
   - −B₁f(n) gives Ω_h = Σ B₁(c/5)Φ_c(σ);
   - the remainder is (1/12)·[5|u|ξ/(1−ξ) + 5|u|·(|u|/W)·ξ/(1−ξ)] per class c, summed over 4 classes. This is ≤ R̄.
4. **Why the class 3/4 leading term cancels.** Σ_h A_h ζ^{hm} e^{Ω_h} = e^{−σ} g(σ). The O(1) term Σ_h A_h ζ^{h·cls} = a(cls) = 0, and the O(e^{−σ}) coefficient is f_∞ = −Σ_h A_h ζ^{h·cls} Σ_c B₁ ζ^{−hc} = −1.
5. **Where e^{−5nu} goes.** The factor e^{−σ} = e^{−5nu} joins the exponent: a₅/u − γ/u + (m − 1/6 − 5n)u. With u = Wζ and X = 1/W this is X[(a₅ − γ(Vζ))/ζ + (m′W²)ζ], where m′ = m − 1/6 − 5n. It equals XΞ_V iff **α(V)X² = m − 1/6 − 5n**.
6. **Stationarity.** Ξ′(0) = i[VQ/5 − a₅ + γ + α] = 0 identically.
7. **Amplitude.** (W/2π)·√(2π/(XΞ₂)) = W^{3/2}/√(2πΞ₂), which is Amp. ✓
8. **Code comparison.** My mpmath main term (my own α, Ξ₂, g) agrees with `lap34_mainchk` to 30 digits at 4 points.

## 2. Term inventory: every term is bounded

| term of the exact representation | bounded in |
|---|---|
| e(m) ≤ e^{mW}E_n(e^{−W}) | `crude`: Te |
| a₅H_h, 4 values of h | `crude`: 4a₅Fm |
| K_h, the e^{a₅/u} part, \|η\| ≤ η₀ | `window` + `gtail` (Gaussian mass outside \|s\| ≤ cc) |
| K_h, the e^{a₅/u} part, η₀ < \|η\| ≤ 2 | `mid` (tiles) + `TAILSKIP` (arg > 40) |
| K_h, the −(1 + a₅/u) part, \|y\| ≤ 2W | Tcr, 1st term: 2η₁(W + a₅) |
| K_h, all of r₁, \|y\| ∈ [2W, 1] | Tcr, 2nd term: \|e^{a₅/u}\| ≤ e^{a₅X/5} |
| K_h, all of r₁, \|y\| > 1 | Tcr, 3rd term: a₅²e^{a₅} |
| E-M remainders R_h inside G | ε_R, in window, mid and TAILSKIP |
| g(Vζ) − g(V) in the window | sup\|g′\| (ball, or Cauchy ε_g(V_a − 1)); tail: 2ε_g |
| ρ = Σ e(j)R(m−j) | Trho |

The Trho check against Thm 3.1: φ(10) = 4, Σ_{j=3}^{K/5} φ(5j) ≤ Σ 4j ≤ 2K₅(K₅ + 1), a_k ≤ a₅/9 for k ≥ 15, and K(ν) ≤ √(4π(a₅X² + VX + 1)) + 2 for every ν ≤ m. E_con = 27.4636 carries the factor e^{W/6}E_n.

The majorant E_n(e^{−W}) ≤ exp(0.8X·Li₂ + 4Li₁) was re-derived (sum over the 4 classes c). I found no term without a bound.

## 3. Monotonicity in X: no counterexample

I evaluated `tile()` at n_min = 821, 900, 1200, 2000, 5000 and 20000, on the tiles [2, 2.02], [3.5, 3.52], [7, 7.02] and [11.98, 12], for both classes. I also evaluated `tailcell(11.99)` at n_min = 821, 1000, 3000 and 20000.

- **Monotone:** every component (win, gt, mid, E₂, k₃, k₁g, ε_R) is non-increasing, and the margin rises.
  - Example [11.98, 12]: win goes 0.1186 → 0.0146 and ε_R goes 0.0424 → 0.0016.
  - Tail: win goes 0.176 → 0.053.
- **X-independent terms:** mid tiles lying inside [η₀(X), η₀(X_min)] are always summed. Their erfc argument cc√(d/Ξ₂) does not depend on X, so the argument is sound.
- **Crude terms:** decay is asserted: X_min > 4.5/c_min.

## 4. Own exact c(m), a third method — every point passes

**Method.** S_n = G₅·E_n, with G₅ = P(q)/P(q⁵), where P is the pentagonal series (nmod_poly inverse). E_n comes from two Euler identities:

- E_n = [Σ_r q^{ar}/(q;q)_r]·[Σ_r (−1)^r q^{br + 5r(r−1)/2}/(q⁵;q⁵)_r], with a = 5n + 1 and b = 5n + 5;
- the terms are evaluated mod 31-bit primes and combined by CRT, with a stability check on 2 extra primes.

This shares no code with `lap34_exact.py` or with the first audit's product trees. Checks:

- it reproduces brute force at n = 60;
- G₅ mod p agrees with `g5_coeffs.pkl` at all 200001 coefficients;
- because G₅ is computed mod p, it is not limited to the pickle's range.

**Bounds used.** The certified bound is taken from the **actual tile** the point falls in (n_min = 821, width 0.02, endpoints reproduced as the driver builds them) or from `tailcell(11.99)`.

| n | points | t | μ | V | \|c/Amp − g\| | certified | factor |
|---|---|---|---|---|---|---|---|
| 821 | 20 | 4.000–172.37 | 0.0024–**0.10500** | 2.0906–18.97 | 1.75e-4 – 6.6e-3 | 0.074–0.185 | 28–579 |
| 2000 | 8 | 4.000–33.5 | 0.001–0.008 | 9.0–29.6 | 4.2e-4 – 4.2e-3 | 0.099–0.185 | 44–347 |
| 5000 | 3 | 4.000 | 0.0004 | 46.8 | 8.3e-4 – 2.7e-3 | 0.185 | 69–225 |

**Result.** All 31 signs are negative, and every point is within its bound.

**Hard points included.**

- Pairs straddling the tile/tail junction: V = 11.9899/11.9901 and 12.0094/12.0099, at n = 821 and n = 2000.
- Class 3 and class 4 side by side at each V.
- The largest admissible μ at n = 821: m = 707739 and 707743, with M = 707743 exact, 65 primes, 1946 bits.

**Note on V near 2.** The smallest reachable V for n ≥ 821 is > 2.09. A(V) < 0.042 − 0.2/n forces V > V* ≈ 2.07, so the tiles on [2, 2.07] are slack.

**Other checks (item 2 support).**

- E-M remainder against direct summation, at n = 821 and 3000, V ∈ {2, 4, 8, 12}, η ∈ {0, 0.3, 1, 2}, and all h: max |R_h|/R̄ = 0.092.
- max |g + 1|/ε_g = 0.233 at the same points.
- A(V) is 0.04425 at V = 2 and 0.04176 at V = 2.0906, which is consistent with §6.5.

## 5. Numeric comparisons and NaN handling — sound

**NaN handling.**

- `_fin` guards the three supremum loops.
- The pass test `not all(x > 0)` rejects NaN.
- `assert dk.lower() > 0` fails on NaN.
- The `> 40` skip is False on NaN, so a NaN tile is kept and ends in a NaN margin, which is a FAIL.

**`amax`.** Its key is the float lower bound, and each candidate it chooses from is itself a valid bound, so any choice (including NaN ordering) is sound or fails.

**The windowed `min(k3 s³, 2)`.** Same argument: both entries are valid bounds.

**Coverage of the balls, verified.**

- The gsup sub-balls cover every GRID tile: 0 gaps.
- The 100 Ξ⁗ balls are contiguous, and they cover every GRID endpoint b ≤ 1.
- η₀ sits at least 8e-6 below the Ξ⁗ cutoff in all 500 tiles.

**Only cosmetic items.**

- (a) `e0 = float(eta0.upper())` rounds to nearest. It lands 1 ulp below the arb upper bound in 231 of 500 tiles, so gpsup covers [0, e0] rather than [0, η₀^arb]. This is harmless: the true η₀ = cc/√(XΞ₂(V)) ≤ cc/√(X·Xi2l), and Xi2l is below the true Ξ₂ by far more than an ulp. Suggested fix: round up (`nextafter`) for form.
- (b) The `> 40` skip compares a rounded float, so the error is at most ulp·erfc(40), which is negligible.
- (c) `.mid()` appears only in the sanity asserts Ξ′(0) ≈ 0 and Im g(V) ≈ 0. Both hold exactly by identity (§1).

**Tile endpoints.** They are floats reused as the next tile's start, so the tiles are contiguous from 2.0 to exactly 12.0, and the tail from 11.99 overlaps. GRID covers [0, 2] exactly.

## 6. Rerun

I ran `python3 lap34_cert.py 2.0 12.0 0.02` on the current (post-fix) script. It finished in 76 s with **500 tiles, 0 failures**, and minimum margins **0.8363 / 0.8361** on [11.98, 12.0]. `tailcell(11.99)` gives **0.81456** for both classes. This confirms the post-fix numbers in §6.4. The repo `lap34_cert.log` is still the pre-fix log; I did not overwrite it. The log caveat in §6.4 stands until the Appendix A.6 commands are rerun into the repo.

## 7. Comparison with the first audit (read after the above)

**Both audits found:**

- the derivation, g, ε_g, the tail strip bounds, coverage (t ≥ 4 ⇒ X > 3V/a₅; V > 2), and monotonicity are correct;
- the exact validation passes.

**The first audit found and I did not:**

- **Defect A** (the float sliver (11.99999999999983, 12)) and **Defect B** (NaN pass-through). Both were already fixed when I started. I confirmed the fixes are present and effective: the tile snap, `_fin`, and `not all(>0)`.

**I did beyond the first audit:**

- A third exact method (G₅·E_n through Euler identities, mod p, not limited to the pickle). It reaches μ = 0.10500 at n = 821 and adds n = 2000 and n = 5000 (V up to 46.8, deep in the tail cell). The first audit's n = 1000 run at μ ≈ 0.105 ran out of memory; this method ran in 7.6 min within the memory cap.
- An explicit check of the Trho constants against Thm 3.1 (φ(5j) ≤ 4j, K(ν) monotone, E_con carrying e^{W/6}E_n). The first audit deferred this to GAP_AUDIT012.
- An E-M remainder check at a second n (3000).
- A quantified float-rounding check of `e0` (§5a, benign), where the first audit said "immaterial" without a count.
- A rerun of the fixed tile driver.

**Minor difference in reasoning.** For the tail mid bound, the first audit used γ(V) + |γ(Vζ)/ζ| ≤ (11/25)Li₂. I found the linear part (γ + VQ/5)ζ has constant real part and cancels, leaving 2·(6/25)Li₂. Either is ≤ the (12/25)Li₂ the code uses.

## Verdict

| item | verdict |
|---|---|
| 1 representation / saddle | correct, re-derived independently |
| 2 term inventory | complete; no unbounded term |
| 3 monotone in X | no counterexample (6 X-levels × 4 tiles + tail) |
| 4 exact stress test | 31/31 negative, 28–579× below bound, including μ = 0.10500 and t = 4.000 |
| 5 numeric soundness | sound; ulp-level cosmetic items only |

**Overall: SOUND.** Recommended housekeeping, not needed for soundness: rerun the Appendix A.6 commands into `lap34_cert.log`, and round `e0` upward.
