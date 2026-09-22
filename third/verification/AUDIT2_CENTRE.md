# AUDIT2_CENTRE: second independent audit of the centre certificate (§7, mu in [0.105, 1/2], n >= 821)

*19 Sept 2026. Written before reading AUDIT_CENTRE.md; the comparison section was added afterwards.*
Scripts: `audit2c_a3.py` (item 3), `audit2c_q5.py` + `audit2c_q5exact.py` (item 5 at n = 300 against exact coefficients), `audit2c_q5big.py` / `audit2c_q5big.log` (item 5 at n = 821, 857, 900), `audit2c_off.py` (off-window sampling).

## 1. Cauchy decomposition, re-derived by hand. Verdict: CORRECT

- c(m) = (1/2pi) ∫ S_n(re^{iθ}) r^{-m} e^{-imθ} dθ. At peak h, θ = 2πh/5 + y/(5n), so dθ = dy/(5n). r^{-m} = e^{mL/(5n)} = e^{2nμL}, and e^{-imθ} = ζ^{-hm} e^{-2inμy}.
- Factor k = 5t + c: 1 − q^k = 1 − ζ^{hc} e^{-(t + c/5)σ/n}. Euler–Maclaurin with f(x) = Φ_c(xσ/n) and a = c/5 gives n p(σ) + C_h(σ) + O(1/n). The B1 term is (c/5 − 1/2)(Φ_c(σ) − Φ_c(0)), which matches C_h.
- The exponent is n p(L) − i n p'(L) y − n a y²/2 − 2inμy. The linear term vanishes iff p'(L) = −2μ.
- Each peak therefore contributes (1/2π)(1/5n) e^{np+2nμL} sqrt(2π/(na)) ζ^{-hm} e^{C_h(L)} = 𝒩 ζ^{-hm} e^{C_h}, with 𝒩 = e^{np+2nμL} sqrt(2π/(na))/(10πn). This matches the proof.
- For real L, C_{5−h} = conj(C_h) and |e^{C_h}| = 1. So the four peaks sum to Σ_{h=1,2} 2cos(A_h − 2πhm/5), with A_h = Im C_h.
- Off-peak set: |(1/2π)∫_off| ≤ (1/2π)·meas·e^{np+2nμL−E} with meas ≤ 2π. Dividing by 𝒩 gives 10πn sqrt(na/2π) e^{−E(n)} = 10π sqrt(a/2π) e^{−6.3}. The factor r^{−m} is the same e^{2nμL} in both pieces, and the n-powers cancel. This matches `centre3_comb.py` (`OFFb = 10π sqrt(Aup/2π) e^{−6.3}`).
- OFF increases with a, so using a_up is the safe direction.
- Numerical confirmation (item 5): the window quadrature, normalised by my own 𝒩, reproduces **exact** c(m)/𝒩 at n = 300 to within 4e-11 for all 35 test points (7 values of μ × 5 classes, including μ = 0.4, 0.45, 0.49 and 0.5). The normalisation, the ζ^{-hm} phases and the main term are therefore right.

## 2. Saddle map μ → L. Verdict: CORRECT (one conservative inaccuracy)

- μ = m/(10n²) is exact and p does not depend on n. So L is the exact root of p'(L) = −2μ for each (n, m), and the linear term vanishes exactly. It is not approximated.
- μ(0) = 1/2: I computed p'(0⁺) = −1 and a(0) = 2/3 from the closed form, so μ = 1/2 gives L = 0 and r = 1. At L = 0, Φ_c is regular for |y| < 2π/5 = 1.257 > 1.1. My test points at μ = 0.5 (m = 5n²) agree with the exact values.
- Handover: μ(2.1) = 0.103772 (mpmath closed form) < 0.105, so L(0.105) = 2.0814.
- Consequently the claimed worst box [2.09, 2.10] is **not needed**: it lies beyond L(0.105). The worst box actually used is [2.08, 2.09], with slack − OFF = 0.0120. The certificate is conservative here, not wrong.
- a > 0 on every box is asserted in the script. I also checked (sampled) that Q''' < 0 on [0.01, 2.1], so p'' is decreasing and μ(L) is monotone.

## 3. Thin margin at L ∈ [2.09, 2.10], recomputed in arb. Verdict: CONFIRMED, conservative

| quantity | script (log) | my arb (independent code) |
|---|---|---|
| a_up on box | 0.14158 (64 t-pieces, loose) | 0.132451 (upper Riemann sum, Q'' decreasing, 4000 pieces); mpmath a(2.09) = 0.132403 |
| OFF_up = 10π sqrt(a_up/2π) e^{−6.3} | 0.008660 | 0.008376 (0.008660 with the script's a_up) |
| M_4 lower on box | 0.122127 | 0.123719 (200 sub-boxes); M_4(2.10) = 0.123727 |
| err_up | 0.103424 | (taken from the script; realised error is about 2e-4, see item 5) |
| slack − OFF | 0.0100 | 0.0119 |

- Rounding directions: M uses `.lower()`, err uses `.upper()` and OFF uses `.upper()`; all are correct. a_up enters OFF, and the Gaussian and strip normalisations, from above.
- Defect (minor, presentation): the comparison slack > OFF is **not asserted in the script**. It is only printed, as 6-decimal midpoints, and read from the log. Parsing all 210 `COMB` lines confirms min(slack − OFF) = 0.010043 on box 2.09 and 0.012003 on box 2.08. Every box is positive, and the rounding error of the printout (≤ 5e-7) is irrelevant. Adding an `assert` on the Arb values is recommended.
- At L ≈ 0, OFF = 0.0190 against slack 0.774, so it is not an issue.

## 4. "Errors decrease in n". Verdict: CONFIRMED

`centre3_comb.py N0 210 64 b b+1` at 5 boxes (err_up / minslack):

| L-box | n = 821 | 1000 | 5000 | 1e5 |
|---|---|---|---|---|
| 0.00 | 0.2094 / 0.774 | 0.1719 / 0.812 | 0.0343 / 0.949 | 0.0017 / 0.982 |
| 0.20 | 0.1011 / 0.727 | 0.0830 / 0.745 | 0.0166 / 0.812 | 0.0008 / 0.827 |
| 1.05 | 0.0635 / 0.291 | 0.0521 / 0.302 | 0.0104 / 0.344 | 0.0005 / 0.354 |
| 1.80 | 0.0872 / 0.077 | 0.0707 / 0.093 | 0.0141 / 0.150 | 0.0007 / 0.163 |
| 2.09 | 0.1034 / 0.0187 | 0.0805 / 0.0416 | 0.0156 / 0.106 | 0.0008 / 0.121 |

- Each component (ρ_1, ρ_2 and the strip) decreases monotonically. The strip goes from 1.6e-3 to 7.8e-296 at L = 2.09.
- Structurally: t_eps ∝ e^{K/n} − 1, t_zt and t_z ∝ 1/n, the tail ∝ n^{−1/2} e^{−cn}, and the strip ∝ sqrt(n) e^{−nD + K/n}, which is decreasing for n ≥ 1/(2D) (asserted in the script).
- OFF is n-independent given E(n) = 1.5 log n + 6.3. Its validity for all n ≥ 821 is the business of the off-peak lemma, which was not re-audited here.
- No term of the combination depends on n in any other way; the saddle L depends only on μ.

## 5. Independent validation by Cauchy quadrature. Verdict: PASS, with huge real headroom

- Method: trapezoid rule in y on the four windows |y| ≤ 1.1, with step 5e-4. The factor phases are reduced exactly (hk mod 5), and the result is normalised by 𝒩.
- **Validated against exact coefficients** at n = 300 (full S_300 by fmpz product tree): 35 points, maximum |difference| 4e-11.
- n = 821, 857, 900; μ ∈ {0.105, 0.15, 0.2, 0.3, 0.4, 0.45, 0.49, 0.5}; 5 classes each, 120 points: **every sign is correct**.
  - The maximum |c/𝒩 − main| is 1.65e-3 (classes 0 and 1). For the thin classes 3 and 4 it is at most 4.7e-4.
  - Certified err is 0.06–0.21, so the certificate's error bound is 100–500 times pessimistic.
  - Thinnest point: n = 821, m = 707739 (class 4, μ = 0.105): c/𝒩 = −0.12588 against main −0.12608.
  - μ ≥ 0.4 is included (0.4, 0.45, 0.49 and 0.5), and there c/𝒩 ≈ −0.76 to −1.0 for the thin classes.
- Off-window check at n = 821 (numpy, 30k random θ plus a dense edge grid): max(log|S| − np(L)) is −417 at L = 0, −164 at L = 1, and −60.0 at L = 2.09. The required bound is −16.36. The maximum sits at the window edge |y| = 1.1.

## Own findings (before comparison)

- **No error found in the mathematics** of the decomposition, the normalisation, the saddle map or the n-uniformity.
- D-a (minor): the slack > OFF comparison is not asserted in code; it is read from printed floats. The numbers still hold, with margin 0.010.
- D-b (cosmetic): the binding box [2.09, 2.10] is outside the needed range, since L(0.105) = 2.0814. The effective worst margin is 0.0120. The loose a_up (64 t-pieces) also costs about 3e-4 in OFF.
- The "thin" margin of 0.010 is thin only relative to the certified err. The realised error is about 2e-4 against slack 0.12, so a rounding-level slip could not hide a wrong sign.

## Comparison with AUDIT_CENTRE.md (read after the above)

- **Agreement.** Both audits independently derive the same 𝒩, the same cancellation of the linear phase, and |OFF| ≤ 10π sqrt(a/2π) e^{−6.3} with the n-powers cancelling. Both find the per-box minimum of 0.0100 at L = 2.09, and both find the saddle and the n-uniformity correct.
- My n = 300 exact numbers match theirs: class 4 at μ = 0.105 gives −0.12554 here against their −0.12555, and the maximum deviation from the main term is 4.5e-3 (class 0) in both.
- **New in this audit:**
  1. Validation at the actual threshold n = 821, 857 and 900 (120 points). This uses a quadrature route that was itself validated to 4e-11 against exact coefficients.
  2. Coverage of μ ≥ 0.4 up to μ = 1/2 (L = 0), which the first audit could not test. It passes with realised error ≤ 5e-4 in the thin classes.
  3. The binding box [2.09, 2.10] lies outside the range actually needed, since L(0.105) = 2.0814. The effective worst margin is 0.0120, and 0.0119 when recomputed with my own M_4 and a_up enclosures.
  4. The slack > OFF comparison is not asserted in code (D-a).
  5. The off-window deficit at n = 821 is at least 60 at L ≈ 2.1, against the 16.36 required.
- **No disagreement.** The first audit's D1 (the uniform bound 0.018791 against slack 0.018703) is confirmed by my numbers: OFF at L = 0 is 0.01904 with a_up.

## Overall verdict

**SOUND (conditional on Lemma 7.3 (B), which was not re-audited here).**
- Per-box min(slack − OFF) is +0.0100 as logged, +0.0120 on the boxes actually needed.
- Every certified error term is monotone non-increasing in n, checked at n = 821, 1000, 5000 and 1e5.
- The real error is 100–500 times smaller than certified.
- No wrong sign occurs in 120 high-precision centre evaluations at n = 821–900, nor in 35 exact checks at n = 300.
- Recommended housekeeping:
  - assert `M_lo − err_up > OFF_up` in Arb inside `centre3_comb.py`;
  - optionally stop the L-grid at 2.09, or restate the worst margin as 0.0120.
