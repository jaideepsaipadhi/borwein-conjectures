# Referee report: SEC3_PROOF.md (Theorem 3.1, explicit expansion of G_5)

I did this independently. I did not read GAP_CITATION.md or any earlier audit. All numerics use my own code (mpmath, 50 or more digits), in the session scratchpad.

**Overall verdict: correct.** I found no mathematical error and no missing term. The numerical constants and their rounding directions are right. The final inequality holds with no violation on n in [1, 5000] (checked exhaustively) and on 1509 further n up to 200000. The remarks below are minor presentational gaps, not errors.

## Per step

1. **Modular transformation (§3.4): correct.**
   - I re-derived (3.2) and (3.3) from Fact 3.5 (the convention hH ≡ −1 mod k, and ω = e^{πi s(h,k)}).
   - 5 | k: the substitution z_1 = z/5 checks out, since 10πz/k² = 2π(z/5)/k_1². The square roots cancel to exactly 1. The exponent is −π/12z + 5π/12z = +π/(3z). The phase is ω(h,k_1)/ω(h,k).
   - 5 ∤ k: the prefactor is √5 (same argument, so principal roots multiply). The exponent is −π/12z + π/60z = −π/(15z). The phase is ω(5h,k)/ω(h,k).
   - The shift e^{−πz/(3k²)} gives C = 2n − 1/3.
   - Independent check: 10 cusps ((1,5), (2,5), (3,10), (7,15), (1,20), (1,1), (1,3), (2,7), (3,4), (4,9)), at z = 1.3+0.4i and z = 0.2−0.35i. The maximum relative error was 4.6e−49 at 50 digits, with H computed as −h^{−1} mod k.

2. **Farey order (§3.2–3.3): correct.**
   - N + 1 = ⌊√(4πn)⌋ + 3 > √(4πn), so (N+1)² > 4πn > 2πC. That is exactly what Cor. 3.4 needs: Re z ≤ 2k²/(N+1)² gives e^{πC Re z/k²} < e.
   - N ≥ 5 holds for all n ≥ 1, so k = 5 is present.
   - Lemma 3.3 checks out: |z'| = k/√(k²+k_1²); the chord/arc geometry is right; the arc length is θ/2 ≤ (π/2)sin(θ/2) by Jordan's inequality, which applies because θ ≤ π; Re z is monotone on the semicircle.

3. **Main term (§3.5): correct.**
   - Lemma 3.7 checks out, including sign and orientation. My independent quadrature on Re w = 1 (k = 10, n = 7) gave 0.0817400249677017667…, equal to P_10 to 30 digits.
   - Completing the arc is fine: on K∖{0}, |e^{π/(3z)}| = e^{π/3} exactly. Together with e from Cor. 3.4 and total arc length ≤ π√2 k/(N+1), this gives (3.5).
   - |a_k| ≤ φ(k) is trivial.
   - I verified the closed form of a(n).

4. **Non-growing cusps (§3.7): correct.**
   - On the chord, w = Re(1/z) ≥ 1, and all three factors are non-increasing in w: |(x';x')| ≤ f(e^{−2πw}), |F(x_5'')| ≤ f(e^{−2πw/5}), and e^{−πw/15}.
   - The chord length 2√2 k/(N+1) and the 1/k² factor give φ(k)/k ≤ 1 over at most N values of k, divided by N+1, which is < 1.
   - So E_ng = 2√2 e K_5, with K_5 as stated.

5. **Missing terms: none found.**
   - The growing-cusp remainder R_{h,k} is bounded in §3.6. It is moved to the chord by Cauchy's theorem; the region between arc and chord lies in Re z > 0.
   - On the chord, w can be as large as about (N+1)²/(2k²). This is harmless: |e^{π/(3z)}| e^{−2πjw} = e^{(π/3−2πj)w} with negative exponent for all j ≥ 1, so the bound at w = 1 holds. This is where B = 1/3 < 2 is essential, and the text says so.
   - No correction term needs extracting.
   - The majorant |g_j| ≤ [y^j] f(y)f(y⁵) is valid coefficientwise.
   - Small k (k = 1, k = 5) are inside the same sums. k = 5 is the principal term; its arc and remainder errors are included in E_arc and E_rem (φ(k)/k ≤ 4/5 for every 5 | k, and there are at most N/5 such k).

6. **Constants: correct, and rounded in the right direction.** At 50 digits:

   | constant | value | stated bound |
   |---|---|---|
   | K_5 | 2.85499540499087… | ≤ 2.8549955 |
   | E_arc | 5.50644686778894… | ≤ 5.506447 |
   | E_rem | 0.00657086323166… | ≤ 0.006571 |
   | E_ng | 21.9505238422352… | ≤ 21.950524 |
   | E_con | 27.4635415732558… | ≤ 27.4636 |

7. **Numerical test of the final inequality: passes.**
   - Method: per-n precision dps = 2√(0.263n)/2.302 + 40, with π built after setting dps, and E_con = 27.4636.
   - n = 1…5000, exhaustively: max ratio LHS/RHS = 0.90451, at n = 4987. There were no violations.
   - n ≤ 400: max ratio 0.76593, at n = 397. This matches the proof's claimed 0.766.
   - 1500 random n in (5000, 200000], plus spot values up to 200000: max ratio 0.9045085, at n = 122832. There were no violations.
   - Structure of the ratios by n mod 5:
     - n ≡ 3, 4: ratio about 1e−53, because s(n) = 0 and a(n) = 0.
     - n ≡ 0: the ratio tends to 0.3455.
     - n ≡ 1: the ratio tends to 0.5590.
     - n ≡ 2: the ratio tends to 0.904508 = |a_10|/φ(10) for that class. The bound is therefore asymptotically tight within a factor of about 1.1. It can never be violated asymptotically, but the margin is small.

8. **Cited standard facts: stated correctly, and applicable as used.**
   - Apostol Thm 5.1 (Rademacher form, hH ≡ −1): confirmed numerically above.
   - Hardy–Wright Thm 30 (k + k' > N for Farey neighbours): correct.
   - DLMF 10.9.19 (Hankel loop for I_ν).
   - Ford-circle / Rademacher path (Apostol Ch. 5).

## Minor gaps (presentational, not errors)

- **(a) Fact 3.6.** Opening the Hankel loop to a vertical line should say explicitly that the connecting arcs vanish: on |s| = R with Re s ≤ c, the integrand is O(e^{c}R^{−2}), since e^{x²/(4s)} is bounded for |s| ≥ R_0. The argument is standard. My quadrature also confirms the vertical-line identity directly.
- **(b) The k = 1 arc.** The arc for h/k = 0/1 is split across τ_0 and τ_0 + 1. It should be noted that after translating by 1 it is a single Rademacher arc z' → z'' on K. This is standard, but it is not stated.
- **(c) Cauchy deformation for the non-growing cusps.** When moving the integral to the chord, one should note that the region between arc and chord avoids z = 0 (the chord lies in Re z ≥ min(Re z', Re z'') > 0). This is stated for §3.6 but only implicit in §3.7.
- **(d) Numerical margin.** For n ≡ 2 (mod 5) the true error sits at about 90% of the bound. Any future tightening of the error budget must not touch the Σ φ(k)P_k term.

**Verdict: Theorem 3.1 is proved as stated, with E_con ≤ 27.4636.**
