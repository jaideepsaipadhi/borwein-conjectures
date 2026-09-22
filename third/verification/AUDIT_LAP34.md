# AUDIT_LAP34: adversarial audit of the classes-3,4 Laplace certificate (lap34_*)

*19 September 2026. Audit scripts are in the session scratchpad:*

- *`aud_g.py`: tests the thin function g;*
- *`aud_nan34.py`: reruns all 500 tiles with NaN guards;*
- *`aud_exact34.py`: computes exact c(m) independently, with fmpz_poly product trees;*
- *`aud_val34.py`: an independent mpmath main term, compared against the certified tile's bound.*

**Overall verdict: SOUND after two small fixes. Both fixes are applied, and the reruns give identical margins.**

- **Defect A (real, tiny).** Floating-point drift left the V-interval (11.99999999999983, 12) out of every tile, and the tail cell started at exactly 12.
- **Defect B (latent).** Two NaN pass-through paths. In the rerun no NaN actually occurs.

No mathematical error was found. The re-centred main term matches exact c(m) to 0.02–0.6% at 36 new points, including n = 821, 1000 and 1200.

## 1. Re-centring and the thin function — correct

I derived this from GAP_CLASS012 B.1–B.2:

- log E_n(ζ^{−h}e^{−u}) = −γ(σ)/u + Ω_h(σ) + R_h, with σ = 5nu.
- Σ_h A_h ζ^{hm} e^{Ω_h(σ)} = e^{−σ} g(σ) exactly, because ζ^{hm} = ζ^{h·cls}.
- The exponent is therefore a₅/u − γ/u + (m − 1/6 − 5n)u.
- Put u = Wζ and impose (m − 1/6 − 5n)W² = α(V). Then the exponent equals X·Ξ_V(η) with the same Ξ_V as cls012, and Ξ′(0) = 0 ⇔ α = a₅ − γ + Vγ′ = a₅ − γ − VQ/5. ✓
- Nothing remains besides G = e^{Vζ} Σ_h A_h ζ^{hm} e^{Ω_h + R_h}, and |G − g(Vζ)| ≤ 4e^{ω̄}e^V(e^{R̄} − 1). ✓

Checks on a(cls), f_∞ and ε_g:

- a(3) = a(4) = 0, and f_∞ = −1 to within 1e-37.
- The proof of ε_g is correct: the constant term cancels, Σ|B₁| = 0.8, and |e^Ω − 1 − Ω| ≤ |Ω|²e^{|Ω|}/2. It is decreasing in V (ξ^{−1} times a series starting at ξ²).
- Numerical test at 6000 points (V ∈ [0.27, 33], |Im s| up to 1e4, both classes), using an independent mpmath product form: max |g + 1|/ε_g = **0.233**, with no violation.
- arb gfun and gder agree with mpmath (direct formula and numerical derivative) to 12 digits.

## 2. log g kept out of the curvature — correct

The window bound

  |e^ψG − g(V)| ≤ e^{k₄s⁴}(V|η| sup|g′| + ε_R) + |g(V)|(|e^{iφ} − 1|e^{k₄s⁴} + e^{k₄s⁴} − 1)

is a valid triangle inequality. It is first-order and needs no smallness assumption.

The supremum of g′:

- It is taken over acb balls s = V_ball(1 + iη_ball) covering the whole window at X_min.
- The alternative Cauchy bound ε_g(V_a − 1)/1 is valid, because ε_g is decreasing and Re s ≥ V_a − 1 = 1 > 0 on the circle.

What remains of the amplitude is exactly G. No other O(1) factor exists:

- the e^{Vζ} is inside g;
- the E-M remainder is inside ε_R;
- (1 + a₅/u), e(m), H_h and ρ are handled in E₂.

The normalisation Amp = W^{3/2}e^{XΞ₀}/√(2πΞ₂) was re-derived. The exact data confirm it (§6).

## 3. Enclosures and NaN handling — one latent defect, fixed

What I verified:

- The window is an upper Riemann sum: e^{−s_a²/2} is taken at the left endpoint and the increasing factor at the right.
- k₃, k₄ and k₁g are upper bounds, and Ξ₂ is a lower bound.
- The erfc arguments are lower bounds.
- The mid-region d_k uses the Taylor bound or a ball evaluation.
- TAILSKIP is valid: d_min covers every tile, and there are ≤ 170 tiles.
- The E-M constant (5W/3)|ζ|ξ[1/(1−ξ) + |ζ|/(1−ξ)²] was re-derived and is correct.
- Crude decay rates: in a tile c ≥ κ − max(a₅/5, a₅/4), and in the tail c ≥ κ − a₅/3 − a₅/4. Both are asserted.

**Defect B.** Two places let a NaN through.

1. `gsup`, `gpsup` and the sup|Ξ⁗| list take `max(best, float(x.upper()))`. In python-flint, an indeterminate ball gives `nan`, and `max(0.0, nan) = 0.0`. A NaN sub-ball would therefore be dropped silently and the supremum underestimated. This is the same class of defect as in B_uniform5.
2. The driver's pass test is `if min(...) <= 0`, which is the same pass-through as in the sibling.

**Fixed** in `lap34_cert.py`:

- a `_fin()` finiteness assert on all three supremum computations;
- the pass test changed to `not all(x > 0)`;
- the worst-margin update made NaN-safe.

**Rerun** (`aud_nan34.py`, all 500 tiles, finiteness asserted everywhere): **0 non-finite values and 0 failures**. The minimum margins are 0.8363 / 0.8361, identical to the log.

Remaining uses of float and mid are harmless:

- `float(...upper())` for the suprema can round down by at most 1 ulp, which is immaterial;
- `.mid()` appears only in sanity assertions.

**Monotonicity in X, term by term** ✓:

- k₃ and k₁g ∝ X^{−1/2}, and k₄ ∝ X^{−1};
- ε_R ∝ W·e^V·ξ(...) at W_max;
- the erfc arguments η_a√(Xd) are non-decreasing, and cc√(d/Ξ₂) does not depend on X;
- the window sets (Ξ⁗ and g′) are taken at η₀(X_min) ⊇ η₀(X);
- the crude terms are C·X^p·e^{−cX} with X_min > 4.5/c, as asserted;
- K(m) grows only polynomially.

## 4. Coverage — correct, except the float gap (Defect A, fixed)

The coverage inequalities hold:

- **t ≥ 4:** m − 1/6 − 5n ≥ 15n + 23/6 > 3VX. This gives X > 3V/α > 3V/a₅, i.e. V < a₅X/3. ✓
- **n ≥ 821:** X = 5n/V ≥ 4105/V ≥ 4105/V_b. ✓
- **μ < 0.105:** A(V) = m′/(25n²) < 0.4μ < 0.042. `cls012_cov.log` gives A ≥ 0.04424 on [0.001, 2] (with α increasing) and A ≥ 0.1998 on (0, 0.001]. Every root therefore has V > 2. Any root works, because the representation is exact for all W. ✓

In all 36 audit points below, X ≥ X_min of the containing tile.

**Defect A.** The driver builds tiles with `v2 = v + 0.02` in floating point. The logged last tile is [11.98, 11.999999999999831], the loop stops at v > 12 − 1e-12, and the tail cell is certified for V ≥ 12. **V ∈ (11.99999999999983, 12) was covered by nothing.**

This has no numerical consequence, but it is a real hole in coverage. It is fixed two ways:

- the driver now snaps the last tile to V_b: the rerun tile [11.98, 12.0] has margin 0.8363 / 0.8361;
- the tail cell was rerun at V_T = 11.99: margin **0.81456** for both classes.

## 5. The tail cell, V ≥ 12 — correct

The strip bounds are valid for all real η:

- Ξ₂ ∈ 2a₅ ± 8ε, |c₃| ≤ a₅ + 8ε and |Ξ⁗| ≤ 24a₅ + 384ε, from Cauchy with radius 1/2 and |γ(Vζ)/ζ| ≤ (12/25)·Li₂(e^{−V/2});
- the mid bound a₅η²/(1+η²) − (12/25)·Li₂ uses γ(V) + |γ(Vζ)/ζ| ≤ (11/25)·Li₂;
- |g(Vζ) − g(V)| ≤ 2ε_g and |g(V)| ≥ 1 − ε_g.

Every term moves the safe way as V grows:

- ε, ε_g, Li₁, Li₂ and e^V(e^{R̄} − 1) all decrease;
- κ increases;
- e^V ≤ e^{a₅X/3} and VX ≤ a₅X²/3 hold because t ≥ 4.

X ≥ min_V max(4105/V, 3V/a₅) = 216.31 for every V ≥ 12. So the cell covers every V ≥ V_T at every admissible X.

The V = 29.6 assert fires only in a degenerate *tile*, which is never used above 12, so it is irrelevant. Rerun: 0.81462 at V_T = 12, and 0.81456 at V_T = 11.99.

## 6. Independent exact validation — passes everywhere

**Method.** Exact c(m) from truncated fmpz_poly product trees (`mul_low`) of ∏_{k ≤ 5n, 5∤k}(1 − q^k). This shares no code with `lap34_exact.py`. It was cross-checked against brute force at n = 60.

**Comparison.** I used my own mpmath implementation, with α computed from Ξ′(0) = 0 via numerical γ′, and g in product form. The certificate's bound is taken from the **actual certified tile** (width 0.02, n_min = 821) or from the tail cell, not from a degenerate tile.

| n | points (classes 3 and 4) | t | μ | V | \|c/Amp − g\| | certified bound | margin factor |
|---|---|---|---|---|---|---|---|
| 821 | 16 | 4.00–172.4 | 0.0024–**0.10500** | 2.091–18.97 | 1.8e-4 – 6.6e-3 | 0.075–0.185 | 28–579 |
| 1000 | 8 | 4.00–96 | 0.002–0.048 | 3.53–20.9 | 3.7e-4 – 6.0e-3 | 0.074–0.185 | 31–231 |
| 1200 | 12 | 4.00–58 | 0.0017–0.024 | 5.18–22.9 | 4.6e-4 – 5.5e-3 | 0.074–0.185 | 34–160 |

**Result.** All 36 points have a negative sign, the error is within the bound, and X ≥ X_min.

**Coverage of the test points.** They include:

- t = 4.00, the smallest admissible m for each n;
- μ = 0.10500 at V = 2.091;
- V from 2.09 to 22.9.

n = 1000 with μ ≈ 0.105 (m ≈ 1.05e6) ran out of memory at 7 GB and was not tested.

## Summary

| item | verdict |
|---|---|
| 1 re-centring, g, ε_g | correct; \|g+1\|/ε_g ≤ 0.233 numerically |
| 2 g′ in window error | correct; nothing O(1) omitted |
| 3 enclosures, EM, monotone in X | correct; NaN pass-through fixed; 0 NaN in rerun |
| 4 coverage | correct except float gap (11.99999999999983, 12): fixed (snap + tail at 11.99, margin 0.8146) |
| 5 tail cell | correct for all V ≥ 12 (now ≥ 11.99), all X ≥ 216.3 |
| 6 exact validation | 36/36 pass, error 28–579× below bound |

**Residual inputs** (unchanged): the §3 expansion bound, and the E_n(e^{−W}) majorant and crude-term lemmas inherited from cls012, which GAP_AUDIT012 audited.
