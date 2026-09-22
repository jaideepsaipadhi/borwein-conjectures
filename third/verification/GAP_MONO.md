# GAP_MONO: B(v,t) over all t >= t_eff (the monotonicity gap of PROOF.md s7/s9)

Scripts: `mono_cert.py` (+ `mono_cert.log`), `mono_Gl.py` (+ `mono_Gl.log`).

## Result

For every v-row, sup over t >= t_eff of B(v,t) is bounded rigorously. Each term of the corrected budget has
an explicit majorant that is non-decreasing in W = a5/(tv) on (0, W_hi], with W_hi = a5/(t_eff v_lo). The
proofs are analytic (below), and their side conditions are checked in ball arithmetic. Each majorant is then
evaluated at W_hi, with v as a ball on 16 sub-cells per row, taking the minimum over c in {2, 2.5, 3, 3.5}.
This is valid because sup_W min_c <= min_c sup_W.

| v-row | t_eff | W_hi | **sup_{t>=t_eff} B** | repr_tile3 at t_eff |
|---|---|---|---|---|
| [2.5, 3.0] | 42.57 | 0.00247 | **0.0208** | 0.081 |
| [3.0, 3.6] | 29.56 | 0.00297 | **0.0240** | 0.093 |
| [3.6, 4.5] | 18.92 | 0.00386 | **0.0304** | 0.139 |
| [4.5, 6.0] | 10.64 | 0.00550 | **0.0427** | 0.252 |
| [6.0, 7.0] | 7.82 | 0.00561 | **0.0435** | 0.148 |
| [7.0, 8.5] | 5.30 | 0.00709 | **0.0555** | 0.226 |
| [8.5, 10] | 3.83 | 0.00808 | **0.0639** | 0.225 |
| [10, 13] | 3.00 | 0.00877 | **0.0702** | 0.377 |
| [13, 16] | 2.25 | 0.00900 | **0.0804** | 0.326 |
| [16, 20] | 2.00 | 0.00822 | **0.0781** | 0.327 |

All side conditions hold with a wide margin: min over terms and cells of kappa/(p W_hi) is 9.4, and 1 is
needed. The new numbers are smaller than repr_tile3's for two reasons:

- The zone-2 sum here uses ratio-1.02 cells, not ratio-1.3 cells. Both are upper sums of a u-decreasing
  integrand.
- Monotone terms are evaluated at the single point W_hi. The earlier runs formed W from a v-ball, which
  blows up the enclosure through dependency.

As a spot check, the zone-2 value by mpmath quad at (v, t, c) = (10, 3, 3) is 0.0333, against the certified
0.0364.

## Exact dependence on W at fixed v

In these variables m = a5/W², N = v/W, and N·w' = t'·v, so **N enters nowhere except through v**. The
thin-factor certificate (eps = 1/N <= 1/1456 < 1/151), D0 (W <= 0.02) and the logEn formulas are uniform in
N. The only N-dependence left is inside C(v), through the coefficients G_l(q); see T2.

Notation: x = e^{-v}, sig = sqrt(W³/(2a5)), win = sqrt(2π)·sig, u = (time variable)/sig.

**T1. Z1 = erfc(c/√2).** Constant.

**T1'. eps_0.** Substitute x = sig·u. Then

  eps_0 = (2/√(2π)) ∫_0^c e^{-u²/2} [cosh(A u³) − 1 + sinh(B u⁴)] du,

with A = sqrt(W/(8a5)) and B = W/(4a5). This is independent of v and increasing in W. **[proved]**

The certified value is an interval upper sum. This fixes a small soundness hole in `B_ball`, whose
x-cells took only the midpoints of ball endpoints.

**T2. Z_2 = C(v)·(2/√(2π))·∫_c^{π/sig(W)} g(u,W) du,** where g = exp(−(u²/2)/(1 + u²W/(2a5))).

- g is increasing in W and decreasing in u, and g <= e^{−a5/W}.
- The upper limit π/sig grows as W decreases. Split the integral at π/sig(W_hi). For W <= W_hi,

  ∫ <= ∫_c^{π/sig(W_hi)} g(u, W_hi) du + (π/sig(W))·e^{−a5/W}.

  The second term equals π·sqrt(2a5)·W^{−3/2}·e^{−a5/W}, which is increasing for W < 2a5/3. **[proved]**
- C(v) = (Σ_{l<=10} Gball_l·x^{l−1} + Cauchy tail)/(K0/2) has no explicit N. Its implicit N-dependence is
  through the G_l(q), and Gball was enclosed only at (n, v) = (40, 4), (80, 4) and (160, 2.5).
  `mono_Gl.py` encloses max|G_l| over the entire box q = e^{−z}, z ∈ [0, 20/1456] × [−50/1456, 50/1456].
  That box contains every N >= 1456 and v <= 20 with tN <= 50, and also the limit N → ∞. The result is
  0.535, 0.417, 0.430, 0.362, 0.125, 0.259, 0.343, 0.406, 0.484, 0.000 for l = 1..10, all below Gball.
  So **C(v) is uniform in N on |arg q| <= 50/N [certified].**

**T3. e^{Dz}·Z_3 = (4/(K0 x))·√(2π)·√(2a5)·W^{−3/2}·exp(−(a5 − Li2(x))/W + 5·Li1(x)).** This is increasing
for W < 2(a5 − Li2(x))/3, which is at least 0.12 for v >= 2.5 and far above W_hi <= 0.009. **[proved; checked
per cell]**

**T4. e^{Dz}·Z_E = D0·e^{−W/6}·(e^{Dz} Z_3) <= D0·(e^{Dz} Z_3).** The bound is monotone by T3. **[proved]**

**T5. e^{Dz}·Z_R.** Fix any t' > 0; B_ball's minimum over t' is at most the value at a fixed t'. Replace the
logmain E_n(e^{−W}) factor by the *lower* bound log E_n(e^{−W}) >= (4/(5W))·Li2(x), which is proved below.
This gives

  e^{Dz} Z_R <= P(v, W)·W^{−3/2}·e^{−κ0/W}·[ECON + Σ_{10<=k<=K(W), 5|k} φ(k)·e^{25a5/(k² t' W)}],

where

  κ0 = 2a5 − a5 t' − (4/(5t'))·Li2(e^{−t'v}) − Li2(x)/5,
  P = sqrt(2a5/(2π))·exp(2Li1(y) + (5t'W/3)·y/(1−y) − log(K0 x) + 5Li1(x)),  y = e^{−t'v}.

P is increasing in W. The number of k-terms, K(W) = int(sqrt(4πa5)/W) + 2, **grows** as W → 0. To remove
it, use 1[k <= K] <= (K/k)³ with K <= c_K/W for k > k1 = 60. Then Σ_{k>k1, 5|k} φ(k)/k³ <= 1/(5k1).

The result is a finite sum of terms const·W^{−p}·e^{−κ_i/W}, with (κ_i, p) = (κ0, 3/2), (κ0 − 25a5/(k² t'),
3/2) for k = 10..60, and (κ0 − 25a5/(65² t'), 9/2). Each term is increasing on W < κ_i/p. **[proved]** The
script checks κ_i/p > W_hi for every term and cell. The best t' is taken from {0.30, 0.35, ..., 1.50}.

*Lower bound used above.* log E_n(e^{−W}) = Σ_j (x^j/j)·h(jW), where
h(y) = (1 − e^{−4y})/((1 − e^{−y})(1 − e^{−5y})). The bound needs h(y) >= 4/(5y):

- For y <= 0.3, the alternating Taylor bounds give
  5y·(4 − 6y + 7y² − 6y³) − 4·(5y − 12.5y² + 125y³/6) = y²(20 − 48.3y − 30y²) >= 0.
- For y >= 0.8, h >= 1 >= 4/(5y).
- On [0.3, 0.8], ball arithmetic gives h − 4/(5y) >= 0.77.

## Findings (soundness notes on the old Z_R and zone 2)

1. **The k-cap at 200 in `B_ball`'s Z_R is unsound** for every Laplace row. K(W_eff) is between 207 and
   about 730. `tailcell2.py` already removed the cap; `repr_tile3.py`/`B_ball` did not. `mono_cert.py`
   has no cap. The effect is negligible: e^{Dz} Z_R <= 1e-3 everywhere.
2. **`B_ball`'s logmain uses the upper-bound formula for log E_n(e^{−W})** where a lower bound is needed. The
   formula overstates the value by about 1.2·Li1(x) + O(W). This is fixed in T5.
3. **Z_2 integrates C(v) over the whole period** (x up to π), but G_l is enclosed only for |arg q| <= 50/N.
   This matters for B_ball and for this certificate alike. It is not a monotonicity issue, but it is an
   unverified hypothesis of the zone-2 bound. **[open, pre-existing]**
4. G_10 enclosures read 0.000 (the known indexing artifact). Gball_10 = 0.6 is used regardless, and
   x^9 < 1e-9, so this is immaterial.

## Consistency of t_eff with coverage (PROOF.md s9)

A point in row [v_lo, v_hi] with t >= 4 has t >= a5·1456/v² >= t_cov(row). The rows are of two kinds:

- **v < 6.** Here t_eff = t_cov > t_min. The values are 42.57, 29.56, 18.92 and 10.64, and t_cov was
  floored to 0.01 using a5.lower(), so t_eff <= t_cov.
- **v >= 6.** t_eff is 7.82, 5.30 and 3.83 (t_cov) on the rows from 6 to 10. It is 3.00, 2.25 and 2.00
  (t_min) on the rows from 10 to 20. Each is <= max(4, t_cov) <= the point's t.

In every row t_eff is therefore at most the smallest t a covered point can have. Using t_eff below 4 is
conservative, since B is monotone. The range t → ∞ (W → 0) is included, so the μ = 0.105 upper end needs no
separate treatment.

## Verdict

**The monotonicity gap is closed for v ∈ [2.5, 20]**, given B as defined (the repr_tile3 budget with the
e^{Dz} corrections). The bound is sup_{t>=t_eff} B <= 0.0804, attained on [13, 16]. The pieces:

- Proved analytically: monotonicity of every term in W, the lower bound on log E_n, and the endpoint handling
  in Z_2.
- Certified: the side conditions, the evaluations at W_hi, and the N-uniformity of C(v) on |arg q| <= 50/N.

Still open, and not part of this gap:

- whether C(v) is valid on zone 2 beyond |arg q| = 50/N (finding 3);
- the K0 = 0.25 thin factor (certified by thin_cert.py) and the other inputs taken as given;
- the v >= 20 tail cell (tailcell2.py, its own joint argument);
- B_uniform5.py (classes 0–2).
