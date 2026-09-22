# GAP_CLASS012: classes 0, 1, 2 away from small n (PROOF.md §11 item 1)

*18 September 2026. Scripts `cls012_*.py`, logs `cls012_*.log`.*

## Summary

| Region | class 0 | class 1 | class 2 | Status |
|---|---|---|---|---|
| **(A) Band**, n > 2000, m < 4N | worst corr/\|s(m)\| = 5.1e-10 (exact, m ≤ 2e5); Q ≤ 5.2e-20 (m > 2e5) | same | same | **closed** (given §3) |
| **(B) Laplace**, n ≥ 291, t ≥ 4, μ < 0.105 (every saddle V ≥ 2, up to V = ∞) | margin ≥ 3.07 (of \|a\| = 3.618) | ≥ 1.69 (of 2.236) | ≥ 0.83 (of 1.382) | **closed** (given §3); new certificate, independent of `B_uniform5.py` |
| **(C) Centre**, μ ≥ 0.105 | rerun done: 0 unresolved, worst 0.849 | worst 0.995 | worst 0.970 | **rerun done, but NOT closed**: the centre certificate has soundness holes that affect every class (§C.2) |

"Margin" in (B) is the certified lower bound for sgn(a)·c(m)/Amp, where Amp > 0 is the main-term scale. The row minimum is taken over all tiles and over all X ≥ X_min.

**Verdict.** Classes 0, 1, 2 are now closed in the band and in the whole Laplace region. That holds uniformly in n, for every v ≥ 2 up to v = ∞, and with no monotonicity-in-t assumption. It relies only on the §3 expansion bound, which is cited, and on elementary lemmas proved below. The centre is **not** closed for any class. The requested sound rerun of `centre_final.py`'s cell routine completed with 0 unresolved cells for classes 0, 1, 2. However, auditing its lemmas while doing so found holes in the certificate itself (§C.2). The most serious is that part of each peak neighbourhood is covered by no term at all. A second is that the tail term is normalised with a lower bound of the curvature, which underestimates it by up to 7.4×. At the handover this second error alone is enough to overturn class 4 (0.053 + 0.082 > 0.126).

---

## A. Band, n > 2000, m < 4N (`cls012_band.py`, `cls012_band.log`)

**Setup.** c(m) = s(m) + Σ_{j≥N} e(j)s(m−j). Parts of E_n are ≥ N.

- For j ∈ [N, 2N), exactly one part fits, so e(j) ≤ 1.
- For j ∈ [2N, 4N), e(j) ≤ 1 + #pairs + #triples ≤ 1 + (N/2 + 1) + ((N+3)²/12 + 1/2) ≤ N² ≤ m². These are the corrected counts from GAP_BAND.md.

With A(I) = Σ_{0≤i≤I} |s(i)|, this gives

  |c(m) − s(m)| ≤ A(m − N) + m²·A(m − 2N) ≤ A(m − N_min) + m²·A(m − 2N_min),

where N_min(m) = max(10001, ⌊m/4⌋ + 1). Here m < 4N is equivalent to N ≥ ⌊m/4⌋ + 1, and A is non-decreasing. So one inequality per m covers every n > 2000 at once, which removes the need for any monotonicity in n.

**Part A, exact, m ≤ 200000.** Every class-0,1,2 coefficient s(m) has sign(a(m)) and |s(m)| > A(m − N_min) + m²·A(m − 2N_min). The one exception is s(7) = 0. That coefficient has m < N, so c(7) = s(7) = 0, which the non-strict conjecture allows. The worst correction ratio is 5.06e-10.

**Part B, m > 200000, arb.** The inputs are:

- |s(i)| ≤ e^{c√i} with c = 2√a₅. This is exact for i ≤ 2e5; beyond that the §3 bound gives |s(i)| ≤ 7.73e-5·e^{c√i}.
- Σ_{i≤I} e^{c√i} ≤ e^{c√I}(1 + 2√I/c), by comparing the sum with an integral.
- P₅(m) ≥ h₀m^{−3/4}e^{c√m}, using I₁ ≥ I_{3/2}.
- R(m) ≤ 1.3·m·e^{c√m/2} + E_con, using Σ_{5|k≤K} φ(k) ≤ K²/10 + K/2 and P_k ≤ P₁₀.

Then

  [corrections + R(m)] / ((5−√5)/2 · P₅(m)) ≤ Q ≤ 5.14e-20.

Each term has the form C·m^p·e^{−κ√m}, which is decreasing once √m > 2p/κ; this is asserted in the script.

**Conclusion.** For every n > 2000 and every m < 4N in classes 0, 1, 2, c(m) has the conjectured sign. This rests on the §3 bound and the exact coefficients.

---

## B. Laplace region (`cls012_lib.py`, `cls012_lap.py`, `cls012_cov.py`; logs `cls012_lap.log`, `cls012_cov.log`)

This is a new, self-contained certificate. It does **not** use `B_uniform5.py`, the thin-factor lemma, the M = 25 Cauchy constant, or the t_eff device. The folded identity of §6 is replaced by an unfolded, exact Bessel–Laplace representation. The twisted E_n is handled by a first-order Euler–Maclaurin expansion with an explicit remainder, which is the correct normalisation from GAP_REPR.md made exact.

### B.1 Exact decomposition

a(i) = Σ_{h=1..4} A_h·ζ^{hi}, with A₁ = e^{iπ/5}, A₄ = e^{−iπ/5}, A₂ = A₃ = 1. Then

  c(m) = e(m) + L(m) + ρ,  L(m) = Σ_h A_h·ζ^{hm}·J_h,  J_h = Σ_{j<m} e(j)·ζ^{−hj}·P₅(m−j),

where |ρ| ≤ Σ_{j<m} e(j)·R(m−j) by §3.

P₅(i) = √(a₅/x)·I₁(2√(a₅x)) with x = i − 1/6. This is the inverse Laplace transform of e^{a₅/u} − 1. Put r₁(u) = e^{a₅/u} − 1 − a₅/u, which is O(|u|⁻²). Then (1/2πi)∫_{(W)} r₁(u)e^{xu} du equals P₅ − a₅ for x > 0 and 0 for x < 0. By Fubini (the integrand converges absolutely because r₁ = O(|u|⁻²)):

  J_h = a₅·H_h + K_h,
  H_h = Σ_{j≤m−1} e(j)·ζ^{−hj},
  K_h = (1/2π) ∫_ℝ r₁(W+iy)·e^{(m−1/6)(W+iy)}·E_n(ζ^{−h}e^{−W−iy}) dy.

This is exact for every W > 0.

### B.2 Euler–Maclaurin lemma (proved; checked in `cls012_em_check.py`)

For θ ∈ (0, 1) and f with f, f′ → 0 at infinity,

  Σ_{t≥n} f(t+θ) = ∫_n^∞ f − B₁(θ)f(n) + R,  |R| ≤ (1/12)(|f′(n)| + ∫|f″|).

*Proof.* On each unit interval, g(θ) − ∫g = ∫k·g′ with k(x) = x − H(x−θ). The mean of k is B₁(θ), and k − B₁ = B̃₁(x−θ). Integrating once more gives kernel (B̃₂(x−θ) − B₂(θ))/2. The boundary terms telescope, |B₂(c/5)|/2 ≤ 1/12, and |B̃₂|/2 ≤ 1/12.

Apply it with f = Φ_c(5u·), where Φ_c(s) = log(1 − ζ^{−hc}e^{−s}) and k = 5t + c. With σ = 5nu, V = 5nW, ξ = e^{−V} and ζ = 1 + iη (so u = Wζ):

  log E_n(ζ^{−h}e^{−u}) = −γ(σ)/u + Ω_h(u) + R_h,

where

- γ(z) = [Li₂(e^{−z}) − Li₂(e^{−5z})/5]/5, using ∫_σ^∞ Q = 5γ(σ) and Q = Σ_c Φ_c = log(1 − e^{−5s}) − log(1 − e^{−s});
- Ω_h(u) = Σ_c B₁(c/5)·Φ_c(σ);
- |R_h| ≤ (5W/3)|ζ|ξ[1/(1−ξ) + |ζ|/(1−ξ)²], using |Φ′| ≤ ξ/(1−ξ) and ∫_ray|Φ″| ≤ |ζ|ξ/(1−ξ)².

Against direct summation, the remainder is 5 to 20 times smaller than the bound in all 12 test cases.

The pairs c ↔ 5−c have opposite B₁ and conjugate roots, so **|e^{Ω_h(W)}| = 1 exactly**.

### B.3 Saddle and scaling

Choose W by

  (m − 1/6)W² = α(V),  α(V) = a₅ − γ(V) − V·Q(V)/5.

Then e^{a₅/u + (m−1/6)u}·E_n(ζ^{−h}e^{−u}) = exp(X·Ξ_V(η) + Ω_h + R_h), where X = 1/W and

  Ξ_V(η) = (a₅ − γ(Vζ))/ζ + α(V)·ζ.

The function Ξ_V depends on V only, not on n or X. It satisfies Ξ′(0) = 0 identically, and its odd Taylor coefficients are purely imaginary because Ξ(−η) is the conjugate of Ξ(η). Put Ξ₂ = −Ξ″(0) > 0 and Amp = W^{3/2}e^{XΞ(0)}/√(2πΞ₂). Then

  c(m) = Amp·[M + Σ_h A_h·ζ^{hm}·e^{Ω_h}·δ_h] + E₂′,  M = Σ_h A_h·ζ^{h·cls}·e^{Ω_h(W)} ∈ ℝ.

As V → ∞, M → a(cls). This M is exactly the main term with the correct normalisation |E_n(ζ₅e^{−W})|.

**Check against exact coefficients** (`cls012_asym_check.py`). c(m)/(Amp·M) = 0.990–0.998 at n = 60 and n = 100, over V = 2.08–5.7 and all three classes. The certified bound, 4δ ≤ 0.55, is well above the true error of 0.2–1%.

### B.4 Error terms

Each error term is either independent of X or non-increasing in X, so evaluating it at X_min covers all X ≥ X_min.

**Peak window** |η| ≤ cc/√(XΞ₂), with cc = 4. Write ψ = X(i·c₃η³ + R₄) + ΔΩ + R_h. The cubic term is a pure phase, so

  |e^ψ − 1| ≤ min(k₃s³, 2)·e^r + e^r − 1,  r = k₄s⁴ + k₁s + R̄,

with k₃ = |c₃|/(Ξ₂^{3/2}√X), k₄ = sup|Ξ⁗|/(24Ξ₂²X) and k₁ = β₁/√(XΞ₂). Here β₁ = 0.8·Vξ/(1−ξ) bounds |dΩ/dη|. The supremum of |Ξ⁗| is enclosed on sub-balls of width 0.01.

**Gaussian tail.** erfc(cc/√2).

**Mid region** η₀ < |η| ≤ η₁ = 2. D(η) = Re(Ξ(0) − Ξ(η)) is bounded below by d_k·η² on 169 fixed η-tiles, using the larger of a Taylor bound, Ξ₂/2 − sup|Ξ⁗|η²/24, and a direct ball evaluation. The contribution is

  e^{2ω̄+R̄} · Σ_k √(Ξ₂/(2d_k)) · erfc(max(cc·√(d_k/Ξ₂), η_k·√(X_min·d_k))),

with ω̄ = 0.8·(−log(1−ξ)).

**Crude terms.** All of these are relative to Amp:

- the part (1 + a₅/u) of r₁ in the window;
- the outer region |y| ∈ [η₁W, 1], with |r₁| ≤ e^{a₅X/(1+η₁²)} + 1 + a₅X/η₁;
- the far region |y| > 1, with |r₁| ≤ (a₅²/2y²)·e^{a₅};
- e(m) ≤ e^{mW}·E_n(e^{−W});
- a₅|H_h| ≤ a₅·e^{(m−1)W}·E_n(e^{−W});
- the R-term, by Rankin at W, with P_k(i) ≤ W(e^{a_kX} − 1)e^{(i−1/6)W}, which follows from x^{k−1}/(k−1)! ≤ e^{xW}/W^{k−1}, and with k ≤ K(m).

All of them carry E_n(e^{−W}) ≤ exp((4/5)X·Li₂(ξ) + 4·Li₁(ξ)). They total at most

  √(2πΞ₂)·X^{3/2}·e^{−κX}·poly(X),  κ = a₅ − γ − (4/5)Li₂(ξ).

The binding exponent is κ − a₅/5 > 0, which is asserted, as is decrease in X for X ≥ X_min. The crude total is ≤ 2.3e-10.

**Condition certified.** sgn(a)·M − 4δ − E₂ > 0, with δ = ε_win + erfc(cc/√2) + ε_mid.

### B.5 Tiles, tail cell, coverage

**Range of X.** On the proof's range, n ≥ 291 gives X = 5n/V ≥ 1455/V. The condition t ≥ 4 gives αX² ≥ 4VX, so X ≥ 4V/α ≥ 4V/a₅. Hence X_min(tile) = max(1455/V_hi, 4V_lo/a₅). The binding case is V ≈ 9.78, X ≈ 148.6.

**Tiles.** V ∈ [2, 12] in 500 tiles of width 0.02: **0 failures**. The minimum margins are:

| class | margin | \|a\| |
|---|---|---|
| 0 | 3.070 | 3.618 |
| 1 | 1.688 | 2.236 |
| 2 | 0.834 | 1.382 |

The worst δ is 0.137 at V = 9.78.

**Tail cell, V ≥ 12, analytic.** Write Ξ = a₅/ζ + a₅ζ + Δ. By Cauchy estimates on the strip |Im η| ≤ 1/2, |γ(Vζ)/ζ| ≤ (12/25)·Li₂(e^{−V/2}). Also |α − a₅| ≤ Li₂(ξ)/5 + Vξ/(5(1−ξ)), and Vξ is decreasing. The tail uses X ≥ 4V/a₅ ≥ 4·12/a₅, and M ≥ |a| − 4(e^{ω̄} − 1). The margins are 3.096 / 1.714 / 0.860 for classes 0 / 1 / 2.

**Coverage** (`cls012_cov.py`). The saddle equation reads A(V) := α(V)/V² = μ′ := (m − 1/6)/(25n²) < 0.4μ.

- A is continuous, tends to 1/5 as V → 0 and to 0 as V → ∞, so a root exists.
- α is increasing, because α′ = V·q(V)/5 with q = 1/(e^V − 1) − 5/(e^{5V} − 1) > 0.
- A ≥ 0.04424 on [0.001, 2] (5858 tiles), and A ≥ 0.1998 on (0, 0.001].

So μ < 0.105 forces every root to satisfy V > 2. The measured root is V ≈ 2.082 at μ = 0.105 for n = 291, 1000 and 10⁶.

The certificate therefore covers every class-0,1,2 point with n ≥ 291, t ≥ 4 and μ < 0.105. This is exactly what PROOF.md §9 assigns to the Laplace region. The certificate also extends into μ slightly above 0.105, which gives an overlap with the centre.

It does not need the open monotonicity-in-t statement of §7, because t enters only through the lower bound X ≥ X_min, and every term is monotone in X.

**Inputs.** The §3 expansion bound for all i ≥ 1. It is cited, and checked to i = 2e5.

---

## C. Centre, μ ≥ 0.105

### C.1 The requested sound rerun

`cls012_centre.py` reuses `centre_final.py`'s `best`/`cell`, which enclose `tail_integral` over the cell ball and bisect down to 1e-4. It was run with coarser starting widths, which is harmless because bisection is automatic. Every row finished with **0 unresolved**; the worst values are the largest err/margin ratio in each row.

| class | L range | starting width | cells | unresolved | worst err/margin | time |
|---|---|---|---|---|---|---|
| 0 | [0.05, 0.12] | 0.0035 | 20 | 0 | 0.3185 | 462 s |
| 1 | [0.05, 0.12] | 0.0035 | 22 | 0 | 0.9404 | 533 s |
| 2 | [0.05, 0.12] | 0.0035 | 22 | 0 | 0.9515 | 503 s |
| 0 | [0.12, 2.081] | 0.03 | 74 | 0 | 0.8490 | 642 s |
| 1 | [0.12, 2.081] | 0.03 | 85 | 0 | 0.9952 | 749 s |
| 2 | [0.12, 2.081] | 0.03 | 86 | 0 | 0.9698 | 323 s |

Logs: `cls012_centre_lo.log`, `cls012_centre_hi.log`. The worst ratios sit near 1 only because `cell` accepts the first route with a ratio below 1; they are not physical margins.

### C.2 Why this does not close the centre (all classes, including 3 and 4)

Auditing the lemmas behind the run found the following holes.

1. **Uncovered region next to every peak.**
   - The peak integral (`eps0_*`) runs only to x₀ = nc/√g, where g is the true curvature. `tail_integral` starts at th₀ = nc/√g_tail, where g_tail = (8/15)/(5ℓr)³ is a much smaller lower bound. It then skips every grid cell whose midpoint lies within th₀ + step of a peak, with step = 2π/2000.
   - So x₀ < |θ − θ_h| < th₀ + step/2 (at least) is bounded by **nothing**. At n = 290 the numbers are:

     | ξ | x₀ | th₀ |
     |---|---|---|
     | 0.1249 | 3.9e-4 | 2.9e-3 |
     | 0.3 | 2.6e-4 | 1.3e-3 |
     | 0.6 | 1.9e-4 | 3.5e-4 |

   - Near 0 the certificate's own tail lemma is too weak to fill the gap. Its quadratic coefficient is about 200× smaller than g, and it carries a factor e^{3C₀} ≈ 12.
   - Integrating the lemma over only the skipped strip gives 1.5e-5 at ξ = 0.1249, 5e-5 at ξ = 0.3, up to 0.66 at ξ = 0.6, and about 400 at ξ = 0.9 (relative; `cls012_centre_strip.py`).
   - The strip is a fixed θ-width. It does not scale with n, so the "n = 290 is worst" argument does not cover it either.
2. **Tail normalised with a lower bound.** `err_common` multiplies the tail by √(g_tail/2π). Relative to the true Gaussian mass √(2π/g), the correct factor is √(g/2π), which is larger. The tail term is therefore underestimated by √(g/g_tail): 7.4× at ξ = 0.1249, 4.9× at ξ = 0.3 and 1.8× at ξ = 0.6.
   - At the handover, eps0 ≤ 0.053 and the corrected tail is 0.082, against margin 0.126 for class 4. **Class 4 fails at the handover** once this is fixed, even before item 1.
   - Classes 0, 1, 2 have a large margin at small ξ, so for them this item alone would not break the certificate.
3. **The tail lemma itself is unproved in the repo.** The bound −log|S(re^{iθ})/S(rζ)| ≥ (1/5)Σ s^k(1 − cos 5kθ) − C₀/|1 − re^{iθ}| with C₀ = 0.8377 (`centre_tail5b.py`, "WK Lemma 9.1–9.2 analogue"), and the companion bound for |1 − re^{iθ}| < 1/3, have no written proof or citation in NOTES.md or PROOF.md.
4. **Euler–Maclaurin constant and total variation.**
   - `data_em` uses |E| ≤ (1/2)(|F(0)| + |F(end)|) + (5/12)·TV. The first-order Euler–Maclaurin form with boundary terms gives (1/2)·TV, since |k − B₁| ≤ 1/2. A narrow bump centred on a sample point violates (5/12)·TV, so the constant is not valid for general functions.
   - The TV is sampled on 400 grid points, which gives a lower estimate, not an enclosure.
5. **Taylor remainder evaluated at the centre.** f₄ is used as the Taylor remainder constant sinh(f₄x⁴/24), but it is computed as the fourth derivative at x = 0, not its supremum over the window.
6. **Argument and modulus at the peaks.**
   - `margin_ballL` uses the arctan law, which is the n → ∞ limit, with no error term.
   - It also assumes |S(rζ)| = |S(rζ²)|, which holds only up to O(1/n).
   - Both errors are O(ℓr) = O(L/n) and can be bounded by the lemma of B.2 applied to the finite product. This is small and easy to add, but it is not in the certificate.
7. **The range L ∈ (0, 0.05) is not in the sound runs.** This is μ from 0.483 to 0.5, i.e. m within about 3% of 5n². The old driver tiled down to L = 10⁻⁶ at n = 290. The endpoint lemma needs n·L ≤ 0.3, so it cannot overlap that tiling for n > 3·10⁵.

**What would close it.** For the near-peak part, the direct analogue of section B applied to the polynomial S_n handles all of it uniformly in n and down to L = 0. The ingredients are the exact profile p(σ) = [5a₅ + Li₂(e^{−5σ})/5 − Li₂(e^{−σ})]/σ, with μ = −p′(L)/2, and the lemma of B.2 with remainder O(L/n). That leaves the far tail (Farey arcs), which needs a **proved** tail lemma (item 3). I did not attempt this, because the tail lemma is the real missing mathematics and it is shared with classes 3 and 4.

---

## Open items after this work

1. **Centre, every class:** items C.2.1–C.2.7. The substantive ones are the uncovered near-peak region, the proof of the tail lemma, the tail normalisation (which breaks class 4 at the handover), and μ > 0.483.
2. **The §3 expansion bound** is still cited (by transplant), not written out, and checked only to i = 2e5. Parts A and B depend on it for all i.

Parts A and B depend on nothing else that is open. They do not use `B_uniform5.py`, K₀, t_eff, or monotonicity in t.

## Files

| File | Contents |
|---|---|
| `cls012_band.py` / `.log` | band, n > 2000 |
| `cls012_lib.py` | γ, Q, α, Ξ and its Taylor series (arb/acb) |
| `cls012_lap.py` / `cls012_lap.log` | V-tiles on [2, 12] and the tail cell V ≥ 12 (`python3 cls012_lap.py 2.0 12.0 0.02`; `python3 cls012_lap.py tail 12`; about 8 s) |
| `cls012_cov.py` / `.log` | saddle coverage, V > 2 whenever μ < 0.105 |
| `cls012_em_check.py` | Euler–Maclaurin expansion and remainder bound against direct sums |
| `cls012_asym_check.py` | Amp·M against exact c(m) at n = 60 and 100 |
| `cls012_centre.py` + logs | the sound-driver rerun of §C.1 |
| `cls012_centre_strip.py` | size of the skipped strip under the centre's own tail lemma |
