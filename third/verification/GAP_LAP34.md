# GAP_LAP34: classes 3, 4 in the Laplace region, rebuilt on the true (thin-factor-shifted) saddle

*19 September 2026. Scripts `lap34_*.py`, logs `lap34_*.log`.*

**Scope.** Classes m ≡ 3, 4 (mod 5), n ≥ 821, t = m/N ≥ 4, μ = m/(10n²) < 0.105. The case n ≤ 820 is exact (PROOF §15).

**Result.** Certified: c_n(m) < 0 at every point of the scope. The certificate is an arb tiling of V ∈ [2, 12] (500 tiles, 0 failures) plus one analytic cell for V ≥ 12. Its worst certified margin is **0.8146** in the tail cell, and **0.836** on the tiles, measured against a main amplitude |g(V)| ≈ 1. Relative to the main term, the certified bound on the error is at most 0.19. The validation below compares this with exact c(m) at 70 points (n = 150, 300, 821, 2000). The true relative error is 0.02% to 1.6%, which is 24 to 435 times below the certified bound, and every sign is negative.

What this rests on:

- the §3 expansion bound (cited and written out, `SEC3_PROOF.md`);
- the exact Bessel–Laplace decomposition and the Euler–Maclaurin lemma of GAP_CLASS012 §B.1–B.2;
- the saddle-coverage lemma of `cls012_cov.py`.

It does **not** use `mono_cert.py`, `zone2_cert.py`, `tail3_cell.py`, `B_uniform5.py`, K_0, the Cauchy constant M, t_eff, or monotonicity in t.

---

## 1. Why the old certificate failed, and what changes

`mono_cert.py` put the saddle at m = a₅/W², the kernel-only saddle. For classes 3 and 4 the amplitude of the twisted sum is *thin*: Σ_h A_h ζ^{h·cls} = a(cls) = 0. The surviving part is proportional to e^{−5n·u}, which comes from the first part of E_n. That factor is not slowly varying. It is an O(N) term in the exponent (the m → m − N shift of GAP_SMALL), and treating it as an amplitude is what put the old main term off by up to e^{13}.

**The fix is exact and simple.** Move e^{−5nu} into the exponent. Everything that remains in the amplitude is then a genuinely O(1), slowly varying function, the *thin function*

  g(s) = e^{s} · Σ_{h=1..4} A_h ζ^{h·cls} exp(Ω_h(s)),  Ω_h(s) = Σ_{c=1..4} B₁(c/5) log(1 − ζ^{−hc} e^{−s}),  s = 5nu = Vζ.

**Proved facts about g:**

1. Expand exp Ω_h = 1 + Ω_h + …. The constant term cancels because a(cls) = 0. The first-order coefficient is

   f_∞ = −Σ_h A_h ζ^{h·cls} Σ_c B₁(c/5) ζ^{−hc} = **−1 exactly** for both classes (checked in arb to 1e-37, `lap34_lib.finf`).

   This is the same f_1(1) = −1 as in `tail3_K0.py`.
2. For every complex s with Re s = V:

   |g(s) − f_∞| ≤ ε_g(V) := 4e^V[0.8(Li₁(ξ) − ξ) + ω̄²e^{ω̄}/2],  ξ = e^{−V}, ω̄ = 0.8 Li₁(ξ).

   This uses Σ_c|B₁| = 0.8, |Ω_h − Ω_h^{(1)}| ≤ 0.8 Σ_{r≥2} ξ^r/r and |e^Ω − 1 − Ω| ≤ |Ω|²e^{|Ω|}/2. The bound ε_g decreases in V, since it is ξ^{−1} times a series in ξ with non-negative coefficients starting at ξ². For example ε_g(12) = 1.8e-5.
3. g(V) is real for real V, because h ↔ 5 − h are conjugate pairs. For instance g(2) = −1.0773 (class 3) and −1.0118 (class 4), and g → −1 as V grows. The certificate uses only the arb enclosure of g(V) on each V-tile, which is ⊂ [−1.13, −0.96].

So in the re-centred variables the "thin factor" becomes an amplitude ≈ −1. There is no e^{−v} loss and no K_0.

## 2. Exact representation and the re-centred saddle

This is the unfolded identity of GAP_CLASS012 B.1, valid for every class:

  c(m) = e(m) + Σ_h A_h ζ^{hm}(a₅ H_h + K_h) + ρ,  K_h = (1/2π) ∫_ℝ r₁(u) e^{(m−1/6)u} E_n(ζ^{−h}e^{−u}) dy,  u = W + iy.

Here |ρ| ≤ Σ_j e(j) R(m−j) by §3. The representation is exact for every W > 0, with no folding and no D₀ constant. The folded identity's D₀ = 1.12601 is therefore not needed.

Euler–Maclaurin (B.2) gives, for every Re u > 0 and with no restriction on arg u,

  log E_n(ζ^{−h}e^{−u}) = −γ(σ)/u + Ω_h(σ) + R_h,  σ = 5nu,  |R_h| ≤ (5W/3)|ζ|ξ[1/(1−ξ) + |ζ|/(1−ξ)²].

This bound was spot-checked against direct summation at n = 821 for η up to 2 (`lap34_emchk.py`): the true remainder is 5 to 20 times below the bound.

**Saddle.** With X = 1/W, V = 5nW and u = Wζ, ζ = 1 + iη, define W by

  **α(V)·X² = m − 1/6 − 5n**,  α(V) = a₅ − γ(V) − V·Q(V)/5.

This is the old saddle equation with m replaced by m − 5n, and it is the stationary point of every O(N) part of a₅/w + (m − 1/6)w + log(−F_m):

- the kernel;
- the twisted-product part −γ/u, which is the R_P / e^{Dz} factor of GAP_SMALL, now inside the exponent exactly;
- the thin factor −5nu.

The exponent becomes

  a₅/u − γ/u + (m − 1/6)u + [Ω_h + log of thin part] = X·Ξ_V(η) + log g(Vζ),

where Ξ_V is exactly the function of `cls012_lib.py`, with Ξ_V′(0) = 0. Hence

  **c(m) = Amp·[g(V) + Δ]**,  Amp = W^{3/2} e^{XΞ_V(0)} / √(2πΞ₂) > 0,  Ξ₂ = −Ξ_V″(0).

The certificate proves |Δ| < −g(V), so c(m) < 0.

**What is and is not in the Gaussian.** The Gaussian's curvature XΞ₂ contains every O(N) contribution of log F_m: kernel, twisted product and thin shift. The remaining log g(Vζ) is O(1). Its contribution to the curvature is V²(log g)″ = O(e^{−V}), which is relative O(1/X). It is *not* put into the Gaussian. Instead it is enclosed rigorously in the window error through a bound on sup|g′|.

So the true saddle of the full Φ differs from ours by O(V|g′/g|/(XΞ₂)) in η, and that difference is paid for inside ε_win. No interval Newton step is needed, because the saddle equation involves only α(V), and the root's existence and location are handled by the coverage lemma (§5).

## 3. Error terms, all relative to Amp, all non-increasing in X

Let s = η√(XΞ₂) and η₀ = cc/√(XΞ₂), with cc = 4.

**Window |η| ≤ η₀.**

  e^{X(Ξ(η) − Ξ(0))} G(η) = e^{−s²/2} e^{ψ} G,  ψ = X(c₃η³ + R₄),

where c₃ is purely imaginary and |R₄| ≤ sup|Ξ⁗|η⁴/24, with the supremum taken over the whole window (sub-balls of width 0.01). Then

  |e^ψ G − g(V)| ≤ e^{k₄s⁴}·(V|η|·sup_win|g′| + ε_R) + |g(V)|·(min(k₃s³, 2)e^{k₄s⁴} + e^{k₄s⁴} − 1).

In this bound:

- G(η) = e^{Vζ} Σ_h A_h ζ^{hm} e^{Ω_h + R_h};
- ε_R = 4e^{ω̄}e^V(e^{R̄} − 1) bounds |G − g(Vζ)|, i.e. the Euler–Maclaurin remainders. Since e^V R̄ = O(W), it is small;
- sup|g′| over s ∈ V(1 + i[−η₀, η₀]) is the smaller of a ball evaluation (on sub-balls of width ≤ 0.05 in Im s) and the Cauchy bound ε_g(V − 1)/1.

The upper Riemann sum is taken over 800 s-cells. The Gaussian truncation contributes |g|·erfc(cc/√2).

**Mid region η₀ < |η| ≤ 2.**

  Re(Ξ(0) − Ξ(η)) ≥ d_k η²

holds on 170 fixed η-tiles, using the better of a Taylor bound and a direct ball evaluation (as in cls012). On each tile, sup|G| ≤ sup|g(Vζ)| (ball evaluation on 4 sub-balls) + ε_R(|ζ|). The contribution is Σ_k sup|G|·√(Ξ₂/(2d_k))·erfc(max(η₀, a_k)√(X d_k)).

Tiles whose erfc argument exceeds 40 are bounded jointly with the crude modulus 4e^{V+ω̄+R̄}·170·√(Ξ₂/2d_min)·erfc(40).

**Beyond |η| > 2, and the non-integral terms** ("crude", E₂). These are:

- the (1 + a₅/u) part of r₁;
- the region |y| ∈ [2W, 1];
- the region |y| > 1;
- e(m);
- a₅|H_h|;
- the R-term by Rankin at the same W, with k ≤ K(m) and K computed from m − 1/6 = αX² + VX.

Each is bounded with |E_n(ζ^{−h}e^{−u})| ≤ E_n(e^{−W}) ≤ exp(0.8X·Li₂(ξ) + 4Li₁(ξ)) and normalised against the same Amp. The gap between the real-axis and twisted products (the e^{Dz} factor) is therefore explicit in the rate

  κ = a₅ − γ(V) − 0.8 Li₂(ξ).

Relative to cls012 there is an extra factor e^V, because Amp now carries e^{−V} through the shifted saddle.

- Tiles use e^{V_hi}.
- The tail uses V ≤ a₅X/3 (from t ≥ 4), so the rate becomes κ − a₅/3.

Every crude term has the form C·X^p·e^{−cX} with p ≤ 4.5 and c ≥ κ − rate − max(a₅/5, a₅/4). Both c > 0 and X_min > 4.5/c are asserted. (Side note: `cls012_lap.py` asserts c ≥ κ − a₅/5, but its Rankin term decays only at κ − a₅/4. Its conclusion is unaffected, since X_min ≫ 4.5/(κ − a₅/4), but the assert should read a₅/4.)

In total E₂ ≤ 1e-8 everywhere.

**Monotonicity in X.** k₃ ∝ X^{−1/2}, k₄ ∝ X^{−1}, the g′ term ∝ X^{−1/2}, ε_R ∝ W, and the erfc arguments grow with X. The window shrinks as X grows, so a supremum over the X_min window covers every larger X. So every term is non-increasing in X, and one evaluation at X_min per V-tile covers all X ≥ X_min, i.e. all n and m mapping to that V-tile. This is the uniformity in n. It needs **no monotonicity in t**.

## 4. Tiles and the tail cell (`lap34_cert.py`, `lap34_cert.log`)

**Range of X.**

- X = 5n/V ≥ 4105/V, since n ≥ 821.
- t ≥ 4 gives m − 1/6 − 5n > 15n = 3VX, so αX > 3V, hence X > 3V/a₅.

Per tile, X_min = max(4105/V_hi, 3V_lo/a₅). For the tail, the minimum over V of max(4105/V, 3V/a₅) is X* = √(15·821/a₅) = 216.3.

| region | cells | failures | min margin −g − \|Δ\| (class 3 / class 4) | worst cell |
|---|---|---|---|---|
| V ∈ [2, 12], width 0.02 | 500 | 0 | **0.8363 / 0.8361** | V ∈ [11.98, 12]: win 0.119, mid 0.005, ε_R 0.042 |
| V ≥ 12 (analytic) | 1 | 0 | **0.8146 / 0.8146** | ε_g = 1.8e-5, win 0.176, mid 0.009, E₂ = 3.5e-8 |

Sample tiles: at V = 2.0 the margins are 0.936 / 0.871, with −g = 1.034 / 0.969 as arb lower bounds on the V-ball. At V = 5.0 they are 0.891 / 0.887, and at V = 10.0 they are 0.853 / 0.853.

**The tail cell** uses the Cauchy-strip bounds of cls012 for Ξ on |Im η| ≤ 1/2:

- Ξ₂ ∈ 2a₅ ± 8ε;
- |c₃| ≤ a₅ + 8ε;
- sup|Ξ⁗| ≤ 24a₅ + 384ε, with ε = (12/25)Li₂(e^{−V_T/2});
- |g(Vζ) − g(V)| ≤ 2ε_g(12), and |g(V)| ≥ 1 − ε_g(12).

Runtime is about 6 min for everything (`python3 lap34_cert.py 2.0 12.0 0.02; python3 lap34_cert.py tail 12`).

## 5. Coverage

Take any (n, m) with n ≥ 821, class 3 or 4, t ≥ 4 and μ < 0.105. Then:

1. m′ := m − 1/6 − 5n > 15n > 0.
2. The saddle equation reads A(V) := α(V)/V² = m′/(25n²) < (m − 1/6)/(25n²) < 0.4μ < 0.042.
3. By `cls012_cov.py` (arb), A ≥ 0.042 on (0, 2]. A is continuous, with A → 1/5 as V → 0 and A → 0 as V → ∞. So a root exists (IVT) and **every root has V > 2**.
4. Pick any root. Then W = V/(5n) makes the representation of §2 exact.
5. The point's X = 5n/V is at least the X_min of its tile (or the tail cell's X*), by §4.

So every point of the scope is covered, with no t_eff and no monotonicity in t. The upper end is also covered: the tail cell holds for all V ≥ 12 up to ∞.

## 6. Validation against exact coefficients (`lap34_exact.py`, `lap34_validate.py`)

Exact c(m) are computed by multi-modular arithmetic with CRT, using as many 61-bit primes as a rigorous size bound requires. They are cross-checked against direct integer expansion at n = 60.

For each point the script computes:

- the true c(m)/Amp;
- the certified main term g(V);
- the true error |c/Amp − g|;
- the certificate's bound |Δ|, from `tile` at the point's own V with nmin = n, so X_min = the point's X.

| n | points | V range | true \|c/Amp − g\| | certified \|Δ\| | all signs − | bound holds |
|---|---|---|---|---|---|---|
| 150 | 16 (t = 4.0 … 28, μ ≤ 0.095) | 2.3 – 8.1 | 1.2e-3 – 1.6e-2 | 0.25 – 0.50 | yes | yes |
| 300 | 18 (t = 4.0 … 56) | 2.3 – 11.5 | 6.0e-4 – 1.1e-2 | 0.14 – 0.28 | yes | yes |
| **821** | 20 (t = 4.0 … 172, μ up to 0.1050) | 2.09 – 18.97 | 1.8e-4 – 6.6e-3 | 0.074 – 0.185 | yes | yes |
| **2000** | 16 (t = 4.0 … 37) | 8.5 – 29.6 | 4.1e-4 – 4.2e-3 | 0.049 – 0.128 | yes | yes |

(The n = 150 and n = 300 rows lie outside the proof's range n ≥ 821. They are included because they are cheap and stress smaller X.)

The error scales like 1/X, halving from n = 150 to n = 300. This is the expected size of the first uncorrected term (the V²g″ and cubic-squared corrections). The certified bound is 24 to 435 times larger than the true error. For V ≥ 12 the bound shown is the tail cell's |Δ| (evaluated with nmin = n), because that is what the certificate uses there. A degenerate *tile* at V = 29.6 trips an internal assert; tiles are not used above V = 12.

Compare the old kernel-only Lc of GAP_SMALL, where I/Lc ranged from 0.077 down to e^{−13}.

## 7. What remains open / caveats

- **Inherited inputs.** The §3 expansion bound, and the Euler–Maclaurin lemma B.2 of GAP_CLASS012 with its remainder constant (proved there, spot-checked here at n = 821). The coverage lemma of `cls012_cov.py`.
- **Not independently audited yet.** Like every certificate in this repo, this is single-author code. The end-to-end validation against exact data (§6) is the main evidence that no normalisation factor is missing. The main term matches exact c(m) to within 1.6% at every tested point, including the points where the old certificate was off by up to e^{13}.
- **What this supersedes.** PROOF §15 item 1 (Laplace region, classes 3–4, n ≥ 821) is closed by this certificate. It replaces `mono_cert.py`, `zone2_cert.py` and `tail3_cell.py` for this region, and makes the question of monotonicity in t moot. The withdrawn §13/§14 entries stay withdrawn; they are simply no longer needed.
- **Still open elsewhere** (unchanged): the centre off-peak lemma for n ≥ 821 (PROOF §15 item 2), and the re-audit of `cls012_lap.py` (§15 item 3). On the latter, this work reused its lemmas and found only the harmless assert constant noted in §3.

## Files

| file | contents |
|---|---|
| `lap34_lib.py` | thin function g(s), g′(s) (acb series), f_∞ |
| `lap34_cert.py` / `lap34_cert.log` | V-tiles on [2, 12] and the tail cell V ≥ 12 |
| `lap34_mainchk.py` | float (mpmath) main term Amp·g(V) at a given (n, m) |
| `lap34_exact.py` / `lap34_exact_{150,300,821,2000}.log` | exact c(m) by multi-modular CRT |
| `lap34_validate.py` / `lap34_validate.log` | true error vs certified bound, and sign |
| `lap34_emchk.py` | Euler–Maclaurin remainder spot check at n = 821, complex u |
