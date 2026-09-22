# AUDIT2_OFFPEAK: second independent audit of the off-peak lemma (B), Lemma 7.3

*19 Sept 2026. I did the review and the numerical checks below before reading AUDIT_OFFPEAK.md or ARB_F_RUN.md.
My scripts are in the session scratchpad (scan.py, brute.py, adv.py, sound.py) and are not part of the certificate.*

## Overall verdict

**Lemma 7.3 is sound for every n >= 821.** I found no defect that changes a margin or opens a real gap.

I found three new defects, all at the level of floating-point endpoints:
- N1: an 8.9e-17 sliver at |y| = 1.1.
- N2: a 5.5e-18 sliver in L at 0.1.
- N3: the float constant 6.3 is about 1.8e-16 low.

All three are closed by a one-line continuity argument, or they are far below the smallest certified slack (0.36). The same continuity argument also makes the first audit's D1 (the sliver at π) moot.

The smallest true slack at any of the 5,570 adversarial Arb points is **43.18**, at n = 821, L = 2.1, |y| = 1.1⁺. It grows roughly linearly in n: 127 at n = 2000, 706 at n = 10⁴ and 7239 at n = 10⁵.

## 1. The mathematical ingredients, checked on paper

**(E) and (E0), Euler–Maclaurin (EM): correct.**
- The factors match: q^k = ζ^{hc} e^{-(t+c/5)σ/n} for k = 5t + c, so the sum over t < n and c = 1..4 is exactly the product S_n.
- The EM terms come out as follows:
  - The B1 terms give C_h.
  - The f′ terms carry |σ|/n.
  - The remainder is at most (1/12)∫|f″| = (|σ|²/(12n))∫₀¹|Φ″(tσ)|dt, because |B̃2|/2 ≤ 1/12.
- Φ_c is analytic on the segment [0, σ]:
  - Re(1 − w) > 0 when |w| ≤ 1 and w ≠ 1, so the principal log is continuous along the segment.
  - Only real parts are used, so the choice of branch does not matter.
- The closed form for p, using Li2 and Σ_c Li2(ζ^c x) = Li2(x⁵)/5, is correct.
- The mean-value enclosures for p(σ) and C(σ) are valid, because an acb box is convex.
- (E0) re-derived:
  - The B1 terms cancel because F does not depend on c.
  - The B2 coefficients sum to Σ|B2(c/5)|/2 = 24/300.
  - The remainder is (|ω||σ|/3)∫_{1/n}^1 |Φ″(tσ)| dt.
  - Ψ = 4 + log(24/625) + O(|ω|). The code's 6|ω| is conservative; 4|ω| would suffice.
  - The bound |Φ″(s)| ≤ 1/(|s|²m²) holds for every s with Re s ≥ 0 and m > 0, not only for |s| ≤ 0.1. That matters because the code uses it up to |s| = 0.31, and it is still valid there.
- **Numerical check.** The residual D(n) = log|S_n| − n·Re p(σ) − Re C satisfies n·D(n) ≈ 0.013, −0.0018 and −0.080 for n = 821…5·10⁴. That is O(1/n) and well inside K/n with K ≈ 0.7.

**(N), the triangle bound: correct as far as I checked.**
- Re-derived:
  - G_j = (1 − e^{−jσ})A_j/(1 − e^{−jω}).
  - |A_j| ≤ Ā_j + 2j|ω|, using |e^{−z} − 1| ≤ |z| for Re z ≥ 0.
  - |A_j| ≤ α(jL/n).
  - Ā_j is 1 or 4 for every h ≠ 0.
  - In Part I, the maximum-modulus bound c1 = 1 + 1/π holds, with the edge values checked. Re(jω) ≤ πL/y ≤ 1.1, so every z stays inside the rectangle.
  - The R1 bound works out to B = 2 + 1.6c1 (h ≠ 0) and 2 + 4c1 (h = 0), with 4πc1|σ|/y ≥ 2πc1|σ|/y.
- I did not re-derive every Part II constant on paper. Instead I compared the certified pieces with the exact j-series:
  - 40 points of (N-a) at n ∈ {821, 1500, 5000}: certified Part I and Part II are each at least the exact partial sum. The minimum excess of the total is 75.1.
  - 30 points of (N-b) at n ∈ {821, 2000, 20000}, including u = 0.25 and L = 0.1: the certified V is at least the exact triangle-bound value, with minimum excess 234.5.
  - Violations in both: 0.
- Symmetry (h, y) → (5 − h, −y), and the reduction to h ∈ {0, 1} with y > 0, are valid.

**(F), chunks and the maximum principle: correct.**
- **Factorisation.** k = 5Tm + k′ and V_m = q^{5Tm}, so |V_m| = e^{−mλ}. The drift r^{k′} equals e^{ik′θ̃}, and μ ≤ μ_up for every n ≥ N0.
- **Leftover factors.** There are at most 4(T − 1) of them, each at most 1 + e^{−(L−λ)}, since TM L/n > L − λ.
- **Lower bound for np(L).** I re-derived it with the tangent at mλ and convexity: −Q′(mλ)λ ≤ Q((m−1)λ) − Q(mλ), and the m = 0 term is at most λ² since |Q′| ≤ 2. It also needs Q > 0 and Mλ ≤ L. Numerically the minimum slack over 300 random (n, L, T) is 0.0013, so the bound is tight at small L but valid.
- **Convexity in s.** 𝓑 is convex and non-increasing in s (Hadamard). The cell bound d_i = min(a1, a2) is the exact minimum of (tangent − chord).
- **Subharmonicity.** sup_{|V|=e^{−s}} log|P| is the sup of a compact continuous family of log|entire| functions, so it is subharmonic in θ̃. The maximum over [θa, θb] × [0, μ_up] is therefore on the boundary, and all four sides are checked.
- **Assembly.**
  - For v1: M ≥ n/T − 1.
  - For v2: this is a left Riemann sum of the non-increasing running minimum. It requires all d_i > 0 up to Lb, which the code checks. That makes it legitimate to truncate at X = La(1 − T/N0) ≤ Mλ and to use 1/λ ≥ n/(T Lb).
  - Slope conditions are present on both routes.
- **Arb evaluator.**
  - 1 − 2ρc + ρ² is convex in ρ, so taking the maximum at the two endpoints is enough.
  - cos_lo is a valid lower bound.
  - Cells come from a dyadic union of exact endpoints.
  - There is no early exit in the current code.

**Lemma M: proved.**
- (1 − ρ)² ≤ (1 − ρ′)² ≤ (ρ/ρ′)(1 − ρ′)², and the sin² term scales exactly by ρ/ρ′.
- Σ_{k≤5n, 5∤k} k = 10n², which gives the factor n(L′ − L).
- Numerically, the minimum slack over 1500 random points is 1.66.

## 2. Coverage in θ

For θ ∈ [0, π] and n ≥ 821, write δ_h = θ − 2πh/5 and y = 5nδ_h:

| set | pieces | θ-set (all n ≥ 821) |
|---|---|---|
| near 0 | (E0): 0 ≤ y ≤ 20; (N, h=0): 20 ≤ y ≤ 0.25n | 0 ≤ θ ≤ 0.05 |
| near 2π/5, 4π/5 | (A): \|y\| ≤ 1.1; (E): 1.1 ≤ \|y\| ≤ 6; (N-a): 6–200; (N-b): 200–0.25n | \|δ\| ≤ 0.05 (both signs) |
| rest | (F): 8 pieces, 163 rectangles | [0.0499, TP−0.0499] ∪ [TP+0.0499, 2TP−0.0499] ∪ [2TP+0.0499, fl⁺(π)] |

- (N)'s outer end is 0.25n/(5n) = 0.05 in θ **for every n**, so the n-scaling cancels exactly. (F)'s fixed start at 0.0499 overlaps it by 10⁻⁴ rad.
- (N-a) and (N-b) meet at y = 200. This needs 0.25n ≥ 200, i.e. n ≥ 800.
- I checked in Arb that the float endpoints are 0.0499 − 9·10⁻²⁰ and 0.04990000000000000644…0.0499000000000001534 from their roots. All of these are below 0.05, so the overlap holds.
- The (F) rectangles tile each piece exactly (consecutive float endpoints are equal), and the pieces abut exactly (the same float expressions are used). The only gaps are the two peak windows. PI_UP > π.
- Symmetry: the reflection θ → 2π − θ covers [π, 2π]. It maps the roots 6π/5 and 8π/5 to 4π/5 and 2π/5, and the root 0 to 2π. The point 6π/5 is 0.63 from π.
- Numerically, h1 ↔ h4 and h2 ↔ h3 gave identical slacks at every n.

**New endpoint slivers (N1, N2).** Both are the same kind of issue as D1.
- **N1.** (E) starts at y = fl(1.1) = 1.1 + 8.9·10⁻¹⁷. The lemma needs every |y| > 1.1, so (1.1, fl(1.1)) is not covered by any certified box.
- **N2.** The L-boxes of (E), (E0) and (N) start at fl(0.1) = 0.1 + 5.5·10⁻¹⁸.
  - For L ∈ [0, 5.5·10⁻¹⁸), Lemma M up to fl(0.1) costs n·fl(0.1). The (N) boxes charge only n·arb('0.1').
  - (E) and (E0) charge arb(0.1) = fl(0.1), which is fine.
- **Fix, one sentence.** For fixed n and L, log|S_n(re^{iθ})| is continuous into [−∞, ∞), and likewise in L for fixed θ. So "≤ −E(n)" passes to closures, and every float-endpoint sliver (N1, N2 and D1) is covered.

## 3. Uniformity in n

Each piece uses "value at N0 plus a slope check", which I verified in the code:
- (E): g > 1.5/N0.
- (E0): same test, and Rem is non-increasing in n.
- (N-a): R1 = A + B log n. Part II is non-increasing in n, because M/n, D/n and log⁺(2MG/D) do not increase and β_k decreases.
- (N-b): G ≥ (b + 2c log N0)/N0.
- (F): μ_up covers every n ≥ N0, and the n-dependent terms are taken at N0.

Tests:
- **Re-running the certificates at larger N0.** Everything is certified at N0 = 2000, 10⁴ and 10⁵.
  - (N) worst V: −480 / −2762 / −28512 for h = 1 (N-a), and −112 / −4557 / −51209 for h = 0 (N-b). The certified value at 821 was −1.8.
  - (E) passes at 2000 and 10⁵.
  - (E0) passes at all three. Its worst V (−0.19, −8.5, −44) is a bisection artefact, because the search stops at the first box that passes.
- **Box bounds at the actual n.** I re-evaluated the box bounds with n substituted (n = 821, 2000, 10⁴, 10⁵) and compared them with exact Arb log|S_n| at points inside tiny boxes, 120 points in all.
  - Bound − truth ≥ 0 everywhere.
  - The minimum excess is 0.00 for (E), which is expected because EM is nearly exact. It is 0.39 for (E0), 103 for (N-a) and 248 for (N-b).
  - Some tiny boxes gave a NaN bound at nt = 32. `certify` handles this by raising nt, and these are not violations.
- **True slack along the worst corner** (L = 2.1, |y| = 1.1⁺): 43.18 (n = 821), 43.25 (822), 55.9 (1000), 127.5 (2000), 706 (10⁴), 7239 (10⁵). It increases monotonically.

## 4. The summary arb_F_summary_full.log

- **Coverage.** I parsed every `piece=` line of arb_F_p*.log and compared its endpoints with `rects_of(p)[r]`, bit for bit.
  - All 163 rectangles are present, OK, and have the right T and positive margin.
  - None FAILED, and there are no sub-rectangle splits.
  - The only duplicate is piece 5 rect 0 (p5.log and p5a.log), which is harmless and has the same margin, 80.09. p5_r1.log is empty.
  - The minimum margin is 47.91 (piece 7 rect 26, ending at PI_UP).
- **Old versus new code.** Every log written before the 09:58 code change was produced by the version with the early exit. That covers 132 rectangles: p0 r0, p1–p4, p5 r0, p5 r27–52, p6 and p7.
  - The exit was `if p > thr: return False`, so it could only reject a cell. Any cell it accepted had passed the same full-product test, so those results remain valid.
  - I re-ran 11 of these rectangles with the current code: p7 r20–26 (including the worst, 47.91), p0 r0 (51.76), p1 r4–5 and p5 r52. The margins are identical, with arbfail = 0. The old runs had arbfail > 0, which is the expected signature of false rejections.
  - The lowmem variant only pre-splits (exact tiling) and clears a cache. No log shows it was used, since every rectangle has exactly one sub line.

## 5. Independent brute-force check

I evaluated log|S_n(re^{iθ})| in Arb (96-bit) as a full product, Σ ½ log((1−ρ)² + 4ρ sin²(kθ/2)), one factor at a time.
- **Points:** 5,570 in total.
  - n = 821 (1,920 points, 12 L values).
  - n = 822, 1000, 2000 (1,120 each, 7 L values).
  - n = 10⁴ (180) and n = 10⁵ (110).
- **L values** included 10⁻⁹, 0.005, 0.05, 0.0999, 0.1, 0.1001, 0.5, 1, 1.5, 2, 2.0999 and 2.1.
- **θ was placed at:**
  - |y| = 1.1⁺, 1.2, 2, 3.5, 6± and 200±, and at u = 0.2499 to 0.2505, for all h = 1..4 and both signs;
  - h = 0 at y = 0, 10⁻⁶, 1, 5, 20± and 0.25n on both sides of 0 and 2π;
  - the (F) endpoints 0.0499, 0.1, TP ± 0.1, 2TP ± 0.1 and 3TP ± 0.05, and 2π − 0.0499;
  - π and π ± 10⁻⁹, π/5 and 3π/5;
  - the numerically located local maxima near π/5, 3π/5, π, TP⁻ and just outside the (N)/(F) handover, for every (n, L).
- In addition, a float scan of every θ ∈ [0, π] outside the windows at y-step 0.2 (n = 821) and 0.5 (n = 2000), for 8 and 4 values of L, found nothing worse than the window edges.

| region (n = 821) | min slack |
|---|---|
| window edge \|y\| = 1.1, h = 1/4 (L = 2.1) | **43.18** |
| window edge, h = 2/3 | 43.25 |
| (F): secondary peaks π/5, 3π/5 | 321.5 |
| near π | 428.1 |
| (N)/(F) handover, 0.05–0.12 | 443.7 |
| h = 0 region | 446.3 |

Violations: 0. Minimum slack: 43.18 at n = 821. The first audit reported true V = −43.2, which agrees.

## 6. Comparison with the first audit

- **Agreements.**
  - The lemma is sound, and every margin is as stated.
  - D1 (the sliver at π) and D2 (the unstated h/sign reduction) are confirmed.
  - D3 (the float sort key) is fixed in the Arb assembly; I checked the exact comparisons.
  - The worst true corner is −43.2, and the Lemma M and EM derivations match.
  - ARB_F_RUN.md, read afterwards, confirms that the early exit was `return False`, which is what I had inferred.
- **New findings in this audit.**
  - N1: the |y| = 1.1 sliver of 8.9·10⁻¹⁷ in y. It is closed by continuity.
  - N2: the L sliver at 0.1. The first audit noticed the arb('0.1') versus fl(0.1) cost difference but not the uncovered interval [0.1, fl(0.1)). It is closed by continuity.
  - N3: `E_req` uses arb(6.3), which is fl(6.3) = 6.3 − 1.8·10⁻¹⁶, in (E), (E0) and (N-a). (N-b) and (F) use arb('6.3'). The effect is negligible next to a slack of at least 0.36. The fix is `arb('6.3')`.
  - Remark: the proof should state the continuity closure once. It makes D1, N1 and N2 non-issues and removes the need for PI_UP.
  - Remark: 132 of the 163 (F) rectangles rest on the old-code runs. Their validity rests on the argument that the early exit could only reject cells, and on 11 re-runs that reproduced them exactly. A full re-run with the current code, about 1.5 h on 2 cores, would remove even this dependence, but I do not consider it necessary.
- **Nothing found against:** the (N)/(F) handover, the symmetry, the treatment of 6π/5 and 8π/5, the tiling of the summary, the subharmonicity step, or the monotonicity in n.
