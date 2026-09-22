# AUDIT_OFFPEAK: adversarial audit of the off-peak lemma (B), plus an Arb version of piece (F)

*19 Sept 2026. Scope: GAP_OFFPEAK.md, offpeak_common/em/em0/tb/chunk/crosscheck.py and their logs.
Pieces (A) (GAP_CENTRE2) and the n <= 820 exact runs are out of scope.*

## Verdict

**Sound after stated (cosmetic) fixes, with one open item: the full Arb run of (F) has not been executed.**
- I found no defect that changes a margin.
- I found four defects, all at the level of bookkeeping or documentation (D1–D4 below).
- The Arb version of (F) reproduces the float margins exactly on 7 of 163 rectangles, and every arb-certified cell there agrees with the float certificate.
- The full Arb run has not been executed. See section 3 for the commands.

## 1. Audit items

**(a) Coverage in theta: no gap except a ~1e-16 sliver at theta = pi (D1).**
- **Handovers.** Each pair of neighbouring pieces meets or overlaps:

  | handover | how they meet |
  |---|---|
  | (A) to (E) | abut at \|y\| = 1.1; the statement excludes \|y\| <= 1.1 with a strict >. |
  | (E) to (N) | abut at \|y\| = 6. |
  | (E0) to (N) | abut at y = 20. |
  | (N) to (F) | overlap. (N) reaches \|theta - 2pi h/5\| = 0.25n/(5n) = 0.05 exactly, for every n. (F) starts at 0.0499 (float 0.0499 and float 2pi/5 ± 0.0499 are within 1e-16 of exact). The overlap is 1e-4 rad, i.e. 0.0005n in y, the same for every n and every h. |
  | (N-a) to (N-b) | meet at y = 200. This needs 0.25n >= 200, i.e. n >= 800, which holds. |

- **Grid endpoints.** Every L-grid ends at float 2.1, which is >= 2.1 exactly (checked with Fraction). The y-grids end exactly at 6 and 20.
- **Symmetry.** |S(conj q)| = |S(q)| holds because S has real coefficients, so theta -> 2pi - theta is valid. It maps [pi, 2pi] onto [0, pi], and h = 3, 4 onto h = 2, 1 with y -> -y.
  - (E) runs h = 1, 2 with both signs of y.
  - (E0) runs y >= 0 only. This is correct, since the conjugate of theta near 0 is theta near 2pi.
  - In [0, pi], (F) excludes neighbourhoods of 0, 2pi/5 and 4pi/5. The root 6pi/5 is 0.63 away from pi, so it needs no exclusion.
- **D1.** The last (F) rectangle ends at the float `math.pi`, which is 1.2e-16 below pi. Conjugation maps (float pi, pi] to itself, so that sliver is not covered by any piece. The fix is to end at `nextafter(pi, 4)`, which offpeak_chunk_arb.py does. There is no numerical consequence.
- **D2 (documentation).** offpeak_tb.py runs only h = 1 and h = 0, with y > 0, although the doc says "h = 0..4". The code is still correct:
  - The TB majorant depends on h only through Abar_j = |sum_c zeta^{hcj}|. This is (1 if 5 ∤ j, 4 if 5 | j) for every h = 1..4, and |A_j| <= Abar_j + 2j|omega| likewise does not depend on h.
  - It depends on y only through |1 - e^{-j sigma}|, |sigma| and |y|, all of which are invariant under y -> -y.
  - So h = 1 with y > 0 covers h = 1..4 and both signs. This argument should be stated in GAP_OFFPEAK section 5.

**(b) Uniformity in n: sound. Every piece uses "value at N0 plus a monotonicity check", never a ball in n.**
- **(E).** V(n) = -n g + Re C + K/n + 1.5 log n + 6.3.
  - g, C and K depend on sigma = (L, y) only.
  - dV/dn = -g - K/n^2 + 1.5/n <= 0 when g >= 1.5/n, and this is checked.
  - I checked that the finite EM identity really is n-free in this form: it is the shifted EM for g(x) = Phi_c(x sigma/n) on [0, n] at the points t + c/5.
    - The B1(c/5) boundary terms give C_h.
    - The f' terms carry |sigma|/n.
    - The remainder is sup|B~2|/2 · int_0^n |g''| = (|sigma|^2/(12 n)) int_0^1 |Phi''(t sigma)| dt.
- **(E0).** Psi, EMb and the 1/(3m^2) term are bounded n-free using |omega| <= |sigma|/N0. The tail term decreases like 1/n. Same slope check.
- **(N-a).** R1 = A + B log n, and Part II is non-increasing in n: y/n · M and y/n · D are n-free up to +y/n terms, and beta_k decreases in n. The check is g >= (B + 1.5)/N0.
- **(N-b).** The explicit -nG + a + b log n + c log^2 n form, with the check G >= (b + 2c log N0)/N0.
- **(F).**
  - mu = L/(5n) <= 2.1/(5N0), so the mu-range of the rectangles covers every n.
  - The Riemann/convexity steps hold for every n.
  - In the assembly, the n-dependent terms T^2 Lb/n and 1 - T/n are evaluated at N0, which is the conservative end. The slope condition W/(T Lb) >= 1.5/N0 (or d_min/T >= 1.5/N0) is checked.

**(c) Lemma M: proved, and the constant is exactly n.**
- |1 - rho e^{ia}|^2 = (1 - rho)^2 + 4 rho sin^2(a/2).
- For rho' <= rho <= 1, (1 - rho)^2 <= (1 - rho')^2 <= (rho/rho')(1 - rho')^2.
- So each factor's squared modulus grows by at most rho/rho' = e^{k(L'-L)/(5n)}.
- Summing over k <= 5n with 5 ∤ k gives (1/2)(L'-L)/(5n) · 10n^2 = n(L' - L), using sum k = 25n(5n+1)/2 - 5n(n+1)/2 = 10n^2.
- For L < 0.1, y is independent of L and p(L) >= p(0.1), because p is decreasing: p' = (Q - p)/L and Q is decreasing. So the charge of 0.1n on boxes touching L = 0.1 is correct.
- A negligible quibble: offpeak_tb uses arb('0.1') for the charge, while the reduction point is the float 0.1, which is 5.5e-18 larger. The cost is about 5e-15, against a slack of at least 1.8.

**(d) (E) past |sigma| = 2pi/5: the claim is correct, and K is included.**
- Phi_c(s) = log(1 - zeta^{hc} e^{-s}) with |zeta^{hc} e^{-s}| = e^{-Re s} < 1. So 1 - w lies in the right half-plane, and the principal log is analytic on Re s > 0. It is also analytic at 0, since zeta^{hc} != 1.
- The EM path t·sigma, t in [0, 1], has Re = tL > 0 for t > 0.
- The nearest singular crossing is at t·y = 2pi m/5, where Re >= (2pi/5)(0.1/6) = 0.021. So |Phi''| is at most about 2e3 and integrable.
- K reaches up to 5.56e3 (log), so K/N0 = 6.8. That term is inside V (`K/n` in box_bound). The worst V of -12.0 already includes it.
- `certify` raises the number of t-pieces when a box gives a non-finite value. In the spot check, 13 tiny boxes at L = 0.1, |y| ≈ 4–6 were NaN at nt = 32 and finite at nt = 128, with the bound above the truth (audit_offpeak_spot_nt.log).

**(e) (F) maximum principle and rounding: sound. The float margins are now superseded by Arb.**
- **Chunk factorisation.** k = 5Tm + k', and V_m = q^{5Tm} has |V_m| = e^{-m lambda}. The drift r^{k'} equals e^{i k' (theta + i mu)}.
- **Tail.** 4T' <= 4(T-1) factors (T' = n - MT), each at most 1 + e^{-(L - lambda)}.
- **n p(L) lower bound.**
  - Midpoint (Hermite–Hadamard) plus convex telescoping gives the subtraction T(lambda + log5 - Q(L))/2.
  - The code subtracts T[lambda + (log5 - Q(L))/2], which is larger, so it is conservative.
  - Q > 0, so dropping the integral over [M lambda, L] is valid.
- **B(s).** B(s) is convex in s by Hadamard three-circles (log-convexity of max modulus in log|V|) and non-increasing in s. The cell bound d_i = min over the endpoints of (tangent of TQ minus chord of t) is the correct minimum of a difference of two linear functions.
- **Subharmonicity.** log|P(V, thetat)| is subharmonic in thetat, since P is entire in thetat. The sup over the compact circle |V| = e^{-s} of a continuous family of such functions is subharmonic. So the maximum principle on [th_a, th_b] x [0, mu_up] holds, and the code checks the bottom edge, the top edge and both sides.
- **Float evaluator.**
  - The cos lower bound c(1 - w^2/2) - |s| w (for c >= 0) and c (for c < 0), together with the 1e-12 margins and the treatment of x as convex in rho, is correct.
  - frexp/ldexp are exact. The ldexp clip at +1070 overflows to inf only when the product is genuinely below the threshold.
  - The crosscheck, however, tested only zero-width cells (dl = et = 0), so it never exercised the width margins. This is weak evidence, but it is now moot.
- **D3.** In the assembly, d_min is chosen with a float key on exact arb values. An ulp-level mis-ordering is possible. It is fixed in the Arb version, which uses exact arb comparisons, and has no effect on the margins.

**(f) Spot validation: 0 violations.**
- The check used exact arb full products at 312 random points, in audit_offpeak_spot.py and audit_offpeak_spot.log.
- n was 821, 1000, 3000 and three random n in [821, 2500].
- L included values in (0, 0.1), 0.1, 2.1 and [1.8, 2.1].
- theta was placed at |y| = 1.1+, 6, 20, 0.25n, dist = 0.0499/0.05, pi, 2pi/3 and pi/2, for all h = 0..4 and both halves [0, pi] and [pi, 2pi].
- Every point satisfied true V := log|S| - n p + E(n) <= the worst certified V of its piece. Every point with a computable piece box bound at that n (E, E0, N-a, N-b) was also below that bound.

| piece | points | max true V | certified worst V | min (box bound − truth) |
|---|---|---|---|---|
| (E) | 94 (+13 re-run at nt = 128) | -43.3 | -12.0 | 4e-4 (EM is nearly exact at L ≈ 2, n = 3000) |
| (E0) | 50 | -521 | -0.36 | 0.38 |
| (N) | 105 | -441 | -1.8 (h = 0) / -23.3 | 107 |
| (F) | 63 | -448 | -47.9 (as a margin) | n/a |

- The worst corner at n = 821 (audit_offpeak_focus821.log) is L = 2.1, |y| = 1.1+, where the true V is -43.2 against a certified -12.0.
- The other maxima at n = 821 are:
  - (E0): -494.
  - (N)/(F) handover at dist 0.0499: -447.
  - theta = pi: -1085.

## 2. Defects found

| # | where | defect | effect | status |
|---|---|---|---|---|
| D1 | offpeak_chunk.py pieces | (F) ends at float pi < pi; the gap (fl(pi), pi] has width 1.2e-16 and is uncovered | none numerically; a formal gap | fixed in the arb version |
| D2 | GAP_OFFPEAK §5 / offpeak_tb.py | (N) is certified only for h ∈ {0, 1} and y > 0; the reduction to them is unstated | none; the reduction is valid (see (a)) | needs one sentence in the doc |
| D3 | offpeak_chunk.assemble | d_min is selected by a float key | at most an ulp | fixed in the arb version |
| D4 | offpeak_crosscheck.py | the float checker was validated only on zero-width cells | evidence gap, not an error | superseded by the arb version |

## 3. Arb version of (F): offpeak_chunk_arb.py

- **Design.**
  - The certificate structure is the same as the float version.
  - Each cell (theta-interval, alpha-arc, mu-interval) is given by integer dyadic addresses. Its arb balls are built by `union` of exact endpoints, so the leaves tile the rectangle exactly.
  - Each factor is bounded by max over rho ∈ {rho_lo, rho_hi} of (1 - 2 rho c_lo + rho^2), where:
    - c_lo = (alpha + k' theta).cos().lower(), with arb.
    - rho_lo and rho_hi are the ends of the arb enclosure of e^{-s - k' mu}.
  - The arb product is compared with arb e^{2t}.
  - The float checker is used only as a search heuristic, to choose splits and target values. A cell counts only if arb passes, and a cell where arb fails is split further.
  - The assembly is arb, with exact comparisons, and the last rectangle ends at a float >= pi.
- **Reproduction.** 7 rectangles, including the two float-worst ones:

  | rectangle | arb margin | float margin |
  |---|---|---|
  | theta ∈ [3.1333, pi] | 47.91 | 47.91 |
  | theta ∈ [0.0499, 0.0599], T = 40 | 51.76 | 51.76 |
  | [0.10, 0.12] | 80.09 | 80.09 |
  | piece 2, rect 0 | 198.80 | 198.80 |
  | piece 2, rect 1 | 209.81 | 209.81 |
  | piece 6, rect 25 | 195.42 | 195.42 |
  | piece 6, rect 26 | 183.28 | 183.28 |

  - Every margin is identical to the float version, because the same targets pass.
  - Only 0–2% of the float-passing cells failed in arb and had to be split. This happened near dist 0.05 and near pi, and shows that arb's enclosures are slightly wider than the float version's.
  - Log: offpeak_chunk_arb_slice.log.
- **Speed.**
  - Per rectangle, the CPU time is 8–31 s. The float version takes about 3 s.
  - The whole run needs about 163 rects × ~18 s ≈ 50 CPU-min. That is about 6x the float version, not the 50–100x estimated before.
  - On this 2-core box it would take about 25–30 min wall time when the box is idle. At the time it was shared with three other jobs (load average 4), which roughly doubles the wall time.
  - That is above the ~20 min limit, so **the full run was not executed**.
- **Run commands.** Each piece, or any rectangle range within a piece, is independent:
  ```
  cd /home/claude/third-borwein
  python3 offpeak_chunk_arb.py list 821          # pieces 0-7: rect counts 6,6,6,6,6,53,53,27
  for p in 0 1 2 3 4 5 6 7; do
    nohup python3 offpeak_chunk_arb.py 821 $p -j 2 > arb_F_p$p.log 2>/dev/null
  done &                                        # ~25-30 min on an idle 2-core box
  # or split further across machines, e.g.: python3 offpeak_chunk_arb.py 821 5 0 27 ; ... 821 5 27 53
  python3 offpeak_chunk_arb.py summarize arb_F_p*.log   # must print "CERTIFIED for all n >= N0 (all rects OK)"
  ```
  The expected result is that every rect is OK, with a minimum margin of 47.91 at the rectangle ending at pi.

## 4. What remains

1. Run the full Arb (F) job above, then replace the float caveat in GAP_OFFPEAK §7.
2. Add the D2 sentence to GAP_OFFPEAK §5.
3. The combination with the near-peak certificate (|OFF| <= 0.0187 against slack 0.019 at n = 821) was not audited here.
