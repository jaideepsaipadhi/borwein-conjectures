# GAP_SMALL: zone 2 in the tail cell, and the missing twisted-E_n factor (classes 3, 4, Laplace)

Scripts: `small_tail.py` (+ `.log`), `small_probe.py`, `small_fast.py` (+ `.log`), `small_window.py`,
`small_tailcorner.py` (+ `.log`), `small_Glan.py` (+ `.log`).

**Summary.** Gap 1 is closed *within the existing framework* by splitting at T0 = 5v. Gap 2 is real, and larger
than stated. The missing factor affects the main term itself, not only Z_2. As a result the whole kernel-centred
Laplace certificate for classes 3 and 4 (`mono_cert.py`, `zone2_cert.py`, `tail3_cell.py`) is **unsound**: its
central claim |L(m) − Lc| ≤ B·Lc is numerically false at points it covers. The sign of L(m) is still correct at
every point tested, but that is no longer proved.

## 1. Gap 1: zone 2 inside the v ≥ 20 tail cell

**Fixed T0 = 50 does not work.** On the tail boundary W = a5/(Tv), the far term's exponent is about
v(1 − T·s(v)), with s = 2500/(2500 + v²) → 0. For T = 2 this exponent is −14.5 at v = 20, −8.8 at v = 40,
**+22.7 at v = 70** and +176 at v = 200. The far part is therefore not dominated for v ≳ 50√T, so T0 must
scale with v.

**Fix: T0 = τv with τ = 5,** i.e. t1 = 5W.

- **Near part.** |arg q| ≤ 5W ≤ 5W0 = 0.03290 < 50/1456 = 0.03434 (t ≥ 2). Also Re z = W ≤ 0.00658 < 20/1456.
  So every near-part q lies in `mono_Gl.py`'s certified box, uniformly in v ≥ 20 and N. The Cauchy constant
  of `zone2_M.py` holds on the sector |Im y| ≤ 10 Re y, and 5 ≤ 10. Hence C(v) is valid on the near part.
- **Far part.** Take s = 25/26. The far term Zfar = (4/(K0x))·√(2π)·√(2a5)·W^{−3/2}·exp(−(a5·s − Li2 x)/W + 5 Li1 x)
  has the e^{Dz}Z3 form with κ = a5·s − Li2(x). The GAP_TAIL joint argument carries over:
  - W0 < κ/p holds;
  - the boundary log-derivative in v is ≤ −0.848 for t ≥ 2 and ≤ −2.771 for t ≥ 4.

  So the supremum is at (20, W0).

| t-range | Z2near | Z2far | **B_tail** (was) |
|---|---|---|---|
| t ≥ 2 | 0.00201 | 1.3e-4 | **0.04548** (0.04535) |
| t ≥ 4 | 0.00120 | 7.3e-21 | **0.02237** (0.02237) |

These values are in arb and all side conditions are OK (`small_tail.log`). **This is valid only as a repair of
the tail cell's zone-2 bookkeeping. Section 2 invalidates the framework it sits in.**

## 2. Gap 2: what the zone-2 integrand is, and what is missing

**Exact structure.** Let P̃(q) = (E_n(ζq)·E_n(ζ⁴q))^{1/2}. It is analytic, and on the real axis P̃ = |E_n(ζ e^{−W})|.
The G of the sign lemma continues analytically as G = F_m/(2P̃). So along w = W + it,

  |F_m(e^{−w})| / |F_m(e^{−W})| = R_P · R_G,  where R_P = |P̃(q)|/P̃(r) and R_G = |G(q)|/|G(r)|.

C(v) bounds only R_G (a series in x with coefficients G_l, divided by K0·x/2). **R_P appears nowhere in
`mono_cert.py`, `zone2_cert.py` or `tail3_cell.py`.** It is not assumed ≤ 1 anywhere; it is simply dropped.
The far part of `zone2_cert.py` is fine, because |F_m| ≤ 4E_n(e^{−W}) together with the e^{Dz} factor
bounds R_P·R_G jointly.

**R_P is not ≤ 1** (`small_probe.py`, 40 digits). Measured log R_P:

| (n, v) | u = 2.9 | u = 11.6 | at the far point shown |
|---|---|---|---|
| (200, 2.5) | 1.73 | 11.3 | 6.6 at N·θ = 50 |
| (300, 4) | 0.27 | 2.3 | — |
| (600, 10) | 4e-4 | — | — |

log R_G stays within ±0.11. Here u = t/σ with σ = √(W³/(2a5)). The kernel's decay dominates R_P in zone 2
(log K = −41 against +11.3), so for Z2 alone the factor could be absorbed.

**The same factor breaks the main term, which is the real problem.** mono_cert takes W to be the kernel-only
saddle (m = a5/W²). It takes the main term to be Lc = e^{a5/W+(m−1/6)W}·|F_m(e^{−W})|·σ/√(2π), and assumes a
window error ε0 computed from the kernel alone (f1 = 0, f3 = 6a5/W⁴). But log F_m = log P̃ + log G + log 2
contributes three things:

1. A linear phase in u of size σ·(log P̃)′.
2. A linear phase σN from the thin factor x = e^{−Nw}. This is the m → m − N shift of the saddle.
3. A curvature change δ = −σ²(log P̃)″, up to 0.44.

None of these is in ε0. Two independent float computations agree to about 0.3%:

- the exact line integral I over |u| ≤ 14 at W;
- the Laplace value at the true saddle W* of a5/w + (m−1/6)w + log(−F_m).

| n | v | t | **I/Lc** | Laplace(W*)/Lc | β_P | δ | row's certified B |
|---|---|---|---|---|---|---|---|
| 291 | 2.5 | 61.3 | **0.0774** | 0.0775 | 1.99 | 0.44 | 0.0208 |
| 291 | 3.0 | 42.6 | 0.511 | 0.512 | 1.24 | 0.33 | 0.0240 |
| 291 | 3.6 | 29.6 | 0.955 | 0.957 | 0.71 | 0.24 | 0.0304 |
| 291 | 4.5 | 18.9 | 1.075 | 1.078 | 0.31 | 0.13 | 0.0427 |
| 291 | 6.0 | 10.6 | 0.969 | 0.972 | 0.08 | 0.05 | 0.0435 |
| 291 | 8.5 | 5.3 | 0.752 | 0.755 | 0.007 | 0.007 | 0.0555 |
| 291 | 10 | 3.83 | **0.587** | 0.590 | 0.002 | 0.002 | 0.0639 |
| 2000 | 2.5 | 421 | **~0 (e^{−13})** | same | 5.20 | 0.44 | 0.0208 |
| 2000 | 4 | 165 | 0.474 | 0.474 | 1.27 | 0.18 | 0.0304 |
| 20000 | 6 | 731 | 0.865 | 0.865 | 0.62 | 0.05 | 0.0435 |
| 1216 | 20 | 4.0 | **0.295** | — | 0 | 0 | tail 0.0224 |
| 608 | 20 | 2.0 | 0.054 | — | 0 | 0 | tail 0.0454 |

The last two rows are from `small_tailcorner.log`, at 30 digits. The rest are from `small_fast.log`.

The certificate asserts |I/Lc − 1| ≤ B. That is **false** in every row with |I/Lc − 1| > B, which is most of
them. The failures have two sources:

- **The thin-factor shift, at small t.** W*/W = √(t/(t−1)). The loss is exp(−(a5/W)(1 − √(1−1/t))²), which is
  0.24 at (v, t) = (20, 4).
- **The P̃ shift, at small v.** At fixed v, β_P grows like √t, so I/Lc → 0 as N → ∞. The certificate cannot
  cover t → ∞ at fixed v ≲ 4 in any form.

I/Lc stays positive and real in every row, so the sign of L(m) looks right, but B < 1 no longer implies it.

**What a fix requires (not done here).** Re-centre the Laplace analysis on the true saddle of
Rh(w) = a5/w + (m − 1/6)w + log P̃(w) − N·w, or of the full a5/w + (m − 1/6)w + log(−F_m). This is the scaled
R = N[G(ζ) + μζ] of NOTES §4, which the mono_cert chain simplified away. Concretely:

- main term e^{Rh(W*)}·|2G(W*)|/√(2π Rh″);
- ε0 rebuilt with f3 and f4 of Rh, and with G's slow variation as f1;
- Z2's integrand exp(Re Rh(w) − Rh(W*)) times C(v), which needs a proved bound on log R_P;
- a re-derived map from (n, m) to (v, W*) for the monotonicity and coverage arguments.

That is a rebuild of GAP_MONO / GAP_ZONE2 / GAP_TAIL, not a small patch, so no per-row B is reported. Any B
computed by adding R_P to Z2 alone would be meaningless while the main term is off by up to e^{13}.

**Smaller caveat.** `Gl_box` (used by `mono_Gl.py`) returns Re(W1ζ^{−m}c_l(q)) at complex q. The analytic
coefficient is (f(q) + conj f(conj q))/2. At 12 sample points on the box edge the two differ by ≤ 0.009 at
l = 1 and ≤ 0.04 at l = 2, and the analytic values stay under Gball (max 0.509, 0.345, …), per `small_Glan.log`.
This is not a certificate, and the enclosure should be redone with the analytic form.

## Verdict

- **Gap 1: closed within the old framework.** T0 = 5v, near part inside the certified G_l box and the Cauchy
  sector, far term dominated uniformly for v ≥ 20. B_tail ≤ 0.0455 for t ≥ 2 and ≤ 0.0224 for t ≥ 4. A fixed
  T0 = 50 genuinely fails for v ≳ 70.
- **Gap 2: a factor is missing, and the defect is not small.**
  - Z2's near part omits R_P = |P̃(q)|/P̃(r), which reaches e^{11}.
  - More seriously, the main term and ε0 ignore log F_m altogether. The kernel-only main term Lc is off by
    factors of 0.59 (v = 10, N = 1456), 0.077 (v = 2.5, N = 1456), 0.054 (tail, t = 2), and e^{−13} or worse
    as N grows at small v.
  - The Laplace certificate for classes 3 and 4 on v ∈ [2.5, 20] and on the tail cell (`mono_cert`,
    `zone2_cert`, `tail3_cell`) is therefore **not valid**. PROOF.md §13/§14 "closed" entries for monotonicity,
    zone 2 and the tail cell should be downgraded to **open pending a re-centred (true-saddle) rebuild**.
- Classes 0–2 (`cls012_lap.py`) were not examined here.
