# GAP_TAIL: the v >= 20 tail cell (classes 3, 4), rebuilt

Scripts: `tail3_K0.py` (thin factor for v >= 20), `tail3_cell.py` (the cell), `tail3_cell.log`.

## Result

| K_0 used | t-range | B_tail (all v >= 20) | Z1 | eps_0 | Z2 | e^Dz Z3 | e^Dz Z_E | e^Dz Z_R |
|---|---|---|---|---|---|---|---|---|
| **0.999 (proved here)** | t >= 2 | **0.04535** | 4.7e-4 | 0.0428 | 0.0020 | 2.8e-5 | 3.2e-5 | 4.5e-6 |
| 0.999 | t >= 4 (what §9 needs) | 0.02237 | 4.7e-4 | 0.0207 | 0.0012 | 3e-22 | 4e-22 | 5e-23 |
| 0.25 (for comparison) | t >= 2 | 0.05156 | | | 0.0080 | 1.1e-4 | 1.3e-4 | 1.8e-5 |
| 0.0265 (for comparison) | t >= 2 | 0.12137 | | | 0.0757 | 1.1e-3 | 1.2e-3 | 1.7e-4 |

All side conditions hold (c = 3.5 is best, min over c in {2, 2.5, 3, 3.5}). The old values of 0.0662 and 0.2981 are
superseded. They were larger mainly because of the coarser zone-2 sum. The corrected main term, e^{Dz} and Z_E
together move B by less than 1e-4 here.

**Coverage correction.** "t >= 2, i.e. all N >= 1456" is false. t = a5·N/v², and t_cov = a5·1456/v² drops below 2
for v > 24.1. It does not need to hold: by §9, a Laplace point has t >= 4 (t < 4 is the band). So the cell must
cover {v >= 20, t >= 4}, and it covers {v >= 20, t >= 2}, which contains that set.

## K_0 for v >= 20 [proved; certified in ball arithmetic]

Ghat_m = Σ_h c_h e^{D_h}, with c_1 = W1 ζ^{-m}, c_2 = ζ^{-2m}, c_3 = ζ^{-3m}, c_4 = conj(W1) ζ^{-4m}, and |c_h| = 1.
Also D_1 = 0 and D_h = Σ_{r>=1} x^r d_r(q,h), where

  d_r = Σ_{c=1..4} a_c s^{c-1} / (r(1 − s^5)),   a_c = ζ^{hcr} − ζ^{cr},   s = q^r

(the `D_coef` of B_uniform5.py).

**The bound |d_r| <= 8/r.** Σ a_c = 0, so Σ a_c s^{c-1} = −(1−s) Σ a_c (1+…+s^{c-2}). Dividing by
1 − s^5 = (1−s)(1+…+s^4), each ratio (1+…+s^{c-2})/(1+…+s^4) lies in [0,1] for s in [0,1], and |a_c| <= 2.

**The expansion.** In the thin classes Σ_h c_h = 0 (checked to 1e-33). Hence

  Ghat = x·f1(q) + Rem,   f1 = Σ_{h>=2} c_h d_1(q,h),
  |Rem| <= 3·[8(Li1(x) − x) + 32·Li1(x)²·e^{8 Li1(x)}].

Rem/x is increasing in x, and it is <= 2.3e-7 at x = e^{-20}.

**The first-order coefficient.** f1 is a closed-form rational function of q, with f1(1) = −1 exactly in both
classes. So |Re Ghat|/e^{-v} → 1, as observed numerically. On the Laplace part of the tail (t >= 2, v >= 20), W <= a5/40,
so q = e^{-W} >= 0.99344. Ball tiling on q in [0.99344, 1] gives |Re f1| >= 0.99994 (class 3) and >= 0.99995 (class 4).
Therefore

  **|Re Ghat_m(W)| >= 0.9999·e^{-v} for all v >= 20, all N (every ε), q >= e^{-a5/40}, classes 3 and 4.**

K_0 = 0.999 is used. It is uniform in N, and the sign of Re Ghat is fixed (negative).

Caveat: on all of q in [0,1] class 4 has f1(0) = 0, the known small-q weakness. The bound therefore depends on the
Laplace restriction q >= e^{-a5/40}. That restriction does hold in this cell.

## Terms and the joint (v, W) argument

Region: v >= 20 and W <= a5/(T0 v), with T0 = 2. Let W0 = a5/40 = 0.006580 and x0 = e^{-20}. Every mono_cert.py
correction is in place: the lower bound log E_n(e^{-W}) >= (4/(5W))Li2(x) in the main term, e^{Dz} on Z3, Z_E and Z_R
(Dz <= Li2(x)/W + 5Li1(x)), Z_E = D0·e^{Dz}Z3 with D0 = 1.12601 (valid for W <= 0.02), no k-cap in Z_R together with
1[k <= K] <= (K/k)³, and interval upper sums in eps_0. The code reuses the functions of mono_cert.py unchanged, with
K_0 overridden.

- **Z1** is constant.
- **eps_0** depends on W only and increases in W (GAP_MONO T1'). Evaluated at W0.
- **Z2** is C(v) times an integral in W only, which increases in W, plus the endpoint term (GAP_MONO T2). C(v) is a
  sum of nonnegative multiples of x^{l-1} and x^{L}, so it decreases in v. Evaluated at v = 20, W0.
- **e^{Dz}Z3 and e^{Dz}Z_E** have the form (const/(K0 x))·W^{-3/2}·e^{-κ(v)/W + 5Li1(x)}, with κ = a5 − Li2(x).
  - At fixed v they increase in W for W < κ/p = 0.1755 (> W0), so the maximum is on the boundary W = a5/(T0 v).
  - On that boundary the exponent is v − (T0/a5)·v·κ(v) + p·log v + 5Li1(x) + const. Its v-derivative is
    1 + p/v − (T0/a5)(κ + vκ') + 5Li1'(x) <= 1 + p/v − (T0/a5)κ(v), because κ' = Li1(x) >= 0 and Li1(e^{-v})
    decreases.
  - With κ(v) >= κ(20) this is <= −0.925 < 0. So the supremum is at (20, W0).
- **e^{Dz}Z_R** at a fixed t' = 0.5 is a finite sum of const·P(v,W)·W^{-p}·e^{-κ_i(v)/W}, with
  κ_i = κ0 − 25a5/(k²t') and p in {3/2, 9/2} (GAP_MONO T5).
  - Each κ_i is non-decreasing in v, since Li2(e^{-t'v}) and Li2(x) decrease.
  - Every v-dependence in P decreases along the boundary (y, x and W = a5/(T0 v) all decrease), except the e^{v}
    from 1/(K0 x).
  - Conditions checked at v = 20 for every term: W0 < κ_i/p (with margin 13.2), and 1 + p/20 − (T0/a5)κ_i(20) <= −0.924.
  - So the supremum is at (20, W0).

## Remaining assumptions (inherited, not introduced here)

1. Zone 2's C(v) uses G_l enclosures certified only for |arg q| <= 50/N, while it is applied over the whole period
   (PROOF §13 open item 2 / GAP_MONO finding 3). mono_Gl's box has Re z <= 20/1456, which contains every tail-cell
   W <= a5/40. So the tail inherits exactly the same open item, and nothing new.
2. Z2 carries no e^{Dz} factor, as in mono_cert. It is a ratio of exact F_m values, not a positive-coefficient majorant.
   GAP_REPR caveat 2 still lists this as unchecked.
3. The inputs taken as given elsewhere: ECON, the Rankin/R-term lemma, D0, and §3.

## Verdict

**Closed, modulo the inherited items above.** The corrected tail cell gives B <= 0.0454 on {v >= 20, t >= 2}, which
contains every Laplace point with v >= 20 (those have t >= 4, where B <= 0.0224). It uses K_0 = 0.999, proved for
v >= 20 by the first-order expansion. The K_0 = 0.25 extrapolation and the unbounded cos(Im T_1) factor are no longer
needed. Even with the old fallback K_0 = 0.0265 the cell gives B <= 0.121.
