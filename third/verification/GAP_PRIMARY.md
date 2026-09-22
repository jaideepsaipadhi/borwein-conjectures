# GAP_PRIMARY: audit of the primary Laplace certificate (B_uniform5.py + uniform_lib5.py, driver tile_c.py)

The driver is `tile_c.py`. `grid5.py` only supplies `certify`; its own `run()` is a diagnostic. `tile_c.py`
bisects each v-cell to depth 5. For each cell it takes the first rung t of the LADDER that passes,
with **one** eps ball [0, a5/(t v_hi^2)].

Audit scripts: `audit_B_nan.py` (B_uniform5 with the NaN guard), `audit_B_mod.py` (the fixes of items
1b–1d), and `audit_runnan.py` (the driver). The rerun of tile_c with W0 = 0.02 reproduces the claim:
73/74/74 tiles, 0 unresolved, worst B 0.969/0.943/0.931 for classes 0/1/2. NOTES reports
141 tiles per class, which corresponds to a finer W0.

## FATAL: zone-2 NaN cells are silently counted as zero

B_uniform5.py lines 79–80:

    e1_=-1e6 if not (e1_==e1_ and e1_<1e6) else e1_

When a zone-2 cell's bound is NaN or at least 1e6, the code replaces it with exp(-1e6). The NaN
does occur, because Ghat(z) evaluates to NaN. `D_coef` divides by S = sum_{i<5r} q^i. With the
whole eps ball [0, emax] in one q ball, and complex z, the ball for S contains 0 once r is about 10
or more. At v = 2.5125 and t = 4.5, S is [+/-53] at r = 10, and D_coef is NaN at r = 30. From
t of roughly 4 to 7 up to t_C = 10v, **every** zone-2 cell is NaN and is dropped.

- Instrumented run over the class-0 tiles: 7168 zone-2 cells were dropped.
- With the guard added (a NaN makes the cell fail): **63/73 (class 0), 64/74 (class 1) and 64/74
  (class 2) of the tiles no longer certify.** They fail with 'zone2 NaN' at t between about 4.3
  and 7.
- Only the roughly 10 tiles nearest v = 10 survive.
- So Z2 is not an upper bound over most of zone 2. The 983-tile claim does not hold as run.

**Repairable in principle.** A thin eps ball at the top, [0.98 emax, emax] on [2.5, 2.525] at t =
11.4, certifies with B = 0.631. But thin balls near eps = 0 fail the zone-2 slope test, because V
and eps enter A = Re(G0z - G0v) as independent balls. Halving the eps ball 5 times made things
worse: 72/73 tiles unresolved. A working repair needs:
- fine eps splitting, with ratio about 1.02 near the top;
- bisection in v;
- a slope test that handles the dependency, for example by evaluating G0z - G0v as one series in v.

I have not done this; it will take many hours of CPU time.

## Item 1: normalisation. No defect of the B_ball kind, but three other errors in Z_R, Z3 and zone 1

**(1a) Main term.** Rh = N[a5/v + Re Tt_1 + mu v] + log|Re Ghat|. Here Tt_1 is the **h = 1 twisted**
log, so the main term is |E_n(z5 e^{-W})| times |Ghat|, which is the *small* quantity.

Z3 and Z_R use real-axis majorants (4 E_n(e^{-W}) and Rankin with E_n(e^{-w'}) via PsiU) minus
G0v = a5/v + Re Tt_1. The gap e^{Dz} is therefore already included. Measured Dz/N = PsiU/v - Re Tt_1:

| v | Dz/N | Dz at Nb |
|---|---|---|
| 4.04 | 0.0044 | 2.2 (Nb = 506) |
| 6.06 | 0.0004 | 0.22 |
| 9.9 | 1e-5 | ~0 |
| 2.61 (class 3) | 0.029 | 15.1 |

Z2 uses the twisted |Ghat(z)/Ghat(W)| and G0 on both sides, so it is consistent.

**(1b) Z_R carries an unexplained factor e^{v'-v}.** In grp(), line 96:
`val = lead + Nb*br + vp - V - log|Re Ghat|`.

Neither the Rankin bound (PROOF §6) nor NOTES §4a(ii) contains an e^{v'-v} factor. For the
optimal v' < v it lowers Z_R by up to e^{-(v-v')}. Removing it at [4.10, 4.18], t = 4 (class 0)
raises Z_R from 3.3 to 10.6, where 3.3 already includes the ×2 of (1c).

The lead 3φ(k) has slack: the measured sum_i P_k(i) e^{-iw} / e^{a_k/w} is at most 1.0. That
slack does not offset this factor.

**(1c) A factor 2 is missing in Z3 and Z_R.** The code uses π N^{3/2} sqrt(G2); it should be
2π N^{3/2} sqrt(G2). The units are e^{Rh} σ_y/(2π) with Lint ≈ sqrt(2π). A zone-3 length of up to
2π, and the 1/(2π) of the coefficient integral, both give a factor 2π.

**(1d) Z_E, the representation error (-1 + eps(w)), is missing from zone 1 and zone 3.** It
appears only in zone 2, as the 3 e^{e2} term. In zone 3 it adds at most D0 Z3, and in zone 1
about D0 e^{-tv}. Both are negligible at N ≥ 1456. D0 = 1.126 is certified only for W ≤ 0.02, but
tiles with t = 4 and v = 2.5 have W = 0.026. The code's constant 3 probably covers this, but it
is not certified.

**Effect of 1b + 1c + overlap + zone-1 margin (audit_B_mod.py) at each tile's own t:**
17/73, 17/74 and 15/74 tiles fail (classes 0/1/2), with worst B = 4.7, 4.9 and 4.5 at v ≈ 3.9–4.4
and t = 4–4.5. Z_R dominates there.

At t_eff = max(t, a5·1456/v_hi^2), the t that PROOF §9 actually needs, Z_R is negligible. All
tiles that pass have B ≤ 0.43, but 41–42 tiles per class fail the zone-2 slope test at the larger
t, because of ball dependency. So "B non-increasing in t" does **not** hold for this code as
implemented.

## Item 2: bounds proved only for real q, used off the axis

- **Tt_1 tail** (6 x^{R+1}/(R+1)): valid off the axis, since |1 - e^{-y}| ≥ 1 - e^{-Re y}. Measured
  max |coef| for r > R is 0.014.
- **D_h tail** (4·tail): the measured sup of r|D_coef| on the walk, over N from 151 to 1e5, is 9.0,
  against the implicit bound of about 24/3 per h. That is marginal but holds. The tail is added to
  Ghat after exponentiation, which is valid only if |e^{D_h}| = O(1).
- **The 1e-25 floor and early break:** the true tail after the break is below 1e-26 in every case
  checked, so it does not matter numerically. It is still not a proved enclosure.
- **PsiU:** the real-axis bound is used only against real-axis majorants, which is legitimate. But
  PsiU with Nlo = 1000 is **false for N below about 500**. Direct log E_n vs bound: 4.119 > 4.077
  (N = 151, v = 2.5); 8.143 > 8.127 (N = 301). The tiles claim N ≥ Nb = 151. This is harmless for
  N ≥ 1456, where the bound holds and the Nb·br monotonicity argument is valid.
- **Zone-2 slope test:** sound as a test. It is also the reason larger t fails (dependency, above).
- **Zone 1 has no Taylor remainder.** G4 is computed and never used. Lint assumes |integrand| equals
  e^{-s^2/2} exactly, with only the cubic phase, and Z1 assumes a Gaussian tail on [c, s_s].
  Measured true zone-1 integral against Lint/2:
  - 1.2341 vs 1.2406 (class 0, v = 2.55, N = 1456);
  - 1.2268 vs 1.2323 (class 3, v = 2.61, N = 522).

  Lint is **not** a lower bound; it is about 0.5% too high. The excess of Re phi over -s^2/2
  reaches 0.63 at N = Nb.
- **Overlap:** when s_s = t_s sqrt(N G2) < c (s_s = 2 at N = Nb in most tiles), Lint integrates up to
  c = 3, over a region that zone 2 also covers. The piece over s_s < |s| < c has to be charged a
  second time. It is not.

## Item 3: coverage

- The eps ball [0, a5/(t v_hi^2)] covers every N ≥ t v_hi^2/a5 in the cell, so no monotonicity is
  needed for the claim at t ≥ c(v).
- All tiles have c(v) ≤ max(4, a5·1456/v_hi^2), so the §9 requirements for v in [2.5, 10] are met
  on paper for classes 0/1/2.
- Zone 3 is t ≥ 10v up to y = π, bounded by Z3, and its β < 0 check passes.
- The signs of Re Ghat(W) are correct: +3.63 for class 0, and -2.21, -1.39, -0.0186, -0.0183 for
  classes 1–4.
- Not covered here: classes 0/1/2 for v in [10, 20]. The §7 tail cell and the derived tiling concern
  classes 3 and 4.

## Verdict

**UNSOUND as run.**
- The main term is correctly normalised: defect (1) of the derived tiling is absent.
- The fatal defect is the NaN-to-zero substitution in zone 2. It invalidates about 86% of the tiles
  in every class.
- Secondary defects: the e^{v'-v} factor in Z_R, the missing factor 2 in Z3 and Z_R, the absent
  zone-1 Taylor remainder and zone overlap, Z_E missing outside zone 2, and PsiU false for N < 500.

**Classes 0, 1 and 2 on v in [2.5, 10] cannot rely on this certificate.** It might be repairable
at t_eff with fine eps splitting, a dependency-aware slope test, the factor fixes and a zone-1
remainder. At t_eff Z_R is negligible, and the thin balls that did pass give B ≤ 0.43–0.63. That
repair has not been done.
