# GAP_ZONE2: the arg-q range of zone 2 (classes 3, 4, Laplace, v in [2.5, 20])

Scripts: `zone2_M.py` (Cauchy constant on the sector), `zone2_cert.py` (+ `zone2_cert.log`, the re-run certificate).

## 1. Range actually used: the whole half-period. It was NOT covered.

In `mono_cert.py` (and `B_ball`), Z2 integrates u = t/sig from c to pi/sig, i.e. t = Im w from c*sig to **pi**.
Since q = e^{-w} (folded identity, t in [-pi, pi]), arg q = t (= Im zeta / N with zeta = N w). So Z2 applies C(v)
for |arg q| up to **pi**, in every row and for every t >= t_eff. G_l were enclosed only on |arg q| <= 50/N,
i.e. at most 50 W_hi / v_lo <= 0.035 rad (e.g. 0.0258 in row [16,20]). The range used is 1-2 orders larger and
**does pass through the fifth roots of unity** arg q = 2pi/5, 4pi/5, where 1 - q^{5r} -> 0 as |q| -> 1 and
the x-series coefficients blow up (G_1 ~ 25-51). The NOTES claim "zone 2 never reaches them (t_C ~ 10v)" referred
to an older zone-2 extent |Im zeta| <= 10v; it does not describe mono_cert/B_ball. The gap was real.

## 2. Fix: split zone 2 at t1 = T0/N, T0 = min(50, 10 v_lo) per cell

**Near part, |t| <= t1.** C(v) is valid there:
- G_l: `mono_Gl.py` box z in [0, 20/1456] x [-50/1456, 50/1456] contains every q = e^{-W-it}, W = v/N, |t| <= T0/N <= 50/N,
  N >= 1456, v <= 20. Max|G_l| <= Gball [certified, existing].
- Cauchy constant (new, `zone2_M.py`): with y = r z, |Im y| <= (T0/v) Re y <= 10 Re y uniformly in N, r. r a_r and
  r b_r are functions A_k(y), B_k(y), k = r mod 5, analytic on the closed sector (removable at 0, poles only at
  2 pi i j/5). Maximum principle (circle |y| = 0.1, the two rays, Re y = 4; crude bound beyond) in arb:
  sup r|a_r| <= 4.535, sup r|b_r| <= 9.003, hence M = (1-rho)^{-A} + (1-rho)^{-B} <= **10.21 < 25** (rho = 0.2).
  (This replaces the float-sampled `Mfix2.py`, and uses the sum formula rather than Mfix2's A/2 + B exponent.)
- The existing integral over u in [c, pi/sig] is kept; it over-estimates the integral over [c, t1/sig].

**Far part, t1 <= |t| <= pi.** No C(v). Use |F_m| <= 4 E_n(e^{-W}) (non-negative coefficients, valid on the
whole line), E_n(e^{-W}) / |E_n(z5 e^{-W})| <= e^{Dz} <= exp(Li2(x)/W + 5 Li1(x)), |Re Ghat(W)| >= K0 x, and the
kernel ratio exp(-a5 t^2/(W(W^2+t^2))) <= exp(-(a5/W) s), s = tau^2/(1+tau^2), tau = t1/W = T0/v >= T0/v_hi:

    Z2far <= (4/(K0 x)) sqrt(2pi) sqrt(2a5) W^{-3/2} exp(-(a5 s - Li2(x))/W + 5 Li1(x)).

Same shape as the Z3 majorant with a5 -> a5 s; non-decreasing in W for W < 2(a5 s - Li2 x)/3, checked per cell
(all OK), so the GAP_MONO evaluation at W_hi remains valid.

## 3. Re-run (`zone2_cert.py`, arb, 16 v-cells per row, min over c in {2,2.5,3,3.5})

| v-row | t_eff | Z2far | **sup_{t>=t_eff} B** | before (mono_cert) |
|---|---|---|---|---|
| [2.5, 3.0] | 42.57 | 4.3e-25 | **0.0208** | 0.0208 |
| [3.0, 3.6] | 29.56 | 8.7e-25 | **0.0240** | 0.0240 |
| [3.6, 4.5] | 18.92 | 3.5e-20 | **0.0304** | 0.0304 |
| [4.5, 6.0] | 10.64 | 1.5e-13 | **0.0427** | 0.0427 |
| [6.0, 7.0] | 7.82 | 3.9e-13 | **0.0435** | 0.0435 |
| [7.0, 8.5] | 5.30 | 1.1e-08 | **0.0555** | 0.0555 |
| [8.5, 10] | 3.83 | 4.1e-06 | **0.0639** | 0.0639 |
| [10, 13] | 3.00 | 2.9e-04 | **0.0705** | 0.0702 |
| [13, 16] | 2.25 | 2.4e-02 | **0.1046** | 0.0804 |
| [16, 20] | 2.00 | 1.2e-01 | **0.1980** | 0.0781 |

All side conditions OK. Worst row [16, 20]: B <= 0.198 (far part dominated by e^{Dz}/x there, tau = 2.5).

## Verdict

The gap was genuine (zone 2 used C(v) up to arg q = pi, through the fifth roots of unity, where it is not
enclosed and is false-in-spirit). With the split it is **closed**: C(v) is now used only on |arg q| <= min(50,10v_lo)/N,
where both G_l and the Cauchy constant are enclosed uniformly in N >= 1456, and the rest of the period is
bounded by the far-field argument. sup_{t>=t_eff} B <= 0.198 in every row of v in [2.5, 20].
Not addressed here (inherited): the structure of C(v) itself (that |Ghat(w)|/|Ghat(W)| and the twisted-E_n
ratio are what Z2 needs), K0, and the v >= 20 tail cell, which uses the same zone-2 bound and must get the
same split (at v >= 20, tau = 50/v shrinks, so the far term there needs its own check).
