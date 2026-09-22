# Adversarial audit of the Section-2 citation (G_5 expansion, p = 5, delta = 1)

Date: 2026-09-18.  Scripts: scratchpad `chk.py` (python-flint arb ball arithmetic, precision set per i
BEFORE constructing B = 1/3; ratio reported as the ball's upper endpoint).

## 0. What the source actually contains

There is **no stated proposition-with-proof for G_p** in sz_joint.tex.  Section "Explicit expansions and
error constants for G_p and Q_12" says the dissection "is carried out as in Section [circle]" and lists
B, C, E_arc, E_rem, E_ng.  The only written proof is Proposition [explicit expansion] for Q_8
(hypothesis 2 < delta <= 4).  So the citation is "Q_8 proof, transplanted", and the audit below
re-derives each step for G_5 at delta = 1 rather than checking hypotheses of a written G_p statement.

## 1. Step-by-step verdicts

| step | content at p = 5, delta = 1 | verdict |
|---|---|---|
| Transformation formula | eta transformation (SZ Prop 7, hh' = -1). delta = 1 is an integer power, so there are NO branch issues; the certified-branch lemma (stated only for p = 7, 11) is **not needed** at delta = 1. Only the k = 5 phases enter, and those are confirmed by the amplitude match s(i)/P_5 -> a(i). The source's own inconsistency ("hh' = -1" vs "hh' = 1 in this display") only affects the correction-term phase, which is not used. | applies |
| Growing cusps = {h/k : 5 \| k} | d = gcd(5,k). d = 5: exponent pi(d^2-p)/(12pz) = pi/(3z) > 0 for every h, and G-hat -> 1. d = 1: exponent -pi(p-1)/(12pz) = -pi/(15z), decays. | applies |
| B, C | B = delta(p-1)/12 = 1/3; e^{2pi n z/k^2} e^{-(p-1)pi z/(12k^2)} gives C = 2n - 1/3. | applies |
| Amplitude modulus, phi(k) | prefactor (p/d)^{1/2} = 1 at d = 5, phase is a root of unity; one term per h coprime to k, so \|a_k\| <= phi(k). | applies |
| Correction term | grows iff delta > 24/(p-1) = 6; delta = 1, so it stays in the remainder and all [.] terms are dropped. Needed monotonicity of e^{pi B w} sum_j phi_j e^{-2 pi j w} in w needs B <= 2 (analogue of Q_8's delta <= 4): B = 1/3. | applies |
| K_p (non-growing) | on chords Re(1/z) >= 1: \|sqrt5\| * e^{-pi w/15} * f(e^{-2pi w/5}) * prod(1+t^m) <= sqrt5 e^{-pi/15} f(e^{-2pi/5}) f(e^{-2pi}) = 2.8549954. Uses (x;x) majorised by f(\|x\|), valid. | applies |
| exp bound | \|e^{pi C z/k^2}\| <= e^{2 pi C/(N+1)^2} <= e^{4 pi n/(N+1)^2} <= e; needs N + 1 >= sqrt(4 pi n). | applies |
| Truncation N | Proof needs only N >= sqrt(4 pi n) - 1 (source states the stronger N >= sqrt(4 pi n)+1; Q_8 version N >= sqrt(4 pi n + 2 pi delta)). K(i) = floor(sqrt(4 pi i)) + 2 > sqrt(4 pi i) + 1 always. Any admissible N works (constants are N-uniform: bounded by N/(N+1) < 1), so "+2" is merely a valid choice. | applies |
| E_arc | two arcs, length <= (pi/2) sqrt2 k/(N+1) each, integrand <= e^{1+pi/3}; sum_{5\|k<=N} phi(k)/k <= (4/5)(N/5); /(N+1) gives 4/25. = 5.506447 | applies |
| E_rem | chord length 2sqrt2 k/(N+1), integrand <= e^{1+pi/3}(f(t)f(t^5) - 1), t = e^{-2pi}; same count. = 0.006571 | applies |
| E_ng | sum over all k <= N, phi(k) <= k values, chord 2sqrt2 k/(N+1) / k^2, integrand <= e K_5: 2sqrt2 e K_5 = 21.950524 | applies |
| Assembly | E_con = 27.4635416 (independently recomputed, mpmath dps 30). The stated 27.4635 is a truncation BELOW the true 27.46354; difference 4e-5 is irrelevant given the numerical margin but the proof should use 27.4636. | applies with a check: round up |

## 2. Hypotheses that could silently fail at p = 5

- "p prime > 3"/"p >= 7": none used; nothing in the derivation uses p >= 7.
- delta threshold 24/(p-1): at p = 5 it is 6; delta = 1 is below, so no correction extraction. Consistent.
- Q_8 template's 2 < delta <= 4: its role (growth + remainder monotonicity) translates to 0 < B <= 2; B = 1/3 fine.
- Branch lemma omitted for p = 5: irrelevant at delta = 1 (integer power). The notes' k = 20 "branch +-1" saga is moot.
- Dedekind-sum/"hh' = +-1" convention: affects phases only, not moduli; k = 5 phases verified numerically.
- Residue-class-restricted lemmas: none used by the G_p route (Lemma "which cusps grow" is Q_8-specific).
- Small n: the proof has no large-n hypothesis (n >= 1 only); the bound is exact-inequality for every n >= 1.

No hypothesis fails or is unverified at p = 5, delta = 1.

## 3. Numerics (arb balls, rigorous per-i evaluation; max ball radius of ratio ~1e-60)

Coefficients cross-checked against a direct product for i <= 3000 (identical).

| range | mode | max ratio by class i mod 5 (0,1,2,3,4) | fail |
|---|---|---|---|
| 1..50 | exhaustive | 0.0196, 0.0244, 0.0296, ~1e-60, ~1e-60 | none |
| 1..300 | exhaustive | 0.238, 0.358, 0.538, ~0, ~0 | none |
| 1..2000 | exhaustive | 0.34548, 0.55899, 0.90445 (i=1987), ~0, ~0 | none |
| 1..20000 | exhaustive | 0.3454915, 0.5590170, **0.9045085 (i = 19987)**, ~0, ~0 | none |
| 20001..60000 | exhaustive | 0.3454915, 0.5590170, 0.9045085 (i = 42607), ~1e-87, ~1e-87 | none |
| 60001..200000 | every 7th i (20000 values, all classes mod 5 and mod 10) | 0.3454915, 0.5590170, 0.9045085, ~1e-109, ~1e-109 | none |

Largest ratio anywhere: 0.9045084971874737 = (5+sqrt5)/8 to 16 digits (limit approached from below).
Smallest failing i: none found in [1, 200000].

Note: the worst class is **i = 2 mod 5** (not 0 mod 5 as the task brief said); the limits are
((5-sqrt5)/2)/4, sqrt5/4, ((5+sqrt5)/2)/4 = 0.34549, 0.55902, 0.90451 -- the k = 10 amplitude over phi(10).
Classes 3, 4: s(i) = 0 and a(i) = 0, ratio is pure rounding (~1e-60).
Small i is the safest regime (ratio < 0.03 for i <= 50): E_con dominates there.

## 4. Verdict

The citation is sound at p = 5, delta = 1: every step of the transplanted Q_8 proof re-derives, all
constants reproduce, the truncation is admissible, and the bound holds with no failure (max ratio
0.9045085 < 1).  Caveats: (i) the source has no written G_p proof -- the argument above IS the proof and
should be written into PROOF.md rather than cited; (ii) use E_con = 27.4636 (round up).
