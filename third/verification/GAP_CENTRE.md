# GAP_CENTRE: rebuilding the centre certificate (mu >= 0.105), status

*18 Sept 2026. Nothing below upgrades any PROOF.md label to [proved] or [certified].*

## Result: the centre is still OPEN for all five classes
No sound certificate was produced in this pass. The old centre runs stay withdrawn.

## What was established (analysis; numerics only where stated)
1. **Exact scaled form near a peak.** With theta = 2 pi h/5 + x, y = 5 n x, sigma = L - i y:
   log S_n = sum_c sum_{t<n} Phi_c((5t+c)u), u = sigma/(5n), Phi_c(s) = log(1 - zeta^{hc} e^{-s}).
   The continuum term is n p(sigma), p(sigma) = (1/sigma) int_0^sigma Q, Q(s) = log((1-e^{-5s})/(1-e^{-s})).
   **p is analytic for |sigma| < 2 pi/5.** So the near-peak Laplace analysis has no singularity as L -> 0 in
   the variable y: item C.2.7 (mu > 0.483) is an artefact of parametrising by theta, not a real obstacle,
   *provided* the window is |y| <= Y1 < 2 pi/5. In that window Phi_c is analytic along the ray (distance
   >= 2pi/5 - Y1 from its singularities), so the finite-sum EM lemma (below) has a remainder O(|u|) = O(1/n)
   that is uniform down to L = 0.
2. **Finite-range EM lemma (proved, corrects C.2.4):**
   sum_{t=0}^{n-1} f(t+c/5) = int_0^n f + B1(c/5)(f(n)-f(0)) + (B2(c/5)/2)(f'(n)-f'(0)) - int_0^n (B~2(y-c/5)/2) f'',
   so |R| <= (11/300)(|f'(0)|+|f'(n)|) + (1/12) int |f''|. The old 5/12·TV constant is invalid; 1/2·TV (first order) or this form is correct.
3. **Structure of |S_n| away from the peaks** (heuristic, from the same continuum at other rationals):
   near theta/2pi = a/(5b) (5 not dividing the numerator part) the exponent is (n/b) p(bL) <= n p(L)/2 for b >= 2;
   near rationals with denominator prime to 5 (including theta = 0) there are valleys.
4. **Numerics (n = 291, 600; `centre2_landscape.py`, float, NOT a certificate):**
   max over |y| > Y of [log|S(re^{i th})| - log|S(r zeta)|]/n:

   | L | Y=0.5 | Y=1 | Y=2 | Y=5 |
   |---|---|---|---|---|
   | 1e-6 | -0.086 | -0.40 | -0.80 | -0.80 |
   | 0.3 | -0.080 | -0.34 | -0.78 | -0.78 |
   | 1.0 | -0.045 | -0.17 | -0.49 | -0.61 |
   | 2.08 | -0.016 | -0.062 | -0.20 | -0.41 |

   The global maximum off the main peaks is the secondary peak at theta = 2pi/10 (and 3pi/5), excess -0.80n at small L
   (exactly half of log 5, matching item 3). So the true off-peak ratio is <= e^{-0.06 n} for |y| >= 1 at every L <= 2.08,
   i.e. <= e^{-18} at n = 291: ample room, the difficulty is only in proving it uniformly in n.

## Why (b), the away-from-peak bound, was not closed
Every uniform route tried fails at a specific point:
- **Triangle inequality on E_n** (|S_n| <= |G_5(q)| E_n(r)): loses exp((n/L)(5Li2(xi) - Li2(xi^5)/5)); works only for xi <~ 0.2 (L >~ 1.6).
- **Modulus bound on the j-series** (1 - q^{5nj} bounded by 1 + xi^j): loses the peak cancellation as xi -> 1.
- **Fixed-length chunk bound + maximum principle** (log|S| <= sum over chunks of max_{|U|=rho_i} log|prod_{k'<=5T}(1-U q^{k'})|):
  uniform in n and correct for theta with 5theta near a/b, b >= 2 (per-t bound 4 log(1+rho^b)/b, e.g. b = 2 gives
  0.388 n vs peak 0.571 n at L = 2.08), but within distance ~1/T of a 5th root it reproduces log 5 and sees no decay;
  and the within-chunk radius drift e^{-TL/n} is not small at n = 291 unless T <~ 10.
- **EM continuum around the peak** extends to |y| < 2pi/5 cleanly; beyond that the ray passes the log-singularities of Q
  (individual near-zero factors). A one-sided (smoothed, log(|.|^2+delta^2)/2) EM costs O(1) per singularity, i.e. O(y),
  which is fine for y <= eps n with small eps — this plausibly bridges |y| in [1, eps n] to the chunk bound at |x| >= eps/5,
  but it was **not** written out or certified.
The missing lemma is therefore precisely: an explicit, n-uniform upper bound for log|S_n(r e^{i theta})| - n p(L)
on the intermediate scale 1 <~ |y| <~ eps n around each 5th root, glued to a chunk bound for |x| >= eps/5.

## Handover
Not changed: no mu_1 is certified. The Laplace certificates (classes 0,1,2 V >= 2; classes 3,4 v >= 2.5) remain the only sound coverage;
the centre (mu >= 0.105) is uncovered for every class.

## Audit items C.2.1-7 status
1, 2, 3 (tail): open (need the lemma above). 4: constant corrected (item 2 above), no rerun. 5, 6: fix is the finite-sum EM of item 1
(derivatives enclosed over the whole window); not implemented. 7: resolved in principle by item 1 (y-parametrisation), not certified.
