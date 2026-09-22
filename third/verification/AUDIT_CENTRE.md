# AUDIT_CENTRE: audit of the centre certificate (GAP_CENTRE2 piece (A) and its combination with (B))

*19 Sept 2026.* **Verdict: sound after one stated fix.** The combination bookkeeping as written is wrong: the
uniform bound on OFF is larger than the minimum slack. The per-L-box combination closes with margin at least 0.0100.
That fix is done and logged (`centre3_comb.py`, `centre3_comb_821.log`).

## Defect found

**D1. The uniform combination fails (GAP_CENTRE2 §3, GAP_OFFPEAK §0).**
- With a = p''(L), |OFF| <= 10 pi n sqrt(n a/2pi) e^{-(1.5 log n + 6.3)} = 10 pi sqrt(a/2pi) e^{-6.3}. The powers of n cancel exactly, so this bound is n-free, and it grows with a.
- The largest a is at L = 0, where a = p''(0) = 2/3 (the box enclosure gives a_up = 0.6843). This gives |OFF| <= **0.018791** with the exact 2/3, and 0.019039 with the box enclosure.
- The certificate's own minimum slack min_s(M_s,lo - err_up) is **0.018703**: class 4, box L in [2.09, 2.10], n = 821.
- 0.0187 is a rounded-down value, and it is not less than 0.018703. The claimed inequality "|OFF| <= 0.0187 <= min(M_s - err)" is false.
- The failure comes from pairing the OFF bound at the L with the largest a (L = 0) against the slack at L = 2.1. Paired at the same L, it passes easily.

**Fix (done).** `centre3_comb.py` is `centre3_peak.py` plus one extra output line per box. For each L-box it computes OFF_up = 10 pi sqrt(a_up/2pi) e^{-6.3} in arb and compares it with min_s(M_s,lo - err_up) on the same box.
- Result at N0 = 821 over all 210 boxes and all 5 classes: the minimum of (slack - OFF) is **+0.0100**, at L = 2.09 (slack 0.01870, OFF 0.00866).
- Using the actual certified (B) slack of at least 0.36 (piece E0; e^{-6.66}) instead, the minimum is +0.0127.
- The run takes about 11 min. Command: `python3 centre3_comb.py 821 210 64 > centre3_comb_821.log`.
- **GAP_CENTRE2 §3 should be restated as a per-L comparison.** The "<= 0.0187" wording should be deleted.

## Checklist

1. **The seven old findings (GAP_CLASS012 C.2).** All are fixed, or moved into (B).
   - C.2.1, strip covered by nothing. **Fixed.**
     - The windows are in n-scaled y: inner |y| <= 0.35, strip 0.35-1.1, (B) |y| > 1.1.
     - The (B) pieces tile with overlap: (E) 1.1-6, (N) 6 to 0.25n, (F) from dist 0.0499, and (E0) at h = 0 up to 20.
     - All pieces scale with n, so nothing is left at a fixed theta-width.
   - C.2.2, tail normalised with a lower curvature. **Fixed** in piece (A).
     - The strip and the Mills tail use sqrt(n a_up/2pi).
     - The Gaussian moments use the pointwise a, or A = a_lo, in the safe direction.
     - D1 is a new instance of the same kind of bookkeeping slip, in the combination step.
   - C.2.3, unproved tail lemma. **Replaced by (B)** (GAP_OFFPEAK). It is not re-audited here, except for its statement and normalisation, which are correct. Caveat: (F) uses float64 intervals.
   - C.2.4, EM constant. **Fixed.**
     - B2(1/5) = B2(4/5) = 1/150 and B2(2/5) = B2(3/5) = -11/150, checked by hand, so the per-c constants |B2|/2 are 1/300 and 11/300.
     - The remainder uses |B~2|/2 <= 1/12, and the code uses exactly these values.
     - The derivative chain is f' = (sigma/n) Phi', and int_0^n |f''| <= (|sigma|^2/n) sup|Phi''|. Both are correct.
     - The main term B1(c/5)(Phi_c(sigma) - Phi_c(0)) gives C_h, which is correct.
   - C.2.5, derivatives taken at the peak. **Fixed.** k1, k2, k3, k4 and sup|Phi''| are enclosed over full (L, y) boxes and over t in [0, 1] (64 t-pieces).
   - C.2.6, arctan law and equal modulus. **Fixed.** e^{C_h(L)} is enclosed exactly, and the finite-n error is inside K/n.
   - C.2.7, L in (0, 0.05). **Fixed.** The L-boxes start at 0, and p'(0) = -1 gives mu = 1/2 exactly (Q'(0) = sum_c zeta^c/(1-zeta^c) = -2).
2. **The EM identity, 11/300 and the per-class values.** These are correct, as detailed under C.2.4 above. The per-c |B2(c/5)|/2 is used in the code, which is sharper than the uniform 11/300 and valid.
3. **Uniformity in n. Correct.** n enters only in these places:
   - K/n, as e^{K/n} - 1.
   - The 1/n moment terms (3k4/(nA^2), k2/(nA), 15k3^2/(n c^3), k1^2/(n c)).
   - The Mills tail, e^{-nAY^2/2}/sqrt(n).
   - The strip, sqrt(n) e^{-nD} for n >= 1/(2D), which is asserted.
   - The number of k-terms, which is exactly n per residue in the finite EM. That identity is exact for every n.

   The saddle L solves p'(L) = -2 mu, which depends on mu only, and the windows are in the scaled y. The OFF bound is n-free, as shown in D1.
4. **The saddle. Correct.**
   - p'(0) = -1 gives mu = 1/2.
   - `centre3_saddle_end.py` encloses p'(2.1) with the real mean-value theorem on 4096 t-pieces, which is valid because Q' is real on reals. This gives mu(2.1) < 0.105.
   - a > 0 is asserted on every box, so mu(L) is strictly decreasing and every mu in [0.105, 1/2] has a unique L in [0, 2.1].
   - Integer m enters only through mu = m/(10n^2) and through the class m mod 5.
5. **Combination.** I re-derived the representation c(m) = (1/2pi) int S_n(re^{i theta}) r^{-m} e^{-im theta} d theta with d theta = dy/(5n):
   - r^{-m} = e^{2 n mu L}.
   - The linear phase cancels at the saddle, because -i n p'(L) y = 2i n mu y.
   - Each of the 4 peaks contributes N = e^{n p + 2 n mu L} sqrt(2pi/(na))/(10 pi n) times zeta^{-hm} e^{C_h}.
   - The off-peak theta-set has measure at most 2pi, which gives exactly |OFF| <= 10 pi n sqrt(n a/2pi) e^{-E(n)}.

   The normalisation is consistent with the slack. The uniform comparison fails (D1). The per-box comparison passes with margin 0.0100.
6. **Exact runs.**
   - Method: exact fmpz shift-subtract from S_0 = 1, truncated at x^{5NB^2+1}. The truncation is exact because each factor only feeds higher degrees.
   - Palindromy: c(10n^2 - m) = c(m), because there are 4n factors (an even number) and the degree 10n^2 is 0 mod 5. So m <= 5n^2 suffices.
   - Every m in [0, 5n^2] is checked.
   - Logs, with per-n line counts matching the ranges:

     | log | lines | result |
     |---|---|---|
     | 291-410 | 120 | DONE, 0 violations |
     | 410-580 | 171 | DONE, 0 violations |
     | 580-820 | 241 | DONE, 0 violations |
     | 820-920 | 101 | DONE, 0 violations |
     | 920-1000 | 60 (n = 920..979), then killed (out of memory) | every line wrong_sign = 0 |

   - The segments abut, and exact coverage is n in [291, 979]. GAP_CENTRE2 §1 is stale: it should say that 820-920 has completed and that 920-1000 was killed at 979.
   - Zeros are allowed, and only strict wrong signs count. Zeros occur only at mu < 2e-3, outside the centre.
7. **Independent validation.** Exact S_300, mpmath saddle and closed-form p. The script compares c(m)/N with sum_h zeta^{-hm} e^{C_h(L)}.
   - Points: mu in {0.105, 0.12, 0.15, 0.2, 0.3}, all 5 classes.
   - |difference| <= 0.0045 in every case, against a certified err of about 0.06-0.10. The signs are correct, for example class 4 at mu = 0.105: exact -0.12555, main -0.12609.
   - This confirms the main term and the normalisation N, which is the factor that matters for D1.
   - Not done: n = 821 exact (the build is about 30+ min), and mu >= 0.4, where the mpmath root-finder went complex. This was a script issue, not a certificate issue.

## Minor remarks (no effect on soundness)
- t_zt uses no sqrt(a_up/a_lo) factor. That is correct: the normalised moments use the same pointwise a.
- The certificate at N0 = 1001 cannot be used together with the exact range, which ends at 979. It is not needed.
- GAP_OFFPEAK §7 says it takes "(B) with 6.3" from GAP_CENTRE2. That constant is adequate only with the per-box comparison (D1).
