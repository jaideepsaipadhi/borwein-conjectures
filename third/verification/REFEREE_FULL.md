# Referee report: THIRD_BORWEIN_PROOF.md (full read)

Scope: I read the whole document first, then spot-checked logs and two scripts. The certificates are treated as black boxes. After writing the findings I looked at the other .md files only to see whether issue M1 had already been caught. It had not.

## Verdict

**Correct with gaps to fill.** The logical architecture is sound: the reduction, the three-region partition, the flow of hypotheses and the use of §3 all check out. There is **one genuine coverage hole**, M1 below. It is tiny and trivially repairable, but it is real, and no audit caught it. There are also several stale or inconsistent passages a journal would require fixed.

## Major issues

**M1. Laplace classes 0–2 leave V in (11.999999999999831, 12) uncovered. This is the same float-drift bug fixed for lap34, but not fixed here.**
- `cls012_lap.py` (lines 135–142) steps with `v2=min(v+w,Vb)`. It has no snap to Vb.
- Replaying the loop gives a last tile of `(11.979999999999832, 11.999999999999831)`. The tail cell is run from `tail 12`, which `cls012_lap_rerun.log` confirms ("tail cell V >= 12").
- A saddle root V in that interval is therefore certified by nothing. §8 item 2 ("V lies in one of the 500 tiles … the tiles are contiguous") is false for classes 0–2 as the script actually ran. B.2 records the drift fix for lap34 only, and GAP_AUDIT012.md does not mention it.
- **Fix:** add the lap34 snap (`if Vb-v2<1e-9: v2=Vb`) or run `cls012_lap.py tail 11.99`, and cite the new log. The tail margins, about 0.86 or more, leave room for this.

(Nothing else rises to "blocks acceptance". The items below are presentation and consistency.)

## Checks performed (all passed unless noted)

1. **Statement and reduction.**
   - The product, the signs and the non-strict inequalities are correct.
   - D = 10n² is correct, and the factor (−1)^{4n} = 1.
   - 10n² − m ≡ −m (mod 5), which swaps classes 1↔4 and 2↔3. Classes 1–4 all have sign ≤ 0, so the pairing is harmless.
   - m = 5n² is its own image. It lies in class 0 with μ = 1/2 and L = 0. It is covered by the centre, via the near-peak L-box [0, 0.01], which is valid down to L = 0.
   - For (F), the statement says L ∈ (0, 2.1], but the script's L-grid is `i/100, i=0..210` and the μ-drift range is [0, μ_max], so L = 0 is included. Only the wording should change (m3).
   - Zeros: all checks treat only strict wrong signs as violations, which is consistent with the conjecture.
2. **Coverage.** The partition is t<4 / (t≥4, μ<0.105) / μ≥0.105, with the boundaries t = 4 and μ = 0.105 assigned correctly.
   - Band, classes 3–4: m<N gives c = 0. N ≤ m < 2N gives c = G(T) ≤ 0 for all n. For 2N ≤ m < 4N, `bandfix_exact` covers n ≤ 2000 (inclusive) and Ψ covers N ≥ 10001, i.e. n ≥ 2001.
   - Band, classes 0–2: exact for n ≤ 2000; Parts A and B for n > 2000. The N_min device is correct, and c(7) = 0 is handled.
   - Laplace, classes 0–2: X ≥ 4V/α ≥ 4V/a₅, using α < a₅ (γ, Q > 0) and m − 1/6 ≥ 4VX + 23/6.
   - Laplace, classes 3–4: X > 3V/a₅. The tail uses X* from α ≤ a₅, and V ≤ a₅X/3.
   - The §6.5 root existence and V > 2 are correct: the right-hand side is below 0.042 and min A = 0.04424.
   - V = 12: a hole for classes 0–2 (M1). For classes 3–4 the tiles snap to 12 and the tail runs from 11.99, with ε_g evaluated at V_T = 11.99 in the script. The doc text "ε_g(12)" in §6.4 is stale (m4).
   - Centre: μ(2.1) < 0.105 and μ(0) = 1/2, so L ∈ [0, 2.1] strictly.
   - θ-coverage of (F) is complete. The pieces are [0.0499, 0.1], [0.1, 1.1566], [1.1566, 1.2067], [1.3065, 1.3566], [1.3566, 2.4133], [2.4133, 2.4634], [2.5632, 2.6133] and [2.6133, π⁺]. That is 6·5 + 53 + 53 + 27 = 163 rectangles. The excluded windows are exactly 2πh/5 ± 0.0499, and (N) reaches 0.05.
   - The (N-a)/(N-b) handover needs n ≥ 800. That is satisfied.
   - Thresholds 291, 821 and 2000/2001 are all consistent with the exact range n ≤ 979.
3. **Hypothesis flow.**
   - Band analytic: needs n > 2000. It is used only there.
   - Laplace 0–2: needs n ≥ 291 and t ≥ 4. Supplied.
   - Laplace 3–4: needs n ≥ 821 and t ≥ 4. Supplied.
   - Centre: needs n ≥ 821 and L ∈ [0, 2.1]. Supplied.
   - Every region uses its representation exactly: §6.1 holds for any W > 0, and the Cauchy formula holds for any r.
4. **§3.**
   - Proved for every ν ≥ 1. N = K(ν) ≥ 5 is correct.
   - I checked the exponent bookkeeping for (3.2) and (3.3), the chord bounds, and the weight counts 4/25 and 1.
   - Uses: §6.1 ρ needs i = m − j ≥ 1, with s(0) handled separately. §5 uses it for i > 200000, and in L1 for 5t ≥ 1000. Both are within the proved range.
5. **Logs spot-checked (13), all matching the document:**
   - `bandfix_exact.log`
   - `cls012_band.log` (5.058e-10; Q ≤ 5.139e-20)
   - `cls012_cov.log`
   - `cls012_lap_rerun.log` (3.0700/1.6879/0.8339; tail 3.0956/1.7136/0.8595)
   - `lap34_cert.log`, refreshed (0.8363/0.8361 on [11.98, 12.0]; tail from 11.99, 0.81456)
   - `centre3_comb_821.log` (final True; minslack 0.018703, OFF 0.008660 at L = 2.09)
   - `offpeak_em_821`, `em0`, `tb`, `chunk_821`
   - `arb_F_summary_full.log` (all 8 pieces OK, 47.91)
   - `exact_1_150` and `exact_150_291`
   - all five `centre3_exact_*` logs (920–1000 has 60 `wrong_sign=0` lines ending at n = 979)
6. **Theorem 3.1 usage.** No section relies on a withdrawn component. I found none.

## Minor issues and presentation

- **m1. Stale text left after the addendum.** Many sentences contradict the addendum and the status table:
  - §7.2, lines 699–706: "`arb_F_summary.log` does not exist"; "closed when …".
  - C.8: "until its Arb re-run completes".
  - §4 caveat, A.1 and the Status table: "no log" for n ≤ 290.
  - §6.4 log caveat, A.6 caveat and the Status row: "predates the fixes".
  - B.1: band classes 3–4 "not independently re-audited".
  - Summary vs body: "four parts" vs "three regions".

  The whole document should be regenerated with the addendum merged in.
- **m2. The wrong file is named as the (F) certificate.** A.8 names `arb_F_summary.log` as the log that must print CERTIFIED. That file actually ends with `(F) arb: INCOMPLETE`. The CERTIFIED line is in `arb_F_summary_full.log`. The addendum's "all 13 logs" should also list which files those are.
- **m3.** Lemma 7.3 and piece (F) are stated for L ∈ (0, 2.1], but they are needed and certified on [0, 2.1]. Fix the statement.
- **m4.** §6.4 says the tail uses ε_g(12). The rerun from 11.99 uses ε_g(11.99), per the script and the log (epsg 1.787e-5). Update the text.
- **m5. Certificates with no log in the directory:**
  - `centre3_saddle_end.py` (μ(2.1) < 0.105; the whole centre range depends on it);
  - `bandfix_analytic.py` (Lemma L1 and Ψ(10001) ≤ 0.198);
  - `bandfix_coeffs.py` (Lemma 5.2).

  Add their logs and must-print lines to Appendix A.
- **m6.** Lemma 5.1 says "P₅ is a power series in its argument with positive coefficients". It is a power series in x = i − 1/6, not in i. The conclusion, that P₅′ is increasing, still holds.
- **m7.** §5.1: the ratio "1.13·10⁻⁴" is quoted, but the log prints only "0.0001". Either print more digits or quote what is printed.
- **m8. Details a journal would want written out:**
  - the Rankin bound for the E_con part of ρ in §6.3 (it is only listed);
  - the derivation of the bound on ∫|Φ_c″| along the ray in §6.2;
  - Part II of piece (N), which is only sketched;
  - an explicit statement that the (F) chunk lemma handles n not divisible by T. The leftover factors are covered by the 4(T−1)·log term, which should be said.
- **m9.** Theorem 8.1 says "exactly one" certificate, but the regions for n ≤ 979 overlap the exact computation. This is harmless, but "at least one" is the claim actually needed.

## Afterwards: did the other files resolve anything?

I searched all .md files for "float drift" and "11.9999". The only hits are in AUDIT_LAP34.md, which covers classes 3–4 only. GAP_AUDIT012.md does not mention it. M1 stays open until the classes 0–2 tail is rerun from 11.99 or the tiling is snapped.
