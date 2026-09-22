# Arb run of piece (F), N0 = 821: debugging and completion

## Cause
`ArbChunk.arb_ok` stopped early with `if p > thr: return False` inside the product loop. The factors
|1 - rho e^{i phi}|^2 can be < 1, so a partial product above the threshold does not mean the full product is
above it. For piece 5 rect 1 at s >= 0.3, the partial product rises above e^{2t} around k' ~ 90 and then drops
below it. Example cell, s = 0.4: the full product is about e^{13.4} and the threshold is e^{17.2}. Every such cell
was rejected, so bb kept subdividing cells that could never pass until it hit maxcells. Then tol doubled up to 8
and the rect was pre-split, which gave the ~1e6-cell / OOM / >580 s behaviour. The float checker has no early
exit, so it accepted all 1344 initial cells for each (s, edge).
Earlier OK results are still valid: the early exit could only reject cells, never accept one.

## Fix (it only makes the check less conservative; no bound is weakened)
1. The early exit is removed. `arb_ok` returns `p < thr` for the full arb product.
2. (Extra tightening, also sound) `cos_lo(ph)`: for the ball m +- r, cos(m+h) >= min(cos m, cos m (1 - r^2/2)) - |sin m| r.
   This is valid because cos h lies in [1 - r^2/2, 1] and |sin h| <= r. The code takes the max of this bound and arb's own cos lower bound.
   On its own it did not fix the blow-up.
Both changes are applied in offpeak_chunk_arb.py and offpeak_chunk_arb_lowmem.py.
After the fix, piece 5 rect 1 needs 61,952 cells (arbfail = 0), exactly as in the float run, with margin 185.04 in 30 s.

## Runs (foreground, 2 parallel)
- arb_F_p0_rest.log: piece 0 rects 1..5, all OK
- arb_F_p5_1-13.log, arb_F_p5_14-19.log, arb_F_p5_20-26.log: piece 5 rects 1..26, all OK, arbfail = 0

## Summary (arb_F_summary_full.log)
`summarize` already merges several logs. It was run on p0, p0_rest, p1..p4, p5a, p5_1-13, p5_14-19, p5_20-26,
p5b, p6 and p7. Every rect of pieces 0..7 is OK and no rect FAILED:
**"(F) arb: CERTIFIED for all n >= N0 (all rects OK)"**.
Minimum margin = **47.91** (piece 7 rect 26, [3.13327, pi+]), the same as the float run.
