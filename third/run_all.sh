#!/bin/bash
# Regenerates every certificate log for the Third Borwein conjecture (Theorem B).
# Logs are written to ../logs relative to certificates/.  Steps can be skipped:
#   SKIP_EXACT=1   skip the exact range n <= 979 (about 2.4 h CPU; the last segment needs ~6 GB RAM)
#   SKIP_F=1       skip the Arb run of off-peak piece (F) (about 1.5 h on 2 cores)
#   JOBS=k         worker processes for piece (F) (default 2)
set -u
cd "$(dirname "$0")/certificates"
L=../logs; mkdir -p $L
step(){ echo "=== $1  $(date +%T)"; }
[ -f g5_coeffs.pkl ] || { step "0 coefficient table s(i), i<=200000"; python3 g5_coef.py; }
if [ -z "${SKIP_EXACT:-}" ]; then
step "A.1-A.7 exact range n<=979"
  python3 centre3_exact.py 1 150    > $L/exact_1_150.log 2>&1
  python3 centre3_exact.py 150 291  > $L/exact_150_291.log 2>&1
  python3 centre3_exact.py 291 410  > $L/centre3_exact_291_410.log 2>&1
  python3 centre3_exact.py 410 580  > $L/centre3_exact_410_580.log 2>&1
  python3 centre3_exact.py 580 820  > $L/centre3_exact_580_820.log 2>&1
  python3 centre3_exact.py 820 920  > $L/centre3_exact_820_920.log 2>&1
  python3 centre3_exact.py 920 1000 > $L/centre3_exact_920_1000.log 2>&1   # expected to be killed after n=979 unless >6 GB RAM
fi
step "A.8-A.11 band";         python3 bandfix_exact.py 291 2000 > $L/bandfix_exact.log 2>&1
                              python3 bandfix_coeffs.py         > $L/bandfix_coeffs.log 2>&1
                              python3 bandfix_analytic.py 10001 > $L/bandfix_analytic.log 2>&1
                              python3 cls012_band.py            > $L/cls012_band.log 2>&1
step "A.12-A.14 Laplace 0-2"; { python3 cls012_lap.py 2.0 12.0 0.02; python3 cls012_lap.py tail 11.99; } > $L/cls012_lap_rerun.log 2>&1
                              python3 cls012_cov.py             > $L/cls012_cov.log 2>&1
step "A.15 Laplace 3-4";      { python3 lap34_cert.py 2.0 12.0 0.02; python3 lap34_cert.py tail 11.99; } > $L/lap34_cert.log 2>&1
step "A.16 saddle endpoint";  python3 centre3_saddle_end.py     > $L/centre3_saddle_end.log 2>&1
step "A.17 near peak + combination"
  { python3 centre3_comb.py 821 210 64 0 105; python3 centre3_comb.py 821 210 64 105 210; } > $L/centre3_comb_821.log 2>&1
step "A.18-A.20 off-peak E, E0, N"
  python3 offpeak_em.py 821 6.0 32  > $L/offpeak_em_821.log 2>&1
  python3 offpeak_em0.py 821 20 32  > $L/offpeak_em0_821.log 2>&1
  python3 offpeak_tb.py 821 0.25    > $L/offpeak_tb_821.log 2>&1
if [ -z "${SKIP_F:-}" ]; then
step "A.21 off-peak F (Arb)"
  for p in 0 1 2 3 4 5 6 7; do python3 offpeak_chunk_arb.py 821 $p -j ${JOBS:-2} > $L/rerun_arb_F_p$p.log 2>&1; done
  python3 offpeak_chunk_arb.py summarize $L/rerun_arb_F_p?.log > $L/rerun_arb_F_summary.log 2>&1
  tail -1 $L/rerun_arb_F_summary.log
fi
step "done"
grep -h "CERTIFIED\|violations 0\|0 failures\|COMBINATION ASSERTED\|True" $L/*.log | head -40
