#!/bin/bash
# Regenerates every certificate and log from scratch. Single core: ~3 hours. Set JOBS for the uniform grid.
set -u
cd "$(dirname "$0")/certificates"
L=../logs; mkdir -p $L
step(){ echo "=== $1  $(date +%T)"; }
step "1 finite range n<=243";            python3 exact_small_n.py 1 200 > $L/01_exact_small_n_a.log 2>&1; python3 exact_small_n.py 201 225 > $L/01_exact_small_n_b.log 2>&1; python3 exact_small_n.py 226 243 > $L/01_exact_small_n_c.log 2>&1
step "2 band lemma";                     python3 band_lemma.py 36000         > $L/02_band_lemma.log 2>&1
                                         python3 band_tail_rigorous.py       > $L/02b_band_tail.log 2>&1
step "3 expansion (sz_joint) audit";     python3 sz_expansion_audit.py       > $L/03_sz_audit.log 2>&1
                                         python3 expansion_check.py 36000    > $L/03b_expansion_check.log 2>&1
step "4 hand lemmas";                    python3 poisson_check.py            > $L/04a_poisson.log 2>&1
                                         python3 hand_checks.py              > $L/04b_hand_checks.log 2>&1
                                         python3 koksma_cauchy_check.py      > $L/04c_koksma_cauchy.log 2>&1
                                         python3 identity_check.py           > $L/04d_identity.log 2>&1
step "5 centre tables";                  rm -f midu_table.txt g3u_table.txt clowu_table.txt
  python3 centre_mid_uniform.py  -0.3 1.75 82 >> midu_table.txt; python3 centre_mid_uniform2.py 1.75 3.05 26 >> midu_table.txt; python3 centre_mid_uniform3.py 3.05 3.4 7 >> midu_table.txt
  python3 centre_g3_uniform.py   -0.3 1.75 82 >> g3u_table.txt;  python3 centre_g3_uniform2.py  1.75 3.05 26 >> g3u_table.txt;  python3 centre_g3_uniform3.py  3.05 3.4 7 >> g3u_table.txt
  python3 centre_clow_uniform.py -0.3 3.4 74 >> clowu_table.txt
step "6 uniform Laplace grid";           rm -f $L/uni_*.log
python3 - > /tmp/chunks_release.txt <<'PY'
def add(a,b,n):
    for i in range(n): print('%.9f %.9f'%(a*(b/a)**(i/n), a*(b/a)**((i+1)/n)))
add(2.95,4.0,10); add(4.0,20.0,16); add(20.0,3000.0,10)
PY
cat /tmp/chunks_release.txt | xargs -P ${JOBS:-1} -L1 bash -c 'python3 run_uniform.py $0 $1 14 > ../logs/uni_$0_$1.log 2>&1'
step "7 closing lemma";                  python3 ray_asymptotic2.py          > $L/07_closing_lemma.log 2>&1
step "8 coverage";                       python3 coverage.py                 > $L/08_coverage.log 2>&1
step "done"; tail -1 $L/08_coverage.log
