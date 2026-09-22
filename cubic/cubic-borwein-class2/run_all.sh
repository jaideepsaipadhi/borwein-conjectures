#!/usr/bin/env bash
# Reproduce every certificate of the proof from scratch.  Single core total: ~2 h (the Laplace grid dominates).
# Usage:  ./run_all.sh            (logs in logs/, summary in logs/SUMMARY.txt)
set -e
cd "$(dirname "$0")"; ROOT=$(pwd); mkdir -p logs
echo "== environment" | tee logs/SUMMARY.txt
python3 --version | tee -a logs/SUMMARY.txt; python3 -c "import flint,mpmath,numpy,scipy,sympy;print('python-flint',flint.__version__,'mpmath',mpmath.__version__,'numpy',numpy.__version__,'scipy',scipy.__version__,'sympy',sympy.__version__)" | tee -a logs/SUMMARY.txt
echo "== 1. finite range n <= 243 (GMP)" | tee -a logs/SUMMARY.txt
( cd exact && gcc -O2 -o cubic_check cubic_check.c -lgmp && ./cubic_check 1 243 > ../logs/finite.log )
python3 - <<'PY' | tee -a logs/SUMMARY.txt
import re
bad=[l for l in open('logs/finite.log') if 'violations' in l and not (' violations=0 ' in l and 'zeros_in_3n_to_D/2=0 ' in l and 'nonzero_below_3n=0' in l)]
n=sum(1 for l in open('logs/finite.log') if 'violations' in l)
print('finite: %d values of n checked, lines with violations: %d'%(n,len(bad)))
PY
cd certificates
echo "== 2. band lemma" | tee -a ../logs/SUMMARY.txt
python3 band_certificate.py 2000 > ../logs/band.log; tail -1 ../logs/band.log | tee -a ../logs/SUMMARY.txt
echo "== 3. centre tables (uniform in N)" | tee -a ../logs/SUMMARY.txt
rm -f midu_table.txt g3u_table.txt clowu_table.txt
python3 centre_mid_uniform.py  -0.3 -0.05 5  >> midu_table.txt
python3 centre_mid_uniform.py  -0.05 1.75 72 >> midu_table.txt
python3 centre_mid_uniform2.py  1.75 3.05 26 >> midu_table.txt
python3 centre_mid_uniform3.py  3.05 3.4 7   >> midu_table.txt
python3 centre_g3_uniform.py   -0.3 1.75 82  >> g3u_table.txt
python3 centre_g3_uniform2.py   1.75 3.05 26 >> g3u_table.txt
python3 centre_g3_uniform3.py   3.05 3.4 7   >> g3u_table.txt
python3 centre_clow_uniform.py -0.3 3.4 148  >> clowu_table.txt
wc -l midu_table.txt g3u_table.txt clowu_table.txt | tee -a ../logs/SUMMARY.txt
echo "== 4. Laplace grid v in [2.95, 3000] (uniform in N)" | tee -a ../logs/SUMMARY.txt
rm -f ../logs/uni_*.log
python3 - > /tmp/chunks.txt <<'PY'
def add(a,b,n):
    for i in range(n): print('%.9f %.9f'%(a*(b/a)**(i/n), a*(b/a)**((i+1)/n)))
add(2.95,4.0,10); add(4.0,20.0,16); add(20.0,3000.0,10)
PY
cat /tmp/chunks.txt | xargs -P ${JOBS:-1} -L1 bash -c 'python3 run_uniform.py $0 $1 14 > ../logs/uni_$0_$1.log 2>/dev/null'
grep -h UNIFORM ../logs/uni_*.log | tee -a ../logs/SUMMARY.txt
echo "== 5. closing lemma v >= 3000" | tee -a ../logs/SUMMARY.txt
python3 ray_asymptotic2.py > ../logs/closing_lemma.log; tail -1 ../logs/closing_lemma.log | tee -a ../logs/SUMMARY.txt
echo "== 6. coverage of all (n,m)" | tee -a ../logs/SUMMARY.txt
python3 coverage.py > ../logs/coverage.log; cat ../logs/coverage.log | tee -a ../logs/SUMMARY.txt
echo "== 7. hand-lemma checks" | tee -a ../logs/SUMMARY.txt
python3 poisson_check.py > ../logs/poisson.log; tail -1 ../logs/poisson.log | tee -a ../logs/SUMMARY.txt
python3 hand_checks.py > ../logs/hand_checks.log; cat ../logs/hand_checks.log | tee -a ../logs/SUMMARY.txt
python3 koksma_cauchy_check.py > ../logs/koksma_cauchy.log; cat ../logs/koksma_cauchy.log | tee -a ../logs/SUMMARY.txt
python3 identity_check.py > ../logs/identity.log 2>&1 || true; grep "rel.err" ../logs/identity.log | tee -a ../logs/SUMMARY.txt
python3 sz_expansion_audit.py > ../logs/sz_audit.log; cat ../logs/sz_audit.log | tee -a ../logs/SUMMARY.txt
python3 expansion_check.py 5000 > ../logs/expansion_check.log; tail -1 ../logs/expansion_check.log | tee -a ../logs/SUMMARY.txt
echo "== done" | tee -a ../logs/SUMMARY.txt
