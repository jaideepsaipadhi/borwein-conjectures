cd "$(dirname "$0")/certificates"
while read lo hi; do f=../logs/uni_${lo}_${hi}.log; if [ -f $f ] && grep -q UNIFORM $f; then continue; fi; python3 run_uniform.py $lo $hi 14 > $f.tmp 2>&1 && mv $f.tmp $f; done < ../chunks.txt
