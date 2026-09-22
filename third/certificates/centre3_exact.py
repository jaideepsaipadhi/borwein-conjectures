# Exact sign check of ALL coefficients of S_n(q)=prod_{k<=5n,5nmid k}(1-q^k), n in [NA,NB].
# S_n is palindromic of degree 10n^2 (c(D-m)=c(m), D=0 mod 5, class s <-> -s), so m<=5n^2 suffices.
# Works mod x^{5NB^2+1}: S_{n+1}=S_n*prod_{k=5n+1..5n+4}(1-x^k), each factor by shift-subtract (exact fmpz).
# Usage: python3 centre3_exact.py NA NB   (keep NB <= ~1.4*NA for memory efficiency)
# Output: one line per n: wrong-sign count (strict), zero count and largest zero index (zeros occur only at tiny m).
import sys, time
from flint import fmpz_poly
NA,NB=int(sys.argv[1]),int(sys.argv[2]); T=5*NB*NB+1
def build(n):   # incremental from S_0 = 1 (memory-light; a mul_low product tree OOMs at n ~ 820 on 8 GB)
    P=fmpz_poly([1])
    for j in range(n):
        for k in range(5*j+1,5*j+5):
            P=P-P.left_shift(k)
        P=P.truncate(T)
    return P
def check(P,n):
    H=5*n*n; bad=[]; zeros.clear()   # index P[m] directly: no full coefficient list (memory)
    for m in range(H+1):
        v=P[m]
        if m%5==0:
            if v<0: bad.append(m)
        elif v>0: bad.append(m)
        if v==0: zeros.append(m)
    return bad
t0=time.time(); P=build(NA); print('built S_%d mod x^%d in %.0fs'%(NA,T,time.time()-t0),flush=True)
n=NA; tot=0; zeros=[]
while True:
    t=time.time(); bad=check(P,n); tot+=len(bad)
    print('n=%d wrong_sign=%d %s zeros=%d max_zero_m=%s mu_of_it=%s (%.1fs, total %.0fs)'%(n,len(bad),bad[:5],len(zeros),max(zeros) if zeros else None,'%.2e'%(max(zeros)/(10*n*n)) if zeros else None,time.time()-t,time.time()-t0),flush=True)
    if n==NB: break
    for k in range(5*n+1,5*n+5):
        P=P-P.left_shift(k)
    P=P.truncate(T)
    n+=1
print('DONE n=%d..%d total violations %d, %.0fs'%(NA,NB,tot,time.time()-t0))
