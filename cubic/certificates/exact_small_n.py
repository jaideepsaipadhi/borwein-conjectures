# Certificate (finite range): for every n in [NA, NB] and every m = 2 mod 3 with m <= 9n^2/2:
#   c_n(m) = 0 if m < 3n  and  c_n(m) < 0 if 3n <= m <= 9n^2/2.     Exact integer polynomial arithmetic (FLINT fmpz_poly).
# Usage: python3 exact_small_n.py NB            (all n <= NB)
#        python3 exact_small_n.py NA NB         (a chunk; chunks are independent)
import sys, time
from flint import fmpz_poly
if len(sys.argv)>2: NA,NB=int(sys.argv[1]),int(sys.argv[2])
else: NA,NB=1,int(sys.argv[1]) if len(sys.argv)>1 else 243
t0=time.time(); bad=[]
P=fmpz_poly([1])
for k in range(1,3*(NA-1)+1):
    if k%3: P=P*fmpz_poly([1]+[0]*(k-1)+[-1])
for n in range(NA,NB+1):
    for k in (3*n-2,3*n-1):
        P=P*fmpz_poly([1]+[0]*(k-1)+[-1])
    C=P**3
    for m in range(2,9*n*n//2+1,3):
        v=int(C[m])
        if (m<3*n and v!=0) or (m>=3*n and v>=0): bad.append((n,m,v)); break
    print('  n=%d ok=%s (%.0fs)'%(n,not bad,time.time()-t0),flush=True)
print('EXACT n in [%d,%d]: violations %d %s (%.0fs)'%(NA,NB,len(bad),bad[:3],time.time()-t0))
