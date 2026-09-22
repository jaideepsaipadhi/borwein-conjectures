# Certificate 2 (band 3n+2 <= m < 6n+2 = 2N): exactly, E_n = 1 + 3 sum_{N<=j<2N, 3 not| j} q^j + O(q^{2N}), hence for m = 3M+2 < 2N
#   c_n(m) = 3 sum_{i=0}^{M-n} g(i),   g(i) = p3(3i) + p3(3i+1),   p3 = [q^i] P_infty^3 .
# (a) identity checked against exact c_n(m) for n <= 60;  (b) g(i) < 0 exactly for i <= IMAX (Jacobi triple product + 3-colour partitions);
# (c) g(i) < 0 for i >= 100 from the sz_joint expansion (band_tail.py).  Since M - n < n, (b)+(c) give c_n(m) < 0 for all n.
import sys, time
from flint import fmpz_poly
IMAX=int(sys.argv[1]) if len(sys.argv)>1 else 36000
J=IMAX+2; t0=time.time()
sig=[0]*(J+1)
for d in range(1,J+1):
    for mm in range(d,J+1,d): sig[mm]+=d
b=[0]*(J+1); b[0]=1
for j in range(1,J+1):
    s=0
    for t in range(1,j+1): s+=sig[t]*b[j-t]
    b[j]=3*s//j
def p3(i):
    s=0; k=0
    while k*(k+1)//2<=i:
        r=i-k*(k+1)//2
        if r%3==0: s+=(-1)**k*(2*k+1)*b[r//3]
        k+=1
    return s
ok=True
P=fmpz_poly([1])
for n in range(1,61):
    for k in (3*n-2,3*n-1): P=P*fmpz_poly([1]+[0]*(k-1)+[-1])
    c=(P**3).coeffs(); N=3*n+1
    for m in range(3*n+2,2*N,3):
        M=(m-2)//3
        if int(c[m])!=3*sum(p3(3*i)+p3(3*i+1) for i in range(M-n+1)): ok=False
print('(a) band identity c_n(m) = 3 sum g(i) for all n<=60, 3n+2<=m<2N:', ok)
negall=all(p3(3*i)+p3(3*i+1)<0 for i in range(IMAX+1))
print('(b) g(i) < 0 exactly for all i <= %d: %s  (%.0fs)'%(IMAX,negall,time.time()-t0))
