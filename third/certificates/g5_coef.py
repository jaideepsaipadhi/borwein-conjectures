import time, pickle
from flint import fmpz_poly
M=200000
t0=time.time()
pent=[0]*(M+1); pent[0]=1
j=1
while True:
    a=j*(3*j-1)//2; b=j*(3*j+1)//2
    if a>M and b>M: break
    sg=-1 if j%2 else 1
    if a<=M: pent[a]+=sg
    if b<=M: pent[b]+=sg
    j+=1
L=M//5
p=[0]*(L+1); p[0]=1
for n in range(1,L+1):
    s=0; j=1
    while True:
        a=j*(3*j-1)//2; b=j*(3*j+1)//2
        if a>n and b>n: break
        sg=-1 if j%2==0 else 1
        if a<=n: s+=sg*p[n-a]
        if b<=n: s+=sg*p[n-b]
        j+=1
    p[n]=s
print("p(n) to %d in %.0fs"%(L,time.time()-t0),flush=True)
part5=[0]*(M+1)
for n in range(L+1): part5[5*n]=p[n]
G=(fmpz_poly(pent)*fmpz_poly(part5)).truncate(M+1)
co=[int(G[i]) for i in range(M+1)]
print("G_5 coeffs in %.0fs; s(0..8)=%s"%(time.time()-t0,co[:9]),flush=True)
pickle.dump(co,open('g5_coeffs.pkl','wb')); print("saved")
