# Exact c(m) = [q^m] P_n(q)^3 for m <= M via multi-modular numpy (each factor (1-q^k)^3 is a vectorised sparse convolution) + CRT.
import numpy as np, sys
from sympy import prevprime
def coeffs(n,M,nprimes=None):
    bits=6*n+64; nprimes=nprimes or bits//59+2
    ps=[]; p=2**60
    while len(ps)<nprimes: p=prevprime(p); ps.append(p)
    res=[]
    for p in ps:
        c=np.zeros(M+1,dtype=np.int64); c[0]=1
        for k in range(1,3*n+1):
            if k%3==0: continue
            old=c.copy()
            c[k:]=(c[k:]-3*old[:-k])%p
            if 2*k<=M: c[2*k:]=(c[2*k:]+3*old[:-2*k])%p
            if 3*k<=M: c[3*k:]=(c[3*k:]-old[:-3*k])%p
        res.append(c)
    return ps,res
def crt_at(ps,res,m):
    x=0; P=1
    for p,c in zip(ps,res):
        r=int(c[m]); t=((r-x)*pow(P,-1,p))%p; x+=P*t; P*=p
    return x-P if x>P//2 else x
if __name__=='__main__':
    n=int(sys.argv[1]); ms=[int(t) for t in sys.argv[2:]]
    ps,res=coeffs(n,max(ms))
    for m in ms: print(m,crt_at(ps,res,m))
