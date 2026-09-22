# End-to-end machine check of the H5 chain:  c(m) = -2 (Im J0 + Im J_eps) + R(m)
#   c(m) exact from P_n^3;  R(m) = sum_j e(j) (p3 - sqrt3 chi P3)(m-j) exact/mpmath;
#   Im J0, Im J_eps: periodic trapezoid on Re w = W, |y| <= pi, integrand from the closed forms S, Delta and the e^{-w/4} Poisson form.
import numpy as np, mpmath as mp, math, sys, time
a=math.pi**2/9
def run(n, ms, K=1<<14, KE=2000):
    t0=time.time(); N=3*n+1; M=max(ms)
    P=[0]*(3*n*n+1); P[0]=1
    for k in range(1,3*n+1):
        if k%3==0: continue
        for _ in range(3):
            for i in range(len(P)-1,k-1,-1): P[i]-=P[i-k]
    p3=[0]*(M+1); p3[0]=1
    for k in range(1,M+1):
        if k%3==0: continue
        for _ in range(3):
            for i in range(M,k-1,-1): p3[i]-=p3[i-k]
    e=[0]*(M+1); e[0]=1
    for Kk in range(3*n+1,M+1):
        if Kk%3==0: continue
        for _ in range(3):
            for i in range(Kk,M+1): e[i]+=e[i-Kk]
    chi=lambda i:(1,-1,0)[i%3]
    def P3(i):
        if i<1: return mp.mpf(0)
        return mp.sqrt(a/(i-mp.mpf(1)/4))*mp.besseli(1,2*mp.sqrt(a*(i-mp.mpf(1)/4)))
    y=-math.pi+2*math.pi*np.arange(K)/K
    for m in ms:
        R=float(sum(e[j]*(p3[m-j]-mp.sqrt(3)*chi(m-j)*P3(m-j)) for j in range(m+1)))
        W=math.sqrt(a/(m-N)); w=W+1j*y; q=np.exp(-w)
        S=np.zeros(K,complex); D=np.zeros(K,complex)
        for r in range(1,80):
            qNr=np.exp(-N*r*w); den=1-np.exp(-3*r*w)
            if r%3: D+=3*math.sqrt(3)/(2*r)*(1 if r%3==1 else -1)*qNr*(1-np.exp(-r*w))/den; S+=-1.5/r*(qNr+qNr*np.exp(-r*w))/den
            else: S+=3.0/r*(qNr+qNr*np.exp(-r*w))/den
        base=np.exp((m-0.25)*w+S+1j*D)
        eps=-np.ones(K,complex)
        for k in range(1,KE):
            eps+=(-1j)**k*(np.exp(a/(w+2j*math.pi*k))-1)+(1j)**k*(np.exp(a/(w-2j*math.pi*k))-1)
        I0=np.mean(np.imag(base*np.exp(a/w))); Ie=np.mean(np.imag(base*eps))
        rec=-2*(I0+Ie)+R
        print('n=%d m=%d c=%d reconstructed=%.6f rel.err=%.2e | -2ImJ0=%.5e -2ImJeps=%.3e R=%.3e  [%.0fs]'%(n,m,P[m],rec,abs(rec-P[m])/max(1,abs(P[m])),-2*I0,-2*Ie,R,time.time()-t0),flush=True)
print('start',flush=True)
run(12,[200,302,404],K=1<<13,KE=2000)
run(20,[500,800,1100],K=1<<14,KE=2000)
