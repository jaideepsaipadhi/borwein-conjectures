# Acceptance test for the centre regime: exact c(m) vs the centre main term 2 Re[ e^{L(u*)-m u*} / sqrt(2 pi L''(u*)) ] (conjugate pair).
import numpy as np, math, sys
exec(open('exact_crt.py').read().split('if __name__')[0])
n=244; N=3*n+1; k=np.array([j for j in range(1,3*n+1) if j%3],dtype=float)
mus=[0.2,0.3,0.4,0.45,0.49]; ms=[]
for mu in mus:
    m=round(mu*N*N); m-=(m-2)%3; ms.append(m)
ps,res=coeffs(n,max(ms))
for m in ms:
    u=2j*np.pi/3-0.5/N
    for _ in range(200):
        z=np.exp(k*u); L1=3*np.sum(-k*z/(1-z)); L2=3*np.sum(-k*k*z/(1-z)**2); u=u-(L1-m)/L2
    z=np.exp(k*u); L=3*np.sum(np.log(1-z)); L2=3*np.sum(-k*k*z/(1-z)**2)
    logM=L-m*u-0.5*np.log(2*np.pi*L2)      # complex log of single-peak contribution
    c=crt_at(ps,res,m)
    lc=(abs(c).bit_length()-60)*math.log(2)+math.log(abs(c)>>(abs(c).bit_length()-60))
    pred_log=logM.real+math.log(2*abs(math.cos(logM.imag)))
    sgn=np.sign(math.cos(logM.imag))
    print('m=%6d mu=%.3f exact sign %s, predicted sign %s, |c|/|main| = %.5f'%(m,m/N**2,'-' if c<0 else '+','-' if sgn<0 else '+',math.exp(lc-pred_log)),flush=True)
