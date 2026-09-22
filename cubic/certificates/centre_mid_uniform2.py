# Uniform-in-N bound for the centre-regime middle stretch  MID/main, valid for every N >= N0 = 733, on xi-balls.
#   Contour Re u = log r, r = e^{-xi/N}; theta = theta* + phi, t = N(theta - 2pi/3) = N(phi + delta), |N delta| <= 0.01.
#   Zone 1 (|t| <= t1):  g(phi) = -int_0^phi (phi-s) Re L''(u*+is) ds,  Re L'' >= N^3 rho_lo(t) - 3N^2 V1   (Koksma per class)
#   Zone 2 (t1..10):     G <= -(3/2) sum g(x_k)(2cos k theta + 1) <= -N gam(t) + C2(t)        (log(1-u) <= -u, Koksma)
#   Zone 3 (10..N*0.051): Dirichlet kernel  G <= -(3/2)(T_lo - (2/3)(1.25 + 0.6673 N/t))
#   G(theta*) >= 0 (theta* maximises G on zone 1, where Re L'' > 0), so g <= G.
# Each zone bound is non-increasing in N (asserted below), so the value at N0 bounds every N >= N0.
from flint import arb, acb, ctx
import sys, math
ctx.prec=80
PI=arb.pi(); I=acb(0,1); om=[(2*PI*I*j/3).exp() for j in (1,2)]
N0=733; t1=0.8; CELLS=120; SHIFT=0.015
def wfun(z): return -z/(1-z)**2
def dwdz(z): return -(1+z)/(1-z)**3
def gfun(x): return x/(1+x+x*x)
def cellballs(a,b,K):
    return [arb((a+(b-a)*(i+0.5)/K),(b-a)/(2*K)) for i in range(K)], (b-a)/K
def rho_data(xi,tb):
    # returns (lower bound of sum_j int F_j, V+S sum over j, upper bound of |sum_j int tau^2 w|) for F_j = tau^2 Re w(z_j)
    taus,h=cellballs(0,1,CELLS)
    lo=arb(0); VS=arb(0); cabs=acb(0)
    for j in range(2):
        tv=arb(0); sup=arb(0)
        for tau in taus:
            z=om[j]*(-xi*tau+I*tau*tb).exp()
            w=wfun(z); F=tau*tau*w
            lo+=arb(F.real.lower())*h
            cabs+=F*h
            dF=2*tau*w+tau*tau*dwdz(z)*z*(-xi+I*tb)
            tv+=arb(abs(dF).upper())*h
            sup=max(sup,arb(abs(F).upper()),key=lambda q:q.upper())
        VS+=tv+sup
    return lo, VS, arb(abs(cabs).upper())
def gam_data(xi,tb):
    taus,h=cellballs(0,1,CELLS)
    lo=arb(0); C=arb(0)
    for j in range(2):
        tv=arb(0); sup=arb(0)
        for tau in taus:
            x=(-xi*tau).exp(); ph=2*PI*(j+1)/3+tau*tb
            F=gfun(x)*(2*ph.cos()+1)
            lo+=arb(F.lower())*h
            dg=-xi*x*(1-x*x)/(1+x+x*x)**2
            dF=dg*(2*ph.cos()+1)-gfun(x)*2*ph.sin()*tb
            tv+=arb(abs(dF).upper())*h
            sup=max(sup,arb(abs(F).upper()),key=lambda q:q.upper())
        C+=tv+sup
    return lo, arb(3)/2*C          # G <= -(3/2)(N/3) lo + (3/2) C  = -(N/2) lo + C'
def mid_uniform(xi, cs=(2.5,3.0,3.5)):
    N=arb(N0)
    # zone-1 data on t-balls covering [0, t1+SHIFT] (symmetric in t: use |t| via both signs)
    K1=16; edges=[(t1+SHIFT)*i/K1 for i in range(K1+1)]
    rho_lo=[]; V1=arb(0); a2abs=None
    for i in range(K1):
        best=None
        for sgn in (1,-1):
            tb=sgn*arb((edges[i]+edges[i+1])/2,(edges[i+1]-edges[i])/2+SHIFT)
            lo,VS,cab=rho_data(xi,tb)
            best=lo if best is None or lo.lower()<best.lower() else best
            V1=max(V1,VS,key=lambda q:q.upper())
            if i==0: a2abs=cab if a2abs is None or cab.upper()>a2abs.upper() else a2abs
        rho_lo.append(arb(best.lower()))
    run=[]; m=None
    for r_ in rho_lo: m=r_ if m is None or r_.lower()<m.lower() else m; run.append(m)
    lo0,_,_=rho_data(xi,arb(0,SHIFT))
    a2lo=arb(lo0.lower())                    # |L''| >= Re L'' >= N^3 a2lo - 3N^2 V1
    Ahi=N**3*a2abs+3*N**2*V1; Alo=N**3*a2lo-3*N**2*V1
    assert Alo.lower()>0
    sig_hi=1/Alo.sqrt()
    kap=[(N**3*r_-3*N**2*V1)/Ahi for r_ in run]
    assert all(k_.lower()>0 for k_ in kap), ('zone1 curvature not positive',xi)
    # y-grid: t(y) = N sigma y <= N sig_hi y ; zone 1 ends at Y with N sig_hi Y = t1
    Y=float((t1/(N*sig_hi)).lower()); M=4000; dy=Y/M
    slope=0.0; g=0.0; gs=[0.0]; ys=[0.0]
    for i in range(M):
        y=i*dy; tt=float((N*sig_hi).upper())*(y+dy)
        idx=min(int(tt/((t1+SHIFT)/K1)),K1-1); kk=float(kap[idx].lower())
        g-= slope*dy + kk*dy*dy/2; slope+=kk*dy
        gs.append(g); ys.append(y+dy)
    tail=math.exp(g)/slope
    out={}
    for c in cs:
        s=0.0
        for i in range(M):
            if ys[i+1]<=c: continue
            s+=math.exp(gs[i] if ys[i]>=c else gs[i])*dy     # left endpoint (g non-increasing)
        out[c]=2*(s+tail)/math.sqrt(2*math.pi)
    # zone 2: t in [t1, 10]
    K2=90; e2=[t1*(10/t1)**(i/K2) for i in range(K2+1)]
    I2=arb(0)
    for i in range(K2):
        tb=arb((e2[i]+e2[i+1])/2,(e2[i+1]-e2[i])/2+SHIFT)
        lo,C=gam_data(xi,tb)
        gam=lo/2                                 # G <= -N gam + C
        assert (gam*N).lower()>0.5, ('zone2 not decreasing',xi,e2[i],gam)
        I2+=2*arb(e2[i+1]-e2[i]+0.02)/N*Ahi.sqrt()/(2*PI).sqrt()*(-N*gam+C).exp()
    # zone 3: t in [10, 0.051 N]
    taus,h=cellballs(0,1,CELLS)
    intg=sum((arb(gfun((-xi*tau).exp()).lower())*h for tau in taus),arb(0))
    Tlo=2*N/3*intg-arb(4)/3
    ex3=-arb(3)/2*(Tlo-arb(2)/3*(arb('1.25')+arb('0.6673')*N/(10-SHIFT)))
    coef=intg-arb('0.0668')                   # exponent ~ -N*coef : need N0*coef > 3/2 for decrease with N^{3/2}
    assert (coef*N).lower()>1.5, ('zone3 not decreasing',xi)
    I3=2*arb('0.052')*Ahi.sqrt()/(2*PI).sqrt()*ex3.exp()
    return {c:out[c]+float(I2.upper())+float(I3.upper()) for c in cs}, float(I2.upper()), float(I3.upper()), [float(k.lower()) for k in (kap[0],kap[-1])]
if __name__=='__main__':
    lo,hi,K=float(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3])
    for i in range(K):
        a_=lo+(hi-lo)*i/K; b_=lo+(hi-lo)*(i+1)/K
        xi=arb((a_+b_)/2,(b_-a_)/2)
        res,I2,I3,kp=mid_uniform(xi)
        print('xi in [%.4f,%.4f]: MIDU c=2.5 %.4f  c=3 %.4f  c=3.5 %.4f | zone2 %.2e zone3 %.2e | kappa %.3f..%.3f'%(a_,b_,res[2.5],res[3.0],res[3.5],I2,I3,kp[0],kp[1]),flush=True)
