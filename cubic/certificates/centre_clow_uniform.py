# Uniform-in-N sign margin for the centre main term, on xi-balls (xi = -N log r), valid for all N >= N0 = 733.
#   phase = alpha(r) - 2pi/3 + b + E1   (mod 2pi), where
#   alpha = arg P_n(r omega)^3 = -pi/6 + arctan( sqrt3 r^{3n} / (2 + r^{3n}) ),  r^{3n} = e^{-xi (1-1/N)}
#   b = -arg L''(u*)/2 with L''(u*) = N^3 (rho_c + eps), rho_c = int tau^2 sum_j w(z_j) at (xi, t=N delta), |eps| <= 3V/N
#   |delta| <= 3 N TV_s / (N^3 rho_lo - 3 N^2 V),   TV_s = TV_tau( tau Im(omega x/(1-omega x)) ),  x = e^{-xi tau}
#   |E1| <= (N^3 |rho_c|_hi + 3 N^2 V) delta^2 / 2
# All N-dependent error terms decrease in N, so their N0 values bound every N >= N0.  Also T_lo for the far region.
from flint import arb, acb, ctx
import sys, math
ctx.prec=80
exec(open('centre_mid_uniform.py').read().split("def mid_uniform")[0])
def rho_complex(xi,tb):
    taus,h=cellballs(0,1,CELLS); tot=acb(0)
    for j in range(2):
        for tau in taus:
            z=om[j]*(-xi*tau+I*tau*tb).exp(); tot+=tau*tau*wfun(z)*h
    return tot
def clow_uniform(xi):
    N=arb(N0)
    lo0,V,cab=rho_data(xi,arb(0,0.05))
    rho_lo=arb(lo0.lower())
    # TV_s
    taus,h=cellballs(0,1,CELLS); TVs=arb(0)
    for tau in taus:
        x=(-xi*tau).exp(); z=om[0]*x
        s_=(z/(1-z)).imag
        ds=(( z/(1-z)**2)*(-xi)).imag          # d/dtau of z/(1-z) = z'/(1-z)^2, z' = -xi z
        TVs+=arb(abs(s_+tau*ds).upper())*h
    delta=3*N*TVs/(N**3*rho_lo-3*N**2*V)
    eta=N*delta
    assert eta.upper()<0.05
    E1=(N**3*cab+3*N**2*V)*delta*delta/2
    rc=rho_complex(xi,arb(0,eta.upper()))
    epsang=(3*V/(N*abs(rc).lower())).asin() if (3*V/(N*abs(rc))).upper()<1 else arb('inf')
    argL2=rc.arg()+arb(0,epsang.upper())
    b=-argL2/2
    R3n=(-xi*(1-arb(0,1)/N)).exp()   # covers 3n/N in [1-1/N0, 1]
    R3n=(-xi*(1-arb(0.5,0.5)/N)).exp()
    alpha=-PI/6+(arb(3).sqrt()*R3n/(2+R3n)).atan()
    phase=alpha-2*PI/3+b+arb(0,E1.upper())
    cosu=phase.cos()
    assert cosu.upper()<0, ('sign not certified',xi,cosu)
    taus2,h2=cellballs(0,1,CELLS)
    intg=sum((arb(gfun((-xi*tau).exp()).lower())*h2 for tau in taus2),arb(0))
    return float(abs(cosu).lower()), float(eta.upper()), float(E1.upper()), float(intg.lower())

def rest_uniform(xi):
    # far region for theta in [0, pi] minus |theta - 2pi/3| < 0.05 (the conjugate side is symmetric), relative to |P(r omega)|^3:
    #  (a) theta >= 0.05 (outside the peak arc): Dirichlet, G <= -(3/2)(T_lo - 2 g(r) kappa_max), kappa_max <= 1/sin(0.025) + 1/sin(0.075)
    #  (b) theta = t/N, t in [0, 40]: G <= -(3/2) sum g_k (2 cos k theta + 1) <= -N int g (1 + 2cos tau t) + 3(TV+S)
    #  (c) t in [40, 0.05 N]: Dirichlet with kappa <= (8/3.0) N / t * 1.01
    # Returns log of the max bound at N0 and asserts each exponent has slope < -3/(2 N0) (so N^{3/2} e^{G} decreases).
    N=arb(N0)
    taus,h=cellballs(0,1,CELLS)
    intg=sum((arb(gfun((-xi*tau).exp()).lower())*h for tau in taus),arb(0))
    Tlo=2*N/3*intg-arb(4)/3
    gr=arb(1)/3
    kap=1/arb('0.025').sin()+1/arb('0.075').sin()
    Ga=-arb(3)/2*(Tlo-2*gr*kap)
    assert (arb(3)/2*2/3*intg*N).lower()>arb(3)/2+ (3*gr*kap).upper(), 'far (a) slope'
    worst=Ga
    K=160; edges=[40*i/K for i in range(K+1)]
    for i in range(K):
        tb=arb((edges[i]+edges[i+1])/2,(edges[i+1]-edges[i])/2)
        lo=arb(0); tv=arb(0); sup=arb(0)
        for tau in taus:
            x=(-xi*tau).exp(); F=gfun(x)*(1+2*(tau*tb).cos())
            lo+=arb(F.lower())*h
            dg=-xi*x*(1-x*x)/(1+x*x+x)**2
            dF=dg*(1+2*(tau*tb).cos())-gfun(x)*2*(tau*tb).sin()*tb
            tv+=arb(abs(dF).upper())*h
            sup=max(sup,arb(abs(F).upper()),key=lambda q:q.upper())
        assert (lo*N).lower()>1.5, ('far (b) profile not positive',xi,edges[i],lo)
        Gb=-N*lo+3*(tv+sup)*2
        worst=max(worst,Gb,key=lambda q:q.upper())
    Gc=-arb(3)/2*(Tlo-2*gr*(arb(8)/3*N/40*arb('1.01')))
    worst=max(worst,Gc,key=lambda q:q.upper())
    return float(worst.upper())

if __name__=='__main__':
    lo,hi,K=float(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3])
    for i in range(K):
        a_=lo+(hi-lo)*i/K; b_=lo+(hi-lo)*(i+1)/K
        XB=arb((a_+b_)/2,(b_-a_)/2)
        c,e,E1,ig=clow_uniform(XB); gr_=rest_uniform(XB)
        print('xi in [%.4f,%.4f]: CLOWU %.5f  RESTLOG %.2f  eta<=%.2e E1<=%.2e intg>=%.5f'%(a_,b_,c,gr_,e,E1,ig),flush=True)
