# Laplace regime, UNIFORM in N >= N0 = 733: certificate over (v, eps) balls, v = N*W the scaled true saddle, eps = 1/N.
# See uniform_lib.py for the scaled closed forms.  Every piece below is shown non-increasing in N (asserted), so the
# value computed with eps in [0, 1/N0] bounds every N >= N0.  mu~ = (m - N - 1/4)/N^2 = -G'(v) is determined by v.
#   window |tau| <= c        : Taylor to order 4 + Cauchy order 5;  G3n ~ sqrt(eps), phi4 ~ eps, G5n ~ eps^{3/2}
#   zone 1 tau in [c, tau_s] : D <= -(tau^2/2) kappa,  kappa = (kap0 - eps kapb)/G2   (t_s = N-free)
#   zone 2 t in [t_s, t_C]   : cells, D <= N A + B with A_hi < -1/(2 N0)
#   zone 3 t >= t_C          : |F| <= E_n(e^-W), |P| <= e^{N a v/(v^2+t_C^2)} + 3;  slope < -3/(2 N0)
#   R(m)                     : sz_joint groups on the real axis;  slope < -7/(2 N0)
from flint import arb, acb, ctx
import sys, math
exec(open('uniform_lib.py').read())
N0=733
def PsiU(x): return 2*(-x).exp()/(1-(-x).exp())
def certify(V, E, cs=(2.0,2.5,3.0), dbg=False):
    Nb=max(733.0, 1.0/float(E.upper()))
    vlo=float(V.lower()); vmid=arb(V.mid())
    se_hi=arb(E.upper()).sqrt(); E_hi=arb(E.upper())
    s0,sb=parts_series(V,E,7)
    c0=[s0.coeffs()[j]*math.factorial(j) for j in range(7)]; cb=[sb.coeffs()[j]*math.factorial(j) for j in range(7)]
    G1=c0[1]+E*cb[1]; G2=c0[2]+E*cb[2]; G3=c0[3]+E*cb[3]; G4=c0[4]+E*cb[4]
    if not G2.lower()>0: return None,'G2'
    G2lo=arb(G2.lower()); G2hi=arb(G2.upper())
    mut=-G1                                   # mu~ for which v is the saddle
    if mut.upper()<0.999*float(E.lower()): return arb(0),'vacuous (m<2N)'
    if mut.lower()>0.0985: return arb(0),'vacuous (mu>=0.098, centre regime)'
    # Cauchy circle data (radius 0.8 vlo) for h = G - a/z (window) and for St, b separately (zone 1)
    rho=arb(0.6*vlo); K=96 if vlo>40 else (128 if vlo>2.9 else 512)
    G0c,bc,Stc=parts_acb(acb(vmid),E,vlo)
    hc=Stc+E*bc
    mxh=arb(0); mxS=arb(0); mxb=arb(0)
    for i in range(K):
        th=arb(2*math.pi*(i+0.5)/K, math.pi/K)
        z=acb(vmid)+rho*acb(th.cos(),th.sin())
        G0z,bz,Stz=parts_acb(z,E,0.4*vlo)
        hz=Stz+E*bz
        mxh=max(mxh,arb(abs(hz-hc).upper()),key=lambda q:q.upper())
        mxS=max(mxS,arb(abs(Stz-Stc).upper()),key=lambda q:q.upper())
        mxb=max(mxb,arb(abs(bz-bc).upper()),key=lambda q:q.upper())
    # base values at v
    G0v,bv,_=parts_acb(acb(V),E,vlo)
    St_v=G0v.real-a/V; b_v=bv.real
    best=None
    # ---- zone 1 constants
    t_s=max(0.1*vlo, 2.0/math.sqrt(Nb*float(G2lo.mid()))); ts=arb(t_s)
    d1=ts+arb(V.rad())
    M6_0=720*a/arb(vlo)**7+720*mxS*rho/(rho-d1)**7
    M6_b=720*mxb*rho/(rho-d1)**7
    kap0=arb(c0[2].lower())-abs(c0[4])*ts*ts/12-M6_0*ts**4/360
    kapb=abs(cb[2])+abs(cb[4])*ts*ts/12+M6_b*ts**4/360
    kappa=(kap0-E_hi*kapb)/G2hi
    if dbg: print('  t_s',t_s,'kap0',kap0,'kapb',kapb,'mxS',mxS,'mxb',mxb,'mxh',mxh,'G2',G2)
    if not kappa.lower()>0: return None,'kappa'
    # ---- zone 2 cells t in [t_s, t_C]
    small_v = vlo < 8
    t_C = 10*vlo if small_v else 3*vlo
    width=min(0.25*max(1.0,vlo/8.0)/max(1.0,160/max(vlo,1.0)/8), t_s/4)
    Z2=arb(0); tt=t_s
    while tt<t_C:
        w_=min(width,t_C-tt)
        z=acb(V,arb(tt+w_/2,w_/2))
        G0z,bz,_=parts_acb(z,E,vlo)
        A=A_diff(V,arb(tt+w_/2,w_/2),E); B=bz.real-b_v
        if not A.upper()<-1.0/(2*Nb): return None,('zone2 slope',round(tt,3),float(A.upper()))
        eps_x=A-a*V/(V*V+arb(tt)**2)
        val=arb((Nb*A+B).upper()).exp()*w_+3*arb((Nb*eps_x+B).upper()).exp()*w_
        Z2+=2*val
        tt+=w_
        # adaptive widening once the cell exponent is hopelessly negative (each cell enclosure stays rigorous)
        if (Nb*A+B).upper()<-80: width=min(width*1.5, max(0.25, 0.05*tt))
    Z2=Z2*arb(Nb).sqrt()*G2hi.sqrt()
    # ---- zone 3 (t >= t_C): crude (valid for all eps) or split (needs eps_lo > 0)
    tC=arb(t_C); Z3=None
    beta=PsiU(V)/V+a*V/(V*V+tC*tC)-a/V-St_v
    if beta.upper()<-3.0/(2*Nb):
        X=Nb*beta-b_v+V-6*(1-(-V).exp()).log()+arb(4).log()
        Z3=arb(Nb)**arb(1.5)*G2hi.sqrt()*PI*X.exp()
    elo=float(E.lower())
    if elo>0:
        Vl=float(V.lower())
        Qmax=sum((3*math.sqrt(3)/(2*r))*math.exp(-(r-1)*Vl)*2*(1/(3*r*Vl*elo)+1) for r in range(1,200) if r%3)
        Smax=sum((3.0/r)*math.exp(-r*Vl)*2*(1/(3*r*Vl)+float(E.upper())) for r in range(1,200))
        Dmax=Qmax*math.exp(-Vl)
        Av=a*V/(V*V+tC*tC)-a/V+2*arb(Smax)
        if Av.upper()<-5.0/(2*Nb):
            Qv=arb(math.exp(b_v.lower()))  # Qt(v)*sinc(v) = e^{b_v}
            X2=Nb*Av+arb(Qmax).log()+arb(Dmax)+arb(Dmax)**2-b_v+arb(4).log()
            # N-monotone: log Qmax grows like log N (1/eps_lo fixed per ball), slope condition -5/(2Nb) covers N^{3/2} and log
            Z3s=arb(Nb)**arb(1.5)*G2hi.sqrt()*PI*X2.exp()
            Z3=Z3s if Z3 is None or Z3s.upper()<Z3.upper() else Z3
    if Z3 is None: return None,('zone3',float(beta.upper()))
    # ---- R(m) groups
    def grp(ak,lead):
        bestv=None
        for tp in (0.15,0.25,0.35,0.5,0.7,1.0,1.4,2.0,3.0):
            vp=arb(tp)*V
            br=mut*vp+ak/vp+PsiU(vp)/vp-a/V-St_v-mut*V
            if not br.upper()<-7.0/(2*Nb): continue
            val=(lead+Nb*br+vp-V-b_v-6*(1-(-vp).exp()).log()+arb(abs(float(vp.upper()-float(V.lower())))/(4*Nb))).upper()
            bestv=val if bestv is None or val<bestv else bestv
        return None if bestv is None else arb(bestv)
    lg=[grp(a/4,arb(6).log()),grp(a/9,arb(18).log()),grp(a/16,arb(12).log()),grp(a/25,(3*(4*PI*arb('0.2')*Nb*Nb+4*(PI*arb('0.2')).sqrt()*Nb+4)).log()),grp(arb(0),arb('50.1501').log())]
    if any(x is None for x in lg): return None,'Rterm slope'
    Rrel=sum((x.exp() for x in lg),arb(0))*arb(Nb)**arb(1.5)*G2hi.sqrt()*PI
    for c in cs:
        cc=arb(c)
        G3n=se_hi*abs(G3)/G2lo**arb(1.5)
        ph4=E*G4/(G2*G2); ph4_lo=arb(ph4.lower()); ph4_hi=arb(ph4.upper())
        d=cc*se_hi/G2lo.sqrt()+arb(V.rad())
        if not d.upper()<0.7*rho.lower(): continue
        M5=120*a/arb(vlo)**6+120*mxh*rho/(rho-d)**6
        G5n=E_hi**arb(1.5)*M5/G2lo**arb(2.5)
        cells=400; Lint=arb(0)
        for i in range(cells):
            t1=arb(c*i/cells); t2=arb(c*(i+1)/cells)
            th3=G3n*t2**3/6+G5n*t2**5/120
            cosl=1-th3*th3/2
            ex_lo=ph4_lo*(t1**4 if ph4_lo.lower()>=0 else t2**4)/24-G5n*t2**5/120
            if cosl.lower()>0: Lc=(-(t2*t2)/2).exp()*cosl*ex_lo.exp()
            else: Lc=(-(t1*t1)/2).exp()*cosl*(ph4_hi*t2**4/24+G5n*t2**5/120).exp()
            Lint+=arb(Lc.lower())*(t2-t1)
        Lint=2*Lint
        if not Lint.lower()>0: continue
        Z1=2*(PI/(2*kappa)).sqrt()*(cc*(kappa/2).sqrt()).erfc()
        weps=3*2*cc*(ph4_hi*cc**4/24+G5n*cc**5/120).exp()*(-(Nb*a*V/(V*V+cc*cc/(Nb*G2lo)))).exp()
        B=(Z1+Z2+Z3+Rrel+weps)/Lint
        if dbg: print('  c',c,'Z1',float(Z1.upper()),'Z2',float(Z2.upper()),'Z3',float(Z3.upper()),'R',float(Rrel.upper()),'Lint',float(Lint.lower()),'B',float(B.upper()))
        if best is None or B.upper()<best.upper(): best=B
    return best,'ok'
def mu_upper(V, E):
    # rigorous upper bound on mu~ = -G'(v) for v in V, eps in E, WITHOUT series (valid for any eps ball):
    #   -G'(v) = a/v^2 - h'(v),  h = St + eps*b holomorphic on |z - v| <= rho;  |h'(v)| <= max_{circle}|h - h(v)| / (rho - rad)
    vlo=float(V.lower()); vmid=arb(V.mid()); rho=arb(0.5*vlo); K=64
    G0c,bc,Stc=parts_acb(acb(V),E,vlo); hc=Stc+E*bc
    mx=arb(0)
    for i in range(K):
        th=arb(2*math.pi*(i+0.5)/K, math.pi/K)
        z=acb(vmid)+rho*acb(th.cos(),th.sin())
        G0z,bz,Stz=parts_acb(z,E,0.45*vlo)
        mx=max(mx,arb(abs(Stz+E*bz-hc).upper()),key=lambda q:q.upper())
    return a/arb(vlo)**2+mx/(rho-arb(V.rad()))
if __name__=='__main__':
    vlo,vhi,K,KE=float(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])
    edges=[vlo*(vhi/vlo)**(j/K) for j in range(K+1)]
    pass
    worst=0; fails=[]
    for j in range(K):
        V=arb((edges[j]+edges[j+1])/2,(edges[j+1]-edges[j])/2)
        ec=min(1.0/N0, 1.0966/(3*edges[j+1]**2))
        eedges=[0.0]+[ec*(N0*ec)**(-i/KE) for i in range(KE+1)] if ec<1.0/N0 else [0.0]+[(1.0/N0)*2.0**(-(KE-i)) for i in range(KE+1)]
        eedges=sorted(set(eedges))
        for i in range(len(eedges)-1):
            E=arb((eedges[i]+eedges[i+1])/2,(eedges[i+1]-eedges[i])/2)
            try: B,st=certify(V,E)
            except AssertionError as e: B,st=None,'assert:'+str(e)
            if B is None or not B.upper()<1: fails.append((round(edges[j],4),i,st if B is None else float(B.upper())))
            else: worst=max(worst,float(B.upper()))
    print('UNIFORM-LAPLACE v in [%.3f,%.3f] (%d x %d balls): failures %d %s ; worst B = %.4f'%(vlo,vhi,K,KE,len(fails),fails[:4],worst),flush=True)
