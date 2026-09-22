# Classes 0,1,2, Laplace region: V-tile certificate (derivation: GAP_CLASS012.md, section B).
# For a V-tile [Va,Vb] and ALL X >= Xmin = max(1455/Vb, 4 Va/a5)  (i.e. n >= 291 and t >= 4), certify
#     sgn(a(cls)) * M  -  4*delta  -  E2  > 0,
#   M      = sum_h A_h z5^{h cls} e^{Omega_h(W)}        (main amplitude; |e^{Omega_h}| = 1; M is real)
#   delta  = eps_win + erfc(cc/sqrt2) + eps_mid          (relative error of each peak integral, h = 1..4)
#   E2     = 4(a5 F + T_crude) + T_e + T_rho             (head sum, e(m), off-peak kernel, R-term; relative)
# Every piece is either X-independent or non-increasing in X, so evaluating at Xmin covers all X >= Xmin.
import sys, math
exec(open('cls012_lib.py').read())
AH={1:acb((PI/5).cos(),(PI/5).sin()),4:acb((PI/5).cos(),-(PI/5).sin()),2:acb(1),3:acb(1)}
def z5(k):
    a=2*PI*(k%5)/5; return acb(a.cos(),a.sin())
B1={1:arb(-3)/10,2:arb(-1)/10,3:arb(1)/10,4:arb(3)/10}
ECON=arb('27.4636')
SGN={0:1,1:-1,2:-1}
def lo(x): return arb(x.lower()) if isinstance(x,arb) else arb(x.real.lower())
def up(x): return arb(x.upper()) if isinstance(x,arb) else arb(abs(x).upper())
def amax(*xs): return max(xs, key=lambda z: float(z.lower()))
def Mval(xi,cls):
    M=acb(0)
    for h in (1,2,3,4):
        Om=acb(0)
        for c in (1,2,3,4): Om+=B1[c]*(1-z5(-h*c)*xi).log()
        M+=AH[h]*z5(h*cls)*Om.exp()
    return M
def crude(Xmin,Xi2,kap,L1,eta1,Y2):
    X=Xmin; W=1/X
    Fm=(2*PI*Xi2).sqrt()*X**arb(1.5)*(-X*kap+4*L1).exp()
    Tcr=Fm/(2*PI)*(2*eta1*(W+a5)+2*Y2*((a5*X/(1+eta1**2)).exp()+1+a5*X/eta1)+a5**2/Y2*(a5/Y2).exp())
    Te=Fm*(W/6).exp()
    K=(4*PI*(a5*X*X+1)).sqrt()+2; K5=K/5
    Trho=Fm*(W/6).exp()*(W*(4*(a5*X/4).exp()+(a5*X/9).exp()*2*K5*(K5+1))+ECON)
    cmin=kap-a5/(1+eta1**2)          # every crude term is C X^p e^{-c X}, p <= 4.5, c >= cmin: decreasing for X > 4.5/cmin
    cmin=arb.min(cmin, kap-a5/4)     # AUDIT FIX: the Rankin term Trho carries e^{a5 X/4}, so its rate is kap - a5/4
    assert cmin.lower()>0 and (Xmin-arb('4.5')/cmin).lower()>0
    return 4*(a5*Fm+Tcr)+Te+Trho
def window(cc,k3,k4,k1,Rb,S=800):
    ew=arb(0)
    for i in range(S):
        sa=cc*i/S; sb=cc*(i+1)/S
        r=(k4*sb**4+k1*sb+Rb).exp()
        ew+=(-sa*sa/2).exp()*(min(k3*sb**3,arb(2),key=lambda z: float(z.upper()))*r+r-1)*(sb-sa)
    return 2*ew/(2*PI).sqrt()
GRID=[i*0.005 for i in range(0,41)]+[0.2+0.01*i for i in range(1,81)]+[1.0+0.02*i for i in range(1,51)]
def tile(Va,Vb,cc=arb(4),eta1=arb(2),Y2=arb(1)):
    Va=arb(Va); Vb=arb(Vb)
    V=arb((Va+Vb)/2,(Vb-Va)/2)
    Xmin=arb(amax(arb(1455)/Vb, 4*Va/a5).lower()); Wmax=1/Xmin
    xi=(-V).exp(); L2=Li2(xi).real; L1=-(1-xi).log()
    al=alpha(V); g0=gam(V).real
    assert al.lower()>0 and (a5-al).lower()>0
    c0=Xi_series(V,arb(0))
    Xi0=c0[0].real; Xi2=-2*c0[2].real
    assert Xi2.lower()>0
    Xi2l=lo(Xi2)
    c1r=up(c0[1]); c3=abs(c0[3].imag)+up(c0[3].real)      # c1 = 0 and c3 is i*real exactly; keep ball slack
    eta0=cc/(Xmin*Xi2l).sqrt()                            # largest window half-width in eta
    X4=[]                                                 # sup|Xi''''| on sub-balls of width 0.01 (even in eta)
    for j in range(100):
        cj=Xi_series(V,arb((j+0.5)*0.01,0.005)); X4.append(float(up(24*cj[4]).upper()))
    def X4sup(b): return arb(max(X4[:max(1,int(math.ceil(b/0.01-1e-9)))]))
    Xi4=X4sup(float(eta0.upper()))
    zmax=(1+eta1*eta1).sqrt()
    Rbar=5*Wmax/3*zmax*xi*(1/(1-xi)+zmax/(1-xi)**2)      # EM remainder of log E_n on |eta| <= eta1
    beta1=arb('0.8')*V*xi/(1-xi); wbar=arb('0.8')*L1    # |d Omega/d eta| <= beta1, |Omega| <= wbar
    # window s = eta sqrt(X Xi2) in [0,cc]: |e^psi - 1| <= min(k3 s^3,2) e^{r} + e^{r} - 1, r = k4 s^4 + k1 s + Rbar
    k3=up(c3/(Xi2l**arb(1.5)*Xmin.sqrt())); k4=up(Xi4/(24*Xi2l**2*Xmin))
    # c1 = Xi'(0) = i(alpha - a5 + gamma(V) - V gamma'(V)) = 0 IDENTICALLY by the definition of alpha (the V-ball
    # enclosure of c1 is wide only through dependency); sanity check at the tile midpoint:
    assert abs(Xi_series(arb(V.mid()),arb(0))[1]).upper()<1e-25
    k1=up(beta1/(Xmin*Xi2l).sqrt()); Rb=up(Rbar)
    eps_win=window(cc,k3,k4,k1,Rb)
    gtail=(cc/arb(2).sqrt()).erfc()
    # mid region: D(eta) = Re(Xi(0)-Xi(eta)) >= d_k eta^2 on each tile
    mid=arb(0); dmin=None
    for a,b in zip(GRID[:-1],GRID[1:]):
        ea=arb(a); eb=arb(b); cands=[]
        if b<=1.0: cands.append(Xi2l/2-X4sup(b)*eb*eb/24)     # Taylor (odd terms are imaginary)
        if a>0:
            e_=arb((ea+eb)/2,(eb-ea)/2)
            cands.append((Xi0-Xi_val(V,e_).real)/(e_*e_))
        dk=amax(*[lo(c) for c in cands])
        assert dk.lower()>0, (a,b,[str(c) for c in cands])
        dmin=dk if dmin is None or dk.lower()<dmin.lower() else dmin
        arg=amax(cc*(dk/Xi2).sqrt(), ea*(Xmin*dk).sqrt())
        mid+=(Xi2/(2*dk)).sqrt()*lo(arg).erfc()
    eps_mid=(2*wbar+Rb).exp()*mid
    delta=eps_win+gtail+eps_mid
    kap=a5-g0-arb(4)/5*L2
    E2=crude(Xmin,Xi2,kap,L1,eta1,Y2)
    out={}
    for cls in (0,1,2):
        Mr=SGN[cls]*Mval(xi,cls).real
        out[cls]=(float((Mr-4*delta-E2).lower()), float(Mr.lower()))
    info=dict(Xmin=round(float(Xmin.mid()),1),Xi2=float(Xi2.mid()),Xi4=float(Xi4.mid()),c3=float(c3.mid()),eta0=float(eta0.mid()),
              eps_win=float(eps_win.upper()),eps_mid=float(eps_mid.upper()),dmin=float(dmin.lower()),E2=float(E2.upper()),delta=float(delta.upper()))
    return out,info

def tailcell(VT,cc=arb(4),eta1=arb(2),Y2=arb(1)):
    """All V >= VT (VT >= 10): uniform analytic bounds, Cauchy estimates on the strip |Im eta| <= 1/2."""
    VT=arb(VT); xiT=(-VT).exp(); L1=-(1-xiT).log(); L2=Li2(xiT).real
    epsg=arb(12)/25*Li2((-VT/2).exp()).real           # sup |gamma(V zeta)/zeta| on |Im eta| <= 1/2
    VQ=VT*xiT/(1-xiT)                                  # V Q(V) <= V xi/(1-xi), decreasing for V >= 1
    Xmin=4*VT/a5; Wmax=1/Xmin                          # X >= 4V/alpha >= 4V/a5 >= 4VT/a5 (>= 1455/V for V >= 9.78)
    Xi2=arb(2*a5,8*epsg); Xi2l=lo(Xi2)
    c3=a5+8*epsg; Xi4=24*a5+384*epsg
    zmax=(1+eta1*eta1).sqrt()
    Rbar=5*Wmax/3*zmax*xiT*(1/(1-xiT)+zmax/(1-xiT)**2)
    beta1=arb('0.8')*VT*xiT/(1-xiT); wbar=arb('0.8')*L1
    k3=up(c3/(Xi2l**arb(1.5)*Xmin.sqrt())); k4=up(Xi4/(24*Xi2l**2*Xmin)); k1=up(beta1/(Xmin*Xi2l).sqrt()); Rb=up(Rbar)
    eps_win=window(cc,k3,k4,k1,Rb)
    gtail=(cc/arb(2).sqrt()).erfc()
    mid=arb(0); dmin=None
    for a,b in zip(GRID[:-1],GRID[1:]):
        ea=arb(a); eb=arb(b); cands=[Xi2l/2-Xi4*eb*eb/24]
        if a>0: cands.append((a5*ea*ea/(1+eb*eb)-arb(12)/25*L2)/(eb*eb))
        dk=amax(*[lo(c) for c in cands]); assert dk.lower()>0,(a,b)
        dmin=dk if dmin is None or dk.lower()<dmin.lower() else dmin
        arg=amax(cc*(dk/up(Xi2)).sqrt(), ea*(Xmin*dk).sqrt())
        mid+=(up(Xi2)/(2*dk)).sqrt()*lo(arg).erfc()
    eps_mid=(2*wbar+Rb).exp()*mid
    delta=eps_win+gtail+eps_mid
    kap=a5-L2                                          # a5 - gamma - (4/5)Li2 >= a5 - Li2(xi)
    E2=crude(Xmin,up(Xi2),kap,L1,eta1,Y2)
    avals={0:(5+arb(5).sqrt())/2,1:arb(5).sqrt(),2:(5-arb(5).sqrt())/2}
    out={}
    for cls in (0,1,2):
        Mr=avals[cls]-4*(wbar.exp()-1)
        out[cls]=(float((Mr-4*delta-E2).lower()),float(Mr.lower()))
    return out,dict(Xmin=float(Xmin.mid()),eps_win=float(eps_win.upper()),eps_mid=float(eps_mid.upper()),E2=float(E2.upper()),delta=float(delta.upper()),dmin=float(dmin.lower()))
if __name__=='__main__':
    if sys.argv[1]=='tail':
        print('tail cell V >=',sys.argv[2],':',tailcell(float(sys.argv[2])),flush=True); sys.exit()
    Va=float(sys.argv[1]); Vb=float(sys.argv[2]); w=float(sys.argv[3])
    worst={0:9,1:9,2:9}; v=Va; cnt=0; wi=None; fails=0
    while v<Vb-1e-12:
        v2=min(v+w,Vb)
        if Vb-v2<1e-9: v2=Vb     # AUDIT (referee): snap the last tile to Vb -- float drift left (11.9999999999998, 12) uncovered
        o,info=tile(v,v2); cnt+=1
        for c in o: worst[c]=min(worst[c],o[c][0])
        if not all(o[c][0] > 0 for c in o): fails+=1; print('FAIL',v,v2,o,info,flush=True)   # NaN-safe per audit
        if wi is None or info['delta']>wi[0]: wi=(info['delta'],v,info)
        v=v2
    print(f'V in [{Va},{Vb}] width {w}: {cnt} tiles, {fails} failures; min certified margin sgn(a)M - 4 delta - E2 by class: '
          + ', '.join(f'{c}: {worst[c]:.4f}' for c in worst)+f'; worst-delta tile V={wi[1]:.3f}: {wi[2]}',flush=True)
