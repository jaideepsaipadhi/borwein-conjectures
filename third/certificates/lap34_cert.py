# lap34_cert.py -- classes 3,4, Laplace region (t >= 4, mu < 0.105), n >= NMIN: certify c(m) < 0.
# Derivation: GAP_LAP34.md.  Exact representation (cls012 B.1): c(m) = e(m) + sum_h A_h z5^{hm}(a5 H_h + K_h) + rho.
# Saddle (re-centred, includes the thin factor's -5n u):  alpha(V) X^2 = m - 1/6 - 5n,  X = 1/W = 5n/V.
# Then  c(m) = Amp [ g(V) + Delta ],  Amp = W^{3/2} e^{X Xi_V(0)}/sqrt(2 pi Xi2) > 0,  and we certify
#     -g(V) - |Delta| > 0,   |Delta| <= eps_win + |g| erfc(cc/sqrt2) + eps_mid + E2,
# for a V-tile [Va,Vb] and ALL X >= Xmin = max(5 NMIN/Vb, 3 Va/a5); every term is non-increasing in X.
import sys, math
def _fin(v):
    v=float(v); assert math.isfinite(v), 'non-finite enclosure (AUDIT: NaN guard)'; return v
exec(open('lap34_lib.py').read())
NMIN=821
def crude(Xmin,Xi2,kap,L1,eta1,Y2,eV,vx,erate):
    """Head sum a5 H_h, off-window/far kernel, e(m), R-term (Rankin), all relative to Amp, summed over h.
       eV = upper bound for e^{V} (tile) or 1 with erate = a5/3 (tail, uses V <= a5 X/3).
       vx = upper bound for V (enters m - 1/6 = alpha X^2 + V X <= a5 X^2 + vx X)."""
    X=Xmin; W=1/X
    Fm=eV*(2*PI*Xi2).sqrt()*X**arb(1.5)*(-X*(kap-erate)+4*L1).exp()
    Tcr=Fm/(2*PI)*(2*eta1*(W+a5)+2*Y2*((a5*X/(1+eta1**2)).exp()+1+a5*X/eta1)+a5**2/Y2*(a5/Y2).exp())
    Te=Fm*(W/6).exp()
    K=(4*PI*(a5*X*X+vx*X+1)).sqrt()+2; K5=K/5
    Trho=Fm*(W/6).exp()*(W*(4*(a5*X/4).exp()+(a5*X/9).exp()*2*K5*(K5+1))+ECON)
    # every crude term is C X^p e^{-cX} (or a sum of such with log-derivative <= p/X - c), p <= 4.5,
    # c >= cmin = kap - erate - max(a5/(1+eta1^2), a5/4, a5/9): decreasing for X > 4.5/cmin
    cmin=kap-erate-amax(a5/(1+eta1**2),a5/4)
    assert cmin.lower()>0 and (Xmin-arb('4.5')/cmin).lower()>0
    return 4*(a5*Fm+Tcr)+Te+Trho
def epsR(W,xi,zmax,wbar,eV):
    """sup |G - g(V zeta)| where G carries the Euler-Maclaurin remainders R_h (cls012 B.2):
       <= e^V sum_h |e^{Omega_h}| (e^{Rbar}-1) <= 4 e^{wbar} e^V (e^{Rbar}-1)."""
    Rb=5*W/3*zmax*xi*(1/(1-xi)+zmax/(1-xi)**2)
    return 4*wbar.exp()*eV*(Rb.exp()-1)
def window(cc,k3,k4,k1g,eR,gabs,S=800):
    """(2/sqrt(2pi)) int_0^cc e^{-s^2/2} [ e^{k4 s^4}(k1g s + eR) + |g|(min(k3 s^3,2) e^{k4 s^4} + e^{k4 s^4} - 1) ] ds
       (integrand increasing in s: upper Riemann sum)."""
    ew=arb(0)
    for i in range(S):
        sa=cc*i/S; sb=cc*(i+1)/S
        r=(k4*sb**4).exp()
        ew+=(-sa*sa/2).exp()*(r*(k1g*sb+eR)+gabs*(min(k3*sb**3,arb(2),key=lambda z: float(z.upper()))*r+r-1))*(sb-sa)
    return 2*ew/(2*PI).sqrt()
GRID=[i*0.005 for i in range(0,41)]+[0.2+0.01*i for i in range(1,81)]+[1.0+0.02*i for i in range(1,51)]
def gsup(V,a,b,cls,sub=4):
    """sup |g(V(1+i eta))| for eta in [a,b] (and by conjugate symmetry [-b,-a]); V a ball."""
    best=0.0
    for j in range(sub):
        ea=a+(b-a)*j/sub; eb=a+(b-a)*(j+1)/sub
        e_=arb((ea+eb)/2,(eb-ea)/2)
        best=max(best,_fin(abs(gfun(acb(V,V*e_),cls)).upper()))
    return arb(best)
def gpsup(V,e0,cls,sub):
    best=0.0
    for j in range(sub):
        ea=e0*j/sub; eb=e0*(j+1)/sub
        e_=arb((ea+eb)/2,(eb-ea)/2)
        best=max(best,_fin(abs(gder(acb(V,V*e_),cls)).upper()))
    return arb(best)
def epsgf(Vp):
    """|g(s) - finf| <= epsg(Re s) for every complex s with Re s = Vp > 0; decreasing in Vp (series in xi with
       nonnegative coefficients, divided by xi).  finf = -1 (checked in tailcell / lap34_lib.finf)."""
    Vp=arb(Vp); x=(-Vp).exp(); L1=-(1-x).log(); wb=arb('0.8')*L1
    return 4*Vp.exp()*(arb('0.8')*(L1-x)+wb**2*wb.exp()/2)
def tile(Va,Vb,cc=arb(4),eta1=arb(2),Y2=arb(1),nmin=NMIN):
    Va=arb(Va); Vb=arb(Vb)
    V=arb((Va+Vb)/2,(Vb-Va)/2)
    Xmin=arb(amax(arb(5*nmin)/Vb, 3*Va/a5).lower()); Wmax=1/Xmin
    xi=(-V).exp(); L2=Li2(xi).real; L1=-(1-xi).log(); wbar=arb('0.8')*L1; eV=Vb.exp()
    al=alpha(V); g0=gam(V).real
    assert al.lower()>0 and (a5-al).lower()>0
    c0=Xi_series(V,arb(0))
    Xi0=c0[0].real; Xi2=-2*c0[2].real
    assert Xi2.lower()>0
    Xi2l=lo(Xi2)
    c3=abs(c0[3].imag)+up(c0[3].real)
    assert abs(Xi_series(arb(V.mid()),arb(0))[1]).upper()<1e-25
    eta0=cc/(Xmin*Xi2l).sqrt()
    X4=[]
    for j in range(100):
        cj=Xi_series(V,arb((j+0.5)*0.01,0.005)); X4.append(_fin(up(24*cj[4]).upper()))
    def X4sup(b): return arb(max(X4[:max(1,int(math.ceil(b/0.01-1e-9)))]))
    Xi4=X4sup(float(eta0.upper()))
    k3=up(c3/(Xi2l**arb(1.5)*Xmin.sqrt())); k4=up(Xi4/(24*Xi2l**2*Xmin))
    out={}
    for cls in (3,4):
        gV=gfun(V,cls); assert abs(gfun(arb(V.mid()),cls).imag).upper()<1e-25   # g(V) real for real V (h <-> 5-h conjugate pairs)
        gV=gV.real; gabs=arb(abs(gV).upper())
        e0=float(eta0.upper())
        g1=gpsup(V,e0,cls,sub=max(4,int(math.ceil(float(Vb.upper())*e0/0.05))))
        g1=amax(-g1,-epsgf(Va-1))*(-1)          # min(ball bound, Cauchy bound eps_g(V-1)/1 on the unit circle)
        k1g=up(Vb*g1/(Xmin*Xi2l).sqrt())
        zw=(1+eta0**2).sqrt()
        eRw=up(epsR(Wmax,xi,zw,wbar,eV))
        eps_win=window(cc,k3,k4,k1g,eRw,gabs)
        gtail=gabs*(cc/arb(2).sqrt()).erfc()
        mid=arb(0); dmin=None
        for a,b in zip(GRID[:-1],GRID[1:]):
            ea=arb(a); eb=arb(b); cands=[]
            if b<=1.0: cands.append(Xi2l/2-X4sup(b)*eb*eb/24)
            if a>0:
                e_=arb((ea+eb)/2,(eb-ea)/2)
                cands.append((Xi0-Xi_val(V,e_).real)/(e_*e_))
            dk=amax(*[lo(c) for c in cands])
            assert dk.lower()>0, (a,b)
            dmin=dk if dmin is None or dk.lower()<dmin.lower() else dmin
            arg=amax(cc*(dk/Xi2).sqrt(), ea*(Xmin*dk).sqrt())
            if float(lo(arg).lower())>40: continue          # skipped here, bounded jointly by TAILSKIP below
            sG=gsup(V,a,b,cls)+up(epsR(Wmax,xi,(1+eb*eb).sqrt(),wbar,eV))
            mid+=sG*(up(Xi2)/(2*dk)).sqrt()*lo(arg).erfc()
        # skipped tiles (arg > 40): bounded jointly below by TAILSKIP
        mid+=TAILSKIP(V,cls,Xmin,Xi2,dmin,cc,eta1,eV,wbar,xi,Wmax)
        kap=a5-g0-arb(4)/5*L2
        E2=crude(Xmin,Xi2,kap,L1,eta1,Y2,eV,Vb,arb(0))
        D=eps_win+gtail+mid+E2
        marg=-gV-D
        out[cls]=(float(marg.lower()),float((-gV).lower()),dict(win=float(eps_win.upper()),gt=float(gtail.upper()),mid=float(mid.upper()),E2=float(E2.upper()),k3=float(k3.upper()),k1g=float(k1g.upper()),eR=float(eRw.upper())))
    return out,dict(Xmin=float(Xmin.mid()),eta0=float(eta0.mid()),dmin=float(dmin.lower()))
def TAILSKIP(V,cls,Xmin,Xi2,dmin,cc,eta1,eV,wbar,xi,Wmax):
    """Tiles whose erfc argument exceeds 40 were skipped.  Their total is <= sup_{|eta|<=eta1}|G| * sum over those
       tiles of sqrt(Xi2/(2dk)) erfc(40) <= Gbar * 170 * sqrt(Xi2/(2 dmin)) * erfc(40)  (<= 170 tiles in GRID).
       Gbar <= e^V * 4 e^{wbar} * e^{Rbar}: crude modulus bound of G (no cancellation used)."""
    Rb=5*Wmax/3*(1+eta1**2).sqrt()*xi*(1/(1-xi)+(1+eta1**2).sqrt()/(1-xi)**2)
    Gbar=eV*4*wbar.exp()*Rb.exp()
    return Gbar*170*(up(Xi2)/(2*dmin)).sqrt()*arb(40).erfc()
def tailcell(VT,cc=arb(4),eta1=arb(2),Y2=arb(1),nmin=NMIN):
    """All V >= VT: X >= max(5 nmin/V, 3V/a5) >= X* = sqrt(15 nmin/a5), and X >= 3 VT/a5; V <= a5 X/3."""
    VT=arb(VT); xiT=(-VT).exp(); L1=-(1-xiT).log(); L2=Li2(xiT).real; wbar=arb('0.8')*L1
    epsg_xi=arb(12)/25*Li2((-VT/2).exp()).real        # Cauchy bound on gamma(V zeta)/zeta, strip |Im eta| <= 1/2
    Xmin=arb(amax((15*arb(nmin)/a5).sqrt(),3*VT/a5).lower()); Wmax=1/Xmin
    Xi2=arb(2*a5,8*epsg_xi); Xi2l=lo(Xi2)
    c3=a5+8*epsg_xi; Xi4=24*a5+384*epsg_xi
    k3=up(c3/(Xi2l**arb(1.5)*Xmin.sqrt())); k4=up(Xi4/(24*Xi2l**2*Xmin))
    # thin function: |g(s) - finf| <= epsg for Re s = V >= VT (finf = -1 exactly, checked), decreasing in V
    epsg=epsgf(VT)
    eV=VT.exp()   # e^V (e^{Rbar}-1) is decreasing in V (Rbar/xi fixed shape): evaluate at VT
    out={}
    for cls in (3,4):
        fi=finf(cls); assert abs(fi+1).upper()<1e-30
        eta0=cc/(Xmin*Xi2l).sqrt()
        eRw=up(epsR(Wmax,xiT,(1+eta0**2).sqrt(),wbar,eV))
        eps_win=window(cc,k3,k4,arb(0),2*epsg+eRw,1+epsg)
        gtail=(1+epsg)*(cc/arb(2).sqrt()).erfc()
        mid=arb(0); dmin=None
        for a,b in zip(GRID[:-1],GRID[1:]):
            ea=arb(a); eb=arb(b); cands=[Xi2l/2-Xi4*eb*eb/24]
            if a>0: cands.append((a5*ea*ea/(1+eb*eb)-arb(12)/25*L2)/(eb*eb))
            dk=amax(*[lo(c) for c in cands]); assert dk.lower()>0,(a,b)
            dmin=dk if dmin is None or dk.lower()<dmin.lower() else dmin
            arg=amax(cc*(dk/up(Xi2)).sqrt(), ea*(Xmin*dk).sqrt())
            sG=1+epsg+up(epsR(Wmax,xiT,(1+eb*eb).sqrt(),wbar,eV))
            mid+=sG*(up(Xi2)/(2*dk)).sqrt()*lo(arg).erfc()
        kap=a5-L2
        E2=crude_tail(Xmin,up(Xi2),kap,L1,eta1,Y2)
        D=eps_win+gtail+mid+E2
        marg=(1-epsg)-D
        out[cls]=(float(marg.lower()),dict(epsg=float(epsg.upper()),win=float(eps_win.upper()),gt=float(gtail.upper()),mid=float(mid.upper()),E2=float(E2.upper()),eR=float(eRw.upper())))
    return out,dict(Xmin=float(Xmin.mid()),dmin=float(dmin.lower()))
def crude_tail(Xmin,Xi2,kap,L1,eta1,Y2):
    """crude terms for V >= VT using e^V <= e^{a5 X/3} and V X <= a5 X^2/3."""
    X=Xmin; W=1/X
    Fm=(2*PI*Xi2).sqrt()*X**arb(1.5)*(-X*(kap-a5/3)+4*L1).exp()
    Tcr=Fm/(2*PI)*(2*eta1*(W+a5)+2*Y2*((a5*X/(1+eta1**2)).exp()+1+a5*X/eta1)+a5**2/Y2*(a5/Y2).exp())
    Te=Fm*(W/6).exp()
    K=(4*PI*(a5*X*X*4/3+1)).sqrt()+2; K5=K/5
    Trho=Fm*(W/6).exp()*(W*(4*(a5*X/4).exp()+(a5*X/9).exp()*2*K5*(K5+1))+ECON)
    cmin=kap-a5/3-a5/4
    assert cmin.lower()>0 and (Xmin-arb('4.5')/cmin).lower()>0
    return 4*(a5*Fm+Tcr)+Te+Trho
if __name__=='__main__':
    if sys.argv[1]=='tail':
        print('tail cell V >=',sys.argv[2],':',tailcell(float(sys.argv[2])),flush=True); sys.exit()
    Va=float(sys.argv[1]); Vb=float(sys.argv[2]); w=float(sys.argv[3])
    worst={3:(9,None),4:(9,None)}; v=Va; cnt=0; fails=0
    while v<Vb-1e-12:
        v2=min(v+w,Vb)
        if Vb-v2<1e-9: v2=Vb   # AUDIT: float drift left (11.99999999999983, 12) uncovered
        o,info=tile(v,v2); cnt+=1
        for c in o:
            if not (o[c][0]>=worst[c][0]): worst[c]=(o[c][0],(v,v2,o[c],info))
        if not all(o[c][0]>0 for c in o): fails+=1; print('FAIL',v,v2,o,info,flush=True)   # AUDIT: NaN-safe
        v=v2
    print(f'V in [{Va},{Vb}] width {w}: {cnt} tiles, {fails} failures', flush=True)
    for c in worst: print(f'  class {c}: min certified margin -g - |Delta| = {worst[c][0]:.4f} at {worst[c][1]}',flush=True)
