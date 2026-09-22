# Laplace-regime certificate at fixed n, built directly on the EXACT closed forms of the L** integrand
# (no continuum phase, no K3/decay/disc constants).  Replaces B_gap4..7, whose phase used mu instead of mu'.
#
#   c(m) = -2 Im J + R(m),   Im J0 = (1/2pi) int_{-pi}^{pi} Re e^{Rh(W+iy)} dy   on ANY line Re w = W > 0,
#   Rh(w) = a/w + (m - 1/4) w + S(w) + log sin Delta(w),     Delta = QD(w) e^{-N w},
#   S  = sum_r s_r (e^{-Nrw} + e^{-(N+1)rw})/(1-e^{-3rw}),  s_r = -3/(2r) (3 not| r), +3/r (3|r)
#   QD = sum_{3 not| r} (3 sqrt3/(2r)) chi(r) e^{-N(r-1)w} (1-e^{-rw})/(1-e^{-3rw}).
#   J = J0 + J_eps with |integrand of J_eps| <= 3 |F(e^-w)| e^{(m-1/4)Re w}   (P(q) = e^{w/4}(e^{a/w} - 1 + eps), |eps| <= 2).
# W is the real saddle of Rh for each m in the m-ball (interval Newton).  With t = y/sigma, sigma = Rh''(W)^{-1/2}:
#   window |t| <= c:  Rh = Rh0 - t^2/2 - i G3 t^3/6 + phi4 t^4/24 + E5,  |E5| <= G5 |t|^5/120 (Cauchy on the non-1/w part)
#   rest   |y| in [T, Y]: exact cell enclosures of D(y) = Re Rh(W+iy) - Rh(W)  (m drops out exactly)
#   crude  |y| in [Y, pi]: |F| <= E_n(e^-W), |P| <= e^{aW/(W^2+Y^2)} + 3
#   R(m):  rigorous group bound (sz_joint expansion, real axis).
# Certified inequality per m-ball:  (rest + eps + crude + Rrel) < sigma * int_{-c}^{c} L(t) dt   (L = pointwise lower bound of Re part).
from flint import arb, acb, arb_series, acb_series, ctx
import sys, math
ctx.prec=100
PI=arb.pi(); a=PI**2/9; S3=arb(3).sqrt()
n=int(sys.argv[1]); N=3*n+1

def chi(r): return 1 if r%3==1 else -1
def Rcut(v): return int(160/max(v,1.0))+6

def pieces(t, R, one):
    """S and QD at t (arb_series / acb / acb_series), truncated at R terms (tails added separately)."""
    S=0*one; QD=0*one
    for r in range(1,R+1):
        den=1-(-3*r*t).exp()
        S+= (arb(3)/r if r%3==0 else -arb(3)/(2*r)) * ((-N*r*t).exp()+(-(N+1)*r*t).exp())/den
        if r%3:
            QD+= (3*S3/(2*r))*chi(r)*(-N*(r-1)*t).exp()*(1-(-r*t).exp())/den
    return S, QD

def tails(Wlo, R, j):
    """bound on |d^j/dw^j| of the dropped terms r>R of S and of QD (generous), for Re w >= Wlo."""
    v=N*Wlo; r1=R+1
    geo=(-(r1)*v/2).exp()/(1-(-v/2).exp())
    inv=(1+3*r1*Wlo)/(3*r1*Wlo)
    grow=(arb(6*N*r1))**j
    tS=6*geo*inv**(j+1)*grow
    tQ=6*(v/2).exp()*geo*inv**(j+1)*grow
    return tS, tQ

def sinc_series(Dl, order, Dabs, scale):
    # sin(x)/x = sum_{k<=6} (-1)^k x^{2k}/(2k+1)!  + remainder;  remainder bound radius |x|^14/15! * (2*scale)^j on coefficient j
    u=Dl*Dl; g=arb_series([1],prec=order); term=arb_series([1],prec=order)
    for k in range(1,7):
        term=term*u; g+= term*((-1)**k/arb(math.factorial(2*k+1)))
    cs=g.coeffs()+[arb(0)]*(order-len(g.coeffs()))
    rem=Dabs**14/arb(math.factorial(15))
    cs=[cs[j]+arb(0,(rem*(2*scale)**j).upper()) for j in range(order)]
    return arb_series(cs,prec=order)
def sinc_acb(D):
    u=D*D; g=acb(1); term=acb(1)
    for k in range(1,7):
        term=term*u; g+=term*((-1)**k/arb(math.factorial(2*k+1)))
    rem=abs(D).upper()**14/math.factorial(15)*math.cosh(abs(D).upper())
    return g+acb(arb(0,rem),arb(0,rem))
def Rh_series(W, m, order):
    t=arb_series([W,1],prec=order)
    R=Rcut(float(N*W.lower()))
    S,QD=pieces(t,R,arb_series([1],prec=order))
    tl=[tails(arb(W.lower()),R,j) for j in range(order)]
    Sc=[S.coeffs()[j]+arb(0,(tl[j][0]/math.factorial(j)).upper()) for j in range(order)]
    Qc=[QD.coeffs()[j]+arb(0,(tl[j][1]/math.factorial(j)).upper()) for j in range(order)]
    S=arb_series(Sc,prec=order); QD=arb_series(Qc,prec=order)
    Dl=QD*(-N*t).exp()
    Dabs=arb(abs(Dl.coeffs()[0]).upper())+arb('1e-30')
    return a/t+(m-N-arb(1)/4)*t+S+QD.log()+sinc_series(Dl,order,Dabs,arb(7*(N+10))).log()

def hval(z, Wlo):
    """h(z) = S(z) + log QD(z) + log(sin Delta/Delta)  (complex), with tails as ball radius."""
    R=Rcut(float(N*Wlo))
    S,QD=pieces(z,R,acb(1))
    tS,tQ=tails(arb(Wlo),R,0)
    S+=acb(arb(0,tS.upper()),arb(0,tS.upper())); QD+=acb(arb(0,tQ.upper()),arb(0,tQ.upper()))
    Dl=QD*(-N*z).exp()
    return S, QD, Dl

def saddle(mb, w0):
    w=arb(w0)
    for _ in range(60):
        s=Rh_series(w,arb(mb.mid()),3); c=s.coeffs()
        wn=float((w-c[1]/(2*c[2])).mid()); wn=min(max(wn,0.5*float(w.mid())),2*float(w.mid())); w=arb(wn)
    r=abs(w)*arb('1e-7')
    for _ in(range(40)):
        X=arb(w.mid(),r.upper())
        c0=Rh_series(w,mb,2).coeffs(); cX=Rh_series(X,mb,3).coeffs()
        Nx=w-c0[1]/(2*cX[2])
        if Nx.lower()>X.lower() and Nx.upper()<X.upper(): return Nx
        r*=3
        if r.upper()>0.2*float(w.mid()): return None
    return None

def logEn_real(w):          # log E_n(e^-w) <= N Psi(Nw)/(Nw) + 6|log(1-e^{-Nw})|, Psi(v) <= 2e^-v/(1-e^-v)
    v=N*w; e=(-v).exp()
    return N*2*e/(1-e)/v-6*(1-e).log()

def certify_ball(mb, w0, cs=(2.0,2.5,3.0,3.5,4.0,4.5), dbg=False):
    W=saddle(mb,w0)
    if W is None: return None,None,'saddle'
    Wlo=float(W.lower()); Wmid=arb(W.mid())
    ser=Rh_series(W,mb,6).coeffs()
    R2=2*ser[2]; R3=6*ser[3]; R4=24*ser[4]
    if not R2.lower()>0: return None,W,'R2'
    sig_lo=1/R2.upper()**arb(0.5) if False else 1/arb(R2.upper()).sqrt(); sig_hi=1/arb(R2.lower()).sqrt()
    G3=abs(R3)*sig_hi**3
    ph4_lo=arb((R4*(sig_lo**4 if R4.lower()>0 else sig_hi**4)).lower()); ph4_hi=arb((R4*(sig_hi**4 if R4.upper()>0 else sig_lo**4)).upper())
    # base values at W (real)
    S0,QD0,D0=hval(acb(W),Wlo)
    logF0=S0.real+QD0.real.log()+sinc_acb(D0).real.log()-N*W        # log(e^S sin Delta) at W
    best=None
    rho=arb(Wlo)*arb('0.8')
    Sc,Qc,Dc=hval(acb(Wmid),Wlo); hc=Sc+Qc.log()+sinc_acb(Dc).log()
    mx=arb(0); K=96
    for i in range(K):
        th=arb(2*math.pi*(i+0.5)/K, math.pi/K)
        z=acb(Wmid)+rho*acb(th.cos(),th.sin())
        Sz,Qz,Dz=hval(z,Wlo*0.2)
        if not Qz.real.lower()>0: return None,W,'Q on circle'
        hz=Sz+Qz.log()+sinc_acb(Dz).log()
        mx=max(mx,arb(abs(hz-hc).upper()),key=lambda q:q.upper())
    # crude-region start Y (independent of c)
    v=N*Wlo; width=min(0.2/(N*max(1.0,6.0/v)), float(sig_lo.lower())/8)
    Y=None; crude=arb(0)
    for kap in (1.5,2,3,5,10,20,50,200):
        Yc=arb(kap)*arb(Wlo)
        if Yc.upper()>math.pi: break
        cr=2*(PI-Yc)*((a*W/(W*W+Yc*Yc)).exp()+3)*(logEn_real(arb(Wlo))-a/W-logF0).exp()
        if cr.upper()<1e-18: Y=Yc; crude=cr; break
    if Y is None: Y=PI
    # exact cells on [T_min, Y], stored once
    base_a=a/W; cells_=[]; y=float((arb(min(cs))*sig_lo).lower()); Yf=float(Y.upper())
    while y<Yf:
        w_=min(width, Yf-y)
        z=acb(W,arb(y+w_/2,w_/2))
        Sz,Qz,Dz=hval(z,Wlo)
        ra=(a/z).real
        Dy=ra-base_a+Sz.real-S0.real+abs(Qz).log()-QD0.real.log()+abs(sinc_acb(Dz)).log()-sinc_acb(D0).real.log()
        cells_.append((y, 2*arb(Dy.upper()).exp()*w_, 2*3*arb((Dy-ra).upper()).exp()*w_))
        y+=w_
    m=mb
    logRh0=a/W+(m-arb(1)/4)*W+logF0
    def grp(ak,lead):
        bestv=None
        for tt in (0.15,0.25,0.35,0.5,0.7,1.0,1.4,2.0,3.0):
            wp=arb(tt)*W
            val=(lead+m*wp+ak/wp+logEn_real(arb(wp.lower()))-logRh0).upper()
            bestv=val if bestv is None or val<bestv else bestv
        return arb(bestv).exp()
    KF=(4*PI*m).sqrt()+2
    Rabs=grp(a/4,arb(6).log())+grp(a/9,arb(18).log())+grp(a/16,arb(12).log())+grp(a/25,(3*(KF*KF/6+KF)).log())+grp(arb(0),arb('50.1501').log())
    Rrel=Rabs*PI
    for c in cs:
        T_hi=arb(c)*sig_hi; T_lo=arb(c)*sig_lo
        d=T_hi+arb(W.rad())
        if not (d.upper()<rho.lower()*0.7):
            if dbg: print('   c',c,'skip window',float(d.upper()),float(rho.lower()))
            continue
        M5=120*a/arb(Wlo)**6+120*mx*rho/(rho-d)**6
        G5=M5*sig_hi**5
        cells=400; Lint=arb(0)
        for i in range(cells):
            t1=arb(c*i/cells); t2=arb(c*(i+1)/cells)
            th3=G3*t2**3/6+G5*t2**5/120
            cosl=1-th3*th3/2
            ex_lo=(ph4_lo*(t1**4 if ph4_lo.lower()>=0 else t2**4)/24-G5*t2**5/120)
            if cosl.lower()>0: Lc=(-(t2*t2)/2).exp()*cosl*ex_lo.exp()
            else: Lc=(-(t1*t1)/2).exp()*cosl*(ph4_hi*t2**4/24+G5*t2**5/120).exp()
            Lint+=arb(Lc.lower())*(t2-t1)
        Lint=2*Lint
        if not Lint.lower()>0:
            if dbg: print('   c',c,'skip Lint',Lint,'G3',G3,'G5',G5)
            continue
        main_lo=sig_lo*Lint
        tl=float(T_lo.lower())
        rest=sum((r_ for (y0,r_,e_) in cells_ if y0+width>=tl),arb(0))
        eps=sum((e_ for (y0,r_,e_) in cells_ if y0+width>=tl),arb(0))
        eps+=3*2*T_hi*(ph4_hi*arb(c)**4/24+G5*arb(c)**5/120).exp()*(-(a*W/(W*W+T_hi*T_hi))).exp()
        B=(rest+eps+crude+Rrel)/main_lo
        if dbg: print('   c',c,'rest',float(rest.upper()),'Rrel',float(Rrel.upper()),'main_lo',float(main_lo.lower()),'B',float(B.upper()))
        if best is None or B.upper()<best.upper(): best=B
    return best, W, 'ok'

if __name__=='__main__':
    mu_a,mu_b,K=float(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4])
    edges=[mu_a*(mu_b/mu_a)**(j/K) for j in range(K+1)]
    worst=0; fails=[]
    for j in range(K):
        m_lo=math.ceil(edges[j]*N*N); m_hi=math.floor(edges[j+1]*N*N)
        if m_hi<m_lo: continue
        mb=arb((m_lo+m_hi)/2,(m_hi-m_lo)/2)
        mup=max(float(mb.mid())-N,1.0)/N**2
        w0=math.sqrt(1.0966/mup)/N
        B,W,st=certify_ball(mb,w0)
        if B is None or not B.upper()<1: fails.append((round(edges[j],6),st if B is None else float(B.upper())))
        else: worst=max(worst,float(B.upper()))
        print('  mu [%.6f,%.6f] m in [%d,%d]: %s'%(edges[j],edges[j+1],m_lo,m_hi, st if B is None else '%.4g'%float(B.upper())),flush=True)
    print('EXACT-LAPLACE n0=%d mu in [%.5f,%.5f], %d balls: failures %d %s ; worst B = %.4f'%(n,mu_a,mu_b,K,len(fails),fails[:4],worst))
