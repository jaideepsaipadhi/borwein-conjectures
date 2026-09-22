# Scaled closed forms of the L** integrand, analytic in eps = 1/N (so one ball eps in [0, 1/N0] covers every N >= N0).
#   zeta = N w.   S(w) = N * St(zeta, eps),  QD(w) = Qt(zeta, eps),  Delta = QD e^{-zeta},
#   St = sum_r s_r e^{-r zeta} (1 + e^{-r zeta eps}) g(3 r zeta eps) / (3 r zeta),   g(x) = x/(1-e^{-x})
#   Qt = sum_{3 not| r} c_r chi(r) e^{-(r-1) zeta} h(r zeta eps),                   h(x) = (1-e^{-x})/(1-e^{-3x}) = g(3x)/(3 g(x))
#   Rh(w) = N * G(zeta) + (m - N - 1/4) zeta / N,   G = a/zeta + St + eps (log Qt + log sinc(Delta)).
from flint import arb, acb, arb_series, acb_series, ctx
import math
ctx.prec=100
PI=arb.pi(); a=PI**2/9; S3=arb(3).sqrt()
# Bernoulli^+ coefficients b_k = B_k^+ / k!  (x/(1-e^{-x}) = sum b_k x^k), |b_k| <= 4/(2pi)^k
from fractions import Fraction
def _bern(K):
    B=[Fraction(1)]
    for mm in range(1,K+1):
        B.append(-sum(Fraction(math.comb(mm+1,k))*B[k] for k in range(mm))/(mm+1))
    B[1]=Fraction(1,2)
    return [arb(Bk.numerator)/arb(Bk.denominator)/math.factorial(k) for k,Bk in enumerate(B)]
KB=24; BC=_bern(KB)
def g_acb(x):
    s=acb(0); p=acb(1)
    for k in range(KB+1): s+=BC[k]*p; p=p*x
    X=abs(x).upper()/(2*math.pi)
    assert X<0.9, 'g argument too large'
    rem=4*X**(KB+1)/(1-X)
    return s+acb(arb(0,rem),arb(0,rem))
def g_series(x0, x1, order):
    # g(x0 + x1 u) as series in u; x0, x1 real balls; truncated Bernoulli polynomial + per-coefficient remainder
    u=arb_series([x0,x1],prec=order); s=arb_series([0],prec=order); p=arb_series([1],prec=order)
    for k in range(KB+1): s+=BC[k]*p; p=p*u
    X0=float(abs(x0).upper()); X1=float(abs(x1).upper())
    assert X0/(2*math.pi)<0.9
    cs=list(s.coeffs())+[arb(0)]*order
    out=[]
    for j in range(order):
        tot=0.0
        for k in range(max(KB+1,j),KB+400):
            lg=math.log(4)+math.lgamma(k+1)-math.lgamma(j+1)-math.lgamma(k-j+1)+(k-j)*(math.log(X0) if X0>0 else (0 if k==j else -1e9))-k*math.log(2*math.pi)
            if lg>-1e8: tot+=math.exp(lg)
        tot+=1e-300
        out.append(cs[j]+arb(0,tot*X1**j*1.01))
    return arb_series(out,prec=order)
def chi(r): return 1 if r%3==1 else -1
def Rcut(v): return int(160/max(v,1.0))+6

def G_acb(z, eps, vlo, extras=False):
    """G(z) for complex z (acb), eps ball in [0, 1/N0]; Re z >= vlo.  Terms whose g/h argument is not small are replaced
    by rigorous magnitude balls:  |g(3 r z eps)/(3 r z)| = |eps/(1-e^{-3 r z eps})| <= 1/(3 r Re z) + eps,
    |h(r z eps)| <= 2 (1/(3 r eps Re z) + 1) -- only used for r with |r z eps| > 0.9 (then 1/eps <= r|z|/0.9)."""
    R=Rcut(vlo); St=acb(0); Qt=acb(0)
    eh=float(eps.upper()); az=float(abs(z).upper()); rez=float(z.real.lower())
    for r in range(1,R+1):
        sr=arb(3)/r if r%3==0 else -arb(3)/(2*r)
        e1=(-r*z).exp()
        if 3*r*az*eh/(2*math.pi) < 0.8:
            St+=sr*e1*(1+(-r*z*eps).exp())*g_acb(3*r*z*eps)/(3*r*z)
        else:
            m_=float(abs(sr).upper())*float(abs(e1).upper())*2*(1/(3*r*rez)+eh)
            St+=acb(arb(0,m_),arb(0,m_))
        if r%3:
            ec=(-(r-1)*z).exp()*(3*S3/(2*r))*chi(r)
            if 3*r*az*eh/(2*math.pi) < 0.8:
                Qt+=ec*g_acb(3*r*z*eps)/(3*g_acb(r*z*eps))
            else:
                m_=float(abs(ec).upper())*2*(r*az/(0.9*3*r*rez)+1)
                Qt+=acb(arb(0,m_),arb(0,m_))
    geo=math.exp(-(R+1)*vlo)/(1-math.exp(-vlo))
    tS=3*geo*2*(1/(3*(R+1)*rez)+eh); tQ=6*math.exp(vlo-(R+1)*vlo)/(1-math.exp(-vlo))*2*(az/(2.7*rez)+1)
    St+=acb(arb(0,tS),arb(0,tS)); Qt+=acb(arb(0,tQ),arb(0,tQ))
    D=Qt*(-z).exp()
    sinc=sinc_acb(D)
    G=a/z+St+eps*(Qt.log()+sinc.log())
    if extras: return G, St, Qt, D, sinc
    return G
def sinc_acb(D):
    u=D*D; g=acb(1); term=acb(1)
    for k in range(1,7):
        term=term*u; g+=term*((-1)**k/arb(math.factorial(2*k+1)))
    rem=abs(D).upper()**14/math.factorial(15)*math.cosh(abs(D).upper())
    return g+acb(arb(0,rem),arb(0,rem))

def G_series(v, eps, order):
    """Taylor series of G at real v (arb ball), eps ball."""
    t=arb_series([v,1],prec=order); one=arb_series([1],prec=order)
    R=Rcut(float(v.lower())); St=0*one; Qt=0*one
    for r in range(1,R+1):
        sr=arb(3)/r if r%3==0 else -arb(3)/(2*r)
        gx=g_series(3*r*v*eps,3*r*eps,order)
        St+=sr*(-r*t).exp()*(1+(-r*eps*t).exp())*gx/(3*r*t)
        if r%3:
            h=g_series(3*r*v*eps,3*r*eps,order)/(3*g_series(r*v*eps,r*eps,order))
            Qt+=(3*S3/(2*r))*chi(r)*(-(r-1)*t).exp()*h
    vlo=float(v.lower()); geo=math.exp(-(R+1)*vlo/2)/(1-math.exp(-vlo/2))
    cs=[]
    Sc=list(St.coeffs())+[arb(0)]*order; Qc=list(Qt.coeffs())+[arb(0)]*order
    for j in range(order):
        grow=(2.0*(R+1)+2)**j
        cs.append((Sc[j]+arb(0,6*geo*grow/(3*vlo)), Qc[j]+arb(0,6*math.exp(vlo-(R+1)*vlo/2)/(1-math.exp(-vlo/2))*grow)))
    St=arb_series([c[0] for c in cs],prec=order); Qt=arb_series([c[1] for c in cs],prec=order)
    D=Qt*(-t).exp()
    u=D*D; sn=arb_series([1],prec=order); term=arb_series([1],prec=order)
    for k in range(1,7):
        term=term*u; sn+=term*((-1)**k/arb(math.factorial(2*k+1)))
    Dabs=float(abs(D.coeffs()[0]).upper())
    rem=Dabs**14/math.factorial(15)
    snc=[c+arb(0,rem*(14*4)**j) for j,c in enumerate(list(sn.coeffs())+[arb(0)]*(order-len(sn.coeffs())))][:order]
    sn=arb_series(snc,prec=order)
    return a/t+St+eps*(Qt.log()+sn.log())

def parts_acb(z, eps, vlo):
    """(G0, b) with G0 = a/z + St, b = log Qt + log sinc(Delta), so G = G0 + eps b."""
    G, St, Qt, D, snc = G_acb(z, eps, vlo, extras=True)
    return a/z+St, Qt.log()+snc.log(), St

def parts_series(v, eps, order):
    t=arb_series([v,1],prec=order); one=arb_series([1],prec=order)
    R=Rcut(float(v.lower())); St=0*one; Qt=0*one
    for r in range(1,R+1):
        sr=arb(3)/r if r%3==0 else -arb(3)/(2*r)
        gx=g_series(3*r*v*eps,3*r*eps,order)
        St+=sr*(-r*t).exp()*(1+(-r*eps*t).exp())*gx/(3*r*t)
        if r%3:
            h=g_series(3*r*v*eps,3*r*eps,order)/(3*g_series(r*v*eps,r*eps,order))
            Qt+=(3*S3/(2*r))*chi(r)*(-(r-1)*t).exp()*h
    vlo=float(v.lower()); geo=math.exp(-(R+1)*vlo/2)/(1-math.exp(-vlo/2))
    Sc=list(St.coeffs())+[arb(0)]*order; Qc=list(Qt.coeffs())+[arb(0)]*order
    St=arb_series([Sc[j]+arb(0,6*geo*(2.0*(R+1)+2)**j/(3*vlo)) for j in range(order)],prec=order)
    Qt=arb_series([Qc[j]+arb(0,6*math.exp(vlo-(R+1)*vlo/2)/(1-math.exp(-vlo/2))*(2.0*(R+1)+2)**j) for j in range(order)],prec=order)
    D=Qt*(-t).exp()
    u=D*D; sn=arb_series([1],prec=order); term=arb_series([1],prec=order)
    for k in range(1,7):
        term=term*u; sn+=term*((-1)**k/arb(math.factorial(2*k+1)))
    Dabs=float(abs(D.coeffs()[0]).upper()); rem=Dabs**14/math.factorial(15)
    snc=list(sn.coeffs())+[arb(0)]*order
    sn=arb_series([snc[j]+arb(0,rem*(56.0)**j) for j in range(order)],prec=order)
    return a/t+St, Qt.log()+sn.log()

def A_diff(V, t, eps):
    """Re[G0(V+it) - G0(V)], G0 = a/z + St, with the V-dependency controlled:
       a-part = -a t^2/(V(V^2+t^2)) (increasing in V: upper at V_hi, lower at V_lo);
       St-part = Re sum_r s_r e^{-rV} [ e^{-irt} k_r(V+it) - k_r(V) ],  k_r(z) = (1+e^{-r z eps}) g(3 r z eps)/(3 r z)."""
    Vl=arb(V.lower()); Vh=arb(V.upper())
    ap=arb.union(-a*t*t/(Vl*(Vl*Vl+t*t)), -a*t*t/(Vh*(Vh*Vh+t*t)))
    vlo=float(V.lower()); R=Rcut(vlo); tot=acb(0)
    z=acb(V,t); eh=float(eps.upper()); az=float(abs(z).upper())
    for r in range(1,R+1):
        sr=arb(3)/r if r%3==0 else -arb(3)/(2*r)
        if 3*r*az*eh/(2*math.pi) < 0.8:
            kz=(1+(-r*z*eps).exp())*g_acb(3*r*z*eps)/(3*r*z)
            kv=(1+(-r*V*eps).exp())*g_acb(3*r*acb(V)*eps)/(3*r*V)
            br=(-acb(0,r)*t).exp()*kz-kv
            tot+=sr*(-r*V).exp()*br
        else:
            m_=float(abs(sr).upper())*math.exp(-r*vlo)*2*2*(1/(3*r*vlo)+eh)
            tot+=acb(arb(0,m_),arb(0,m_))
    geo=math.exp(-(R+1)*vlo)/(1-math.exp(-vlo))
    tS=2*3*geo*2*(1/(3*(R+1)*vlo)+eh)
    return ap+tot.real+arb(0,tS)
