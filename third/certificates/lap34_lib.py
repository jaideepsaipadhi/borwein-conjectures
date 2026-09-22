# lap34_lib.py -- library for the classes-3,4 Laplace certificate (GAP_LAP34.md).
# Built on cls012_lib.py (gamma, Q, alpha, Xi_V and its Taylor series).  New ingredient: the THIN FUNCTION
#     g(s) = e^{s} * sum_{h=1..4} A_h z5^{h cls} exp(Omega_h(s)),   Omega_h(s) = sum_c B1(c/5) log(1 - z5^{-hc} e^{-s})
# (s = 5 n u = V zeta).  For classes 3,4, sum_h A_h z5^{h cls} = a(cls) = 0, so the bracket is O(e^{-s}) and the
# factor e^{s} (the "-N w" part of log F_m, i.e. the m -> m - 5n saddle shift) is moved into the exponent.
exec(open('cls012_lib.py').read())
AH={1:acb((PI/5).cos(),(PI/5).sin()),4:acb((PI/5).cos(),-(PI/5).sin()),2:acb(1),3:acb(1)}
def z5(k):
    a=2*PI*(k%5)/5; return acb(a.cos(),a.sin())
B1={1:arb(-3)/10,2:arb(-1)/10,3:arb(1)/10,4:arb(3)/10}
ECON=arb('27.4636')
def lo(x): return arb(x.lower()) if isinstance(x,arb) else arb(x.real.lower())
def up(x): return arb(x.upper()) if isinstance(x,arb) else arb(abs(x).upper())
def amax(*xs): return max(xs, key=lambda z: float(z.lower()))
def gfun(s,cls):
    """g(s) for complex ball s (Re s > 0)."""
    s=acb(s); y=(-s).exp(); tot=acb(0)
    for h in (1,2,3,4):
        Om=acb(0)
        for c in (1,2,3,4): Om+=B1[c]*(1-z5(-h*c)*y).log()
        tot+=AH[h]*z5(h*cls)*Om.exp()
    return s.exp()*tot
def gder(s,cls):
    """g'(s) via a length-2 Taylor series (ball s)."""
    ctx.cap=2
    ss=acb_series([acb(s),acb(1)]); y=(-ss).exp(); tot=acb_series([acb(0)])
    for h in (1,2,3,4):
        Om=acb_series([acb(0)])
        for c in (1,2,3,4): Om=Om+((1-y*z5(-h*c)).log())*B1[c]
        tot=tot+Om.exp()*(AH[h]*z5(h*cls))
    r=ss.exp()*tot
    return r.coeffs()[1]
def finf(cls):
    """first-order coefficient: g(s) -> finf as Re s -> infinity."""
    return -sum((AH[h]*z5(h*cls)*sum((B1[c]*z5(-h*c) for c in (1,2,3,4)),acb(0)) for h in (1,2,3,4)),acb(0))
