# Library for the classes-0,1,2 Laplace-region certificate (GAP_CLASS012.md section B).
# Scaled variables: W = saddle abscissa, X = 1/W, V = 5nW, xi = e^{-V}, zeta = 1 + i eta (u = W zeta).
#   gamma(z) = [Li2(e^{-z}) - Li2(e^{-5z})/5]/5,  Q(z) = log(1-e^{-5z}) - log(1-e^{-z}) = -5 gamma'(z)
#   alpha(V) = a5 - gamma(V) - V Q(V)/5                     (saddle: (m - 1/6) W^2 = alpha(V))
#   Xi_V(eta) = (a5 - gamma(V zeta))/zeta + alpha(V) zeta   (exponent / X on the contour)
from flint import arb, acb, acb_series, arb_series, ctx
import math
ctx.prec=128
PI=arb.pi(); a5=2*PI**2/75; I=acb(0,1)
def Li2(z): return acb(z).polylog(2)
def gam(z):
    z=acb(z); return (Li2((-z).exp())-Li2((-5*z).exp())/5)/5
def Qf(z):
    z=acb(z); return (1-(-5*z).exp()).log()-(1-(-z).exp()).log()
def alpha(V):
    return (a5-gam(V)-V*Qf(V)/5).real
def Xi_series(V,etac,K=5):
    """Taylor coefficients (in eps) of Xi_V(etac+eps), etac may be a ball."""
    ctx.cap=K
    zc=acb(1,0)+I*acb(etac)
    zs=acb_series([zc,I])                      # zeta(eps)
    e1=(-(acb(V))*zs).exp(); e5=(-5*acb(V)*zs).exp()
    Qs=(1-e5).log()-(1-e1).log()               # Q(V zeta)
    gs=(Qs*(-I*acb(V)/5)).integral()           # gamma(V zeta) - gamma(V zc)
    gs=gs+acb_series([gam(acb(V)*zc)])
    al=acb(alpha(V))
    Xs=(acb_series([a5])-gs)*zs.inv()+al*zs
    return [Xs.coeffs()[k] if k<len(Xs.coeffs()) else acb(0) for k in range(K)]
def Xi_val(V,eta):
    z=acb(1,0)+I*acb(eta); return (a5-gam(acb(V)*z))/z+alpha(V)*z
