# offpeak_common.py -- shared ball-arithmetic helpers for the off-peak lemma (B).
# Notation (GAP_CENTRE2): q = zeta^h e^{-sigma/(5n)}, sigma = L - i y, y = 5n(theta - 2 pi h/5),
# Phi_c(s) = log(1 - zeta^{hc} e^{-s}) (principal), Q = sum_{c=1..4} Phi_c,
# p(sigma) = (1/sigma) int_0^sigma Q = (Li2(e^{-5 sigma})/5 - Li2(e^{-sigma}) + 4 zeta(2)/5)/sigma   (h != 0),
# p0(sigma) = (1/sigma) int_0^sigma 4 log(1-e^{-s}) ds = 4 (Li2(e^{-sigma}) - zeta(2))/sigma      (h = 0).
# All closed forms use the principal Li2, valid because |e^{-s}| < 1 on the whole path (Re sigma = L > 0).
from flint import arb, acb, ctx
ctx.prec = 80
PI = arb.pi()
Z = [acb.exp_pi_i(acb(2*k)/5) for k in range(5)]
ZETA2 = PI**2/6

def Li2(z):
    return z.polylog(2)

def p_h(sig):            # h != 0
    return (Li2((-5*sig).exp())/5 - Li2((-sig).exp()) + 4*ZETA2/5)/sig

def p_0(sig):            # h == 0
    return 4*(Li2((-sig).exp()) - ZETA2)/sig

def p_real_point(L):     # p(L) for real L>0 (point value, decreasing in L)
    L = arb(L)
    return ((Li2(acb((-5*L).exp()))/5 - Li2(acb((-L).exp())) + 4*ZETA2/5)/L).real

def Phi(h, c, s):
    return (1 - Z[(h*c) % 5]*(-s).exp()).log()

def dPhi(h, c, s):
    w = Z[(h*c) % 5]*(-s).exp()
    return w/(1 - w)

def ddPhi(h, c, s):
    w = Z[(h*c) % 5]*(-s).exp()
    return -w/(1 - w)**2

def B2(x):
    x = arb(x)
    return x*x - x + arb(1)/6

def ball(lo, hi):          # arb enclosing [lo, hi] (lo, hi python floats/ints/Fractions)
    from fractions import Fraction
    lo = Fraction(lo); hi = Fraction(hi)
    m = (lo+hi)/2; r = (hi-lo)/2
    mid = arb(m.numerator)/m.denominator
    return mid + arb(0, (arb(r.numerator)/r.denominator).upper()) + arb(0, mid.rad())

def sigbox(Llo, Lhi, ylo, yhi):   # sigma = L - i y  over the box
    return acb(ball(Llo, Lhi), -ball(ylo, yhi))

def E_req(n):            # required deficit 1.5 log n + 6.3
    return 1.5*arb(n).log() + arb('6.3')   # AUDIT2 N3: exact decimal, not the float

def absl(z):              # rigorous lower bound (arb) for |z| over an acb rectangle
    def ml(x):
        m = abs(arb(x.mid())) - arb(x.rad())
        return m if m > 0 else arb(0)
    a = ml(z.real); b = ml(z.imag)
    return arb((a*a + b*b).lower()).sqrt().lower() if (a*a+b*b) > 0 else arb(0)

def ddPhi_abs_ub(h, c, s):  # upper bound of |Phi_c''(s)| = |w|/|1-w|^2 over the ball s
    w = Z[(h*c) % 5]*(-s).exp()
    lo = absl(1 - w)
    if not (lo > 0):
        return arb('inf')
    return (abs(w).upper()/(lo*lo)).upper()

def dPhi_abs_ub(h, c, s):
    w = Z[(h*c) % 5]*(-s).exp()
    lo = absl(1 - w)
    if not (lo > 0):
        return arb('inf')
    return (abs(w).upper()/lo).upper()
