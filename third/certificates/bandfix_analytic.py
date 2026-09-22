# Analytic band lemma for large n (classes 3,4 mod 5, 2N <= m < 4N), and lemma L1 (g(t) < 0 for large t).
# All bounds in arb ball arithmetic; only UPPER ends of upper bounds and LOWER ends of lower bounds are used.
# Inputs taken as given: the cited expansion |s(i) - a(i)P5(i)| <= Rgrp(i) with E_con = 27.4636 (i >= 1).
# Standard Bessel facts used:  (B1) I_nu(y) is decreasing in nu >= 0 for fixed y > 0;
#   (B2) I_{1/2}(y) = sqrt(2/(pi y)) sinh y,  I_{5/2}(y) = sqrt(2/(pi y)) [(1+3/y^2) sinh y - (3/y) cosh y].
# Hence  I_1(y) <= e^y/sqrt(2 pi y)  and  I_2(y) >= I_{5/2}(y) >= e^y/sqrt(2 pi y) * kap(y),
#   kap(y) = 1 - 3/y - 2 e^{-2y}   (y >= 3), increasing in y.
import sys, pickle, math
from flint import arb, ctx
ctx.prec = 200
PI = arb.pi(); a = 2 * PI**2 / 75; ECON = arb('27.4636'); sa = a.sqrt()
LOG = lambda x: arb(x).log()

def kap(Y): return 1 - 3 / Y - 2 * (-2 * Y).exp()

# ---- check (B2)-based lower bound numerically against arb's own I_2 at a few points (sanity only) ----
for yy in (5, 20, 100, 400):
    y = arb(yy); env = y.exp() / (2 * PI * y).sqrt() * kap(y)
    assert (y.bessel_i(2) - env).lower() > 0 and (y.bessel_i(2) - y.bessel_i(arb(5)/2)).lower() > 0

# ---- alpha: |s(i)| <= alpha * exp(2 sqrt(a i)) for all i >= 0 ----
s = pickle.load(open('g5_coeffs.pkl', 'rb')); M = len(s) - 1
alpha1 = arb(0)
for i in range(M + 1):
    if s[i]:
        v = abs(s[i]) * (-2 * (a * i).sqrt()).exp()
        if v.upper() > alpha1.upper(): alpha1 = arb(v.upper())
# i > M (cited bound): |s(i)| <= 3.62 P5 + Phi(K) P10 + E_con with P5 <= (2a/y) e^y/sqrt(2 pi y),
#   P10 <= (a/y) e^{y/2}/sqrt(pi y), Phi(K) <= K^2/10 + K/2, K <= c y + 3 (c = sqrt(pi/a)), y = 2 sqrt(a(i-1/6)) <= 2 sqrt(a i).
# Divided by e^y, each term is decreasing in y for y >= y_M (y^{-3/2}; poly(y) y^{-3/2} e^{-y/2} with deg 2 -> decreasing for y > 1;
# e^{-y}), so evaluate at y_M.
yM = 2 * (a * (arb(M + 1) - arb(1) / 6)).sqrt(); c = (PI / a).sqrt(); KM = c * yM + 3
alpha2 = arb('3.62') * (2 * a / yM) / (2 * PI * yM).sqrt() + (KM**2 / 10 + KM / 2) * (a / yM) * (-yM / 2).exp() / (PI * yM).sqrt() + ECON * (-yM).exp()
alpha = arb(max(alpha1.upper(), alpha2.upper()))
print('alpha1 (exact, i<=%d) <= %s ; alpha2 (cited, i>%d) <= %s ; alpha = %s' % (M, alpha1.str(5), M, alpha2.str(5), alpha.str(5)))

# ---- Lemma L1:  3 Rgrp(X+2) <= P5'(X)  for all X >= X0 ----
# P5'(X) = 4a^2 I2(Y)/Y^2 >= 4a^2 kap(Y) e^Y / (Y^2 sqrt(2 pi Y)),  Y = 2 sqrt(a(X-1/6)).
# Rgrp(X+2) <= Phi(K) sqrt(a/(4X)) e^{y10}/sqrt(2 pi y10) + E_con,  y10 = sqrt(a(X+11/6)) <= Y/2 + sqrt(a/(X-1/6)),
#   y10 >= Y/2, K <= sqrt(4 pi (X+2)) + 2.
# ratio <= RA(X) + RB(X) with
#   RA = 3 (K^2/10+K/2) sqrt(a/(4X)) e^{sqrt(a/(X-1/6))} Y^2 sqrt(2 pi Y) e^{-Y/2} / (sqrt(pi Y) 4 a^2 kap(Y))
#   RB = 3 E_con Y^2 sqrt(2 pi Y) e^{-Y} / (4 a^2 kap(Y)).
# Monotonicity for X >= X0: K^2/10+K/2 <= C1 X (C1 := value/X at X0, since (K^2/10+K/2)/X is decreasing: K^2/X
# = (sqrt(4pi(X+2))+2)^2/X decreasing), sqrt(a/4X) ~ X^{-1/2}, Y^2 ~ X, 1/kap decreasing, e^{sqrt(a/(X-1/6))} decreasing:
# RA <= const * X^{3/2} e^{-Y/2}; d/dX log = 3/(2X) - (1/2) sqrt(a/(X-1/6)) < 0 iff X > 9/a-ish (X >= 40 suffices).
# RB ~ X^{5/4} e^{-Y}: likewise decreasing.  So ratio(X) <= RA(X0)+RB(X0) for X >= X0 after bounding the non-monotone
# pieces by their X0 values -- we evaluate exactly the monotone dominating form below.
def L1_bound(X0):
    X = arb(X0); Y = 2 * (a * (X - arb(1) / 6)).sqrt(); K = (4 * PI * (X + 2)).sqrt() + 2
    C1 = (K**2 / 10 + K / 2) / X
    RA = 3 * C1 * X * (a / (4 * X)).sqrt() * (a / (X - arb(1) / 6)).sqrt().exp() * Y**2 * (2 * PI * Y).sqrt() * (-Y / 2).exp() / ((PI * Y).sqrt() * 4 * a**2 * kap(Y))
    RB = 3 * ECON * Y**2 * (2 * PI * Y).sqrt() * (-Y).exp() / (4 * a**2 * kap(Y))
    return RA + RB
for X0 in (1000, 5000, 14997):
    print('L1: sup_{X>=%d} 3 Rgrp(X+2)/P5\'(X) <= %s' % (X0, L1_bound(X0).str(5)))
assert L1_bound(1000).upper() < 1

# ---- Band majorant for N >= N0 ----
# ratio (|b2|+|b3|)/|g(T)| <= Psi(N) = beta * alpha * (12 a)^{5/4} sqrt(2 pi) e^{delta0} / (16 a^2 kap(Ymin)) * N^{9/4} (log N)^2 e^{-cc sqrt N}
cc = 2 * sa * (arb(3).sqrt() - arb(2).sqrt())
def Psi(N0, N=None):
    N0 = arb(N0); N = N0 if N is None else arb(N)
    L0 = N0.log()
    eps1 = (1 + 2 / (3 * (2 / a).sqrt() * N0.sqrt() * L0))**2 - 1
    Wfac = (N0 + arb(3) / 2 + (N0 + 3)**2 / 12) / (N0**2 / 12)      # Wmax <= Wfac N^2/12, Wfac decreasing
    eps3 = 2 * Wfac / 12 / (N0 * L0**2)
    beta = (9 / a) * (1 + eps1) + eps3
    delta0 = 4 * sa / (3 * N0 - 4).sqrt()
    Ymin = 2 * (a * (N0 - 4)).sqrt()
    C = beta * alpha * (12 * a)**(arb(5) / 4) * (2 * PI).sqrt() * delta0.exp() / (16 * a**2 * kap(Ymin))
    return C * N**(arb(9) / 4) * N.log()**2 * (-cc * N.sqrt()).exp()
if __name__ == '__main__':
    for N0 in (10001, 12001, 14001, 15001, 16001, 20001, 50001):
        print('Psi(N0=%d) <= %s' % (N0, Psi(N0).str(5)))
    # D < N requirement: 3 sqrt(2N/a) log N < N  and monotone decrease: cc sqrt N / 2 > 9/4 + 2/log N
    N0 = int(sys.argv[1]) if len(sys.argv) > 1 else 15001
    ok = (3 * (2 * arb(N0) / a).sqrt() * LOG(N0) < N0) and ((cc * arb(N0).sqrt() / 2 - arb(9) / 4 - 2 / LOG(N0)).lower() > 0)
    P = Psi(N0)
    print('N0=%d: D<N and Psi decreasing: %s;  Psi(N0) <= %s' % (N0, ok, P.str(5)))
    assert ok and P.upper() < 1
    print('BAND (analytic): (|b2|+|b3|)/|g(T)| < %.4g for all N >= %d, 2N <= m < 4N, classes 3,4' % (float(P.upper()), N0))
