# offpeak_em0.py -- piece (E0) of the off-peak lemma (B): theta near 0 (h = 0), |y| <= YE0, y = 5 n theta.
# Certified (ball arithmetic over full boxes) for all n >= N0, L in [0.1, 2.1]; L in [0,0.1] via the
# monotonicity lemma (cost 0.1 n, charged on the first L-box).  |S(conj q)| = |S(q)| gives y -> -y, so y >= 0.
#
# With Phi(s) = log(1 - e^{-s}), omega = sigma/n, all four classes c share Phi, and
#   log S_n = sum_c Phi(c omega/5) + sum_c sum_{t=1}^{n-1} Phi((t + c/5) omega).
# Finite EM (Lemma A1) on t = 1..n-1 (F(x) = Phi((x+1) omega), n-1 terms), sum_c B1(c/5) = 0, gives
#   log|S_n| <= n Re p0(sigma) + Psi(omega) + sum_c |B2(c/5)|/2 |omega| (|Phi'(sigma)| + |Phi'(omega)|)
#               + (|omega||sigma|/3) int_{1/n}^1 |Phi''(t sigma)| dt,
#   p0(sigma) = 4 (Li2(e^{-sigma}) - zeta(2))/sigma,
#   Psi(omega) = Re[ sum_c Phi(c omega/5) - (4/omega) int_0^omega Phi ]
#              = 4 + log(24/625) + sum_c Re l(c omega/5) - 4 Re[(1/omega) int_0^omega l],  l(s) = log((1-e^{-s})/s).
# Analytic facts (proved in GAP_OFFPEAK.md): for Re s >= 0, |s| <= 0.1:  |(1-e^{-s})/s - 1| <= |s| e^{|s|}/2,
#   hence |l(s)| <= -log(1 - |s|e^{|s|}/2) <= |s|, and |s Phi'(s)| = |e^{-s}| / |(1-e^{-s})/s| <= 1/(1-|s|e^{|s|}/2).
# So Psi <= 4 + log(24/625) + 6|omega|, |omega| <= |sigma|/N0.
# int_{1/n}^{t0}: |Phi''(s)| <= 1/(|s|^2 m^2), m = 1 - r e^r/2, r = t0 |sigma|  ->  contributes <= 1/(3 m^2).
# int_{t0}^1: t-pieces with ball enclosures, multiplied by |sigma|^2/(3n) <= |sigma|^2/(3 N0).
import sys, time
from offpeak_common import *

N0 = int(sys.argv[1]) if len(sys.argv) > 1 else 821
YE0 = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0
NT = int(sys.argv[3]) if len(sys.argv) > 3 else 32
LLO, LHI = 0.1, 2.1
n = arb(N0)
ER = E_req(N0)
T0 = arb(1)/64
LOG24_625 = arb(24).log() - arb(625).log()

def Phi0(s):   return (1 - (-s).exp()).log()
def dPhi0_ub(s):
    w = (-s).exp(); lo = absl(1 - w)
    return (abs(w).upper()/lo).upper() if lo > 0 else arb('inf')
def ddPhi0_ub(s):
    w = (-s).exp(); lo = absl(1 - w)
    return (abs(w).upper()/(lo*lo)).upper() if lo > 0 else arb('inf')
def dPhi0(s):
    w = (-s).exp(); return w/(1 - w)

def box_bound(Llo, Lhi, ylo, yhi, nt=NT):
    sig = sigbox(Llo, Lhi, ylo, yhi)
    sc = acb(sig.real.mid(), sig.imag.mid()); dsig = sig - sc
    chi = 1 if Llo <= LLO else 0
    asig = abs(sig).upper()
    om = arb(asig)/n                         # upper bound of |omega| for all n >= N0
    assert om < arb(0.1)
    # p0 mean-value form: p0'(s) = int_0^1 t Q0'(t s) dt, Q0 = 4 Phi.
    pd = acb(0); I = arb(0)
    r = T0*asig; m = 1 - r*r.exp()/2
    assert m > 0
    # piece 0: t Phi'(t sig) = (1/sig) x/(e^x - 1), x = t sig, |x| <= r0 = |sig|/nt, and
    # |x/(e^x-1)| <= 1/(1 - r0 e^{r0}/2)  (since |(e^x-1)/x - 1| <= r0 e^{r0}/2).
    r0 = arb(asig)/nt; b0 = 1/(1 - r0*r0.exp()/2)
    assert b0 > 0
    pd += 4/(sig*nt)*acb(arb(0, b0.upper()), arb(0, b0.upper()))
    for i in range(1, nt):
        t = ball(i/nt, (i+1)/nt)
        ts = acb(t)*sig
        pd += (arb(i+1)**2 - arb(i)**2)/(2*nt*nt)*4*dPhi0(ts)
    # int_{t0}^1 |Phi''(t sig)| dt over pieces of [t0,1]
    for i in range(nt):
        a = float(T0) + (1-float(T0))*i/nt; b = float(T0) + (1-float(T0))*(i+1)/nt
        t = ball(a, b)
        I += ddPhi0_ub(acb(t)*sig)*(arb(b)-arb(a))
    P = p_0(sc) + pd*dsig
    g = p_real_point(Lhi) - P.real - arb(0.1)*chi
    Psi = 4 + LOG24_625 + 6*om
    EMb = arb(24)/300*(1/(1 - om*om.exp()/2) + om*dPhi0_ub(sig))
    Rem = 1/(3*m*m) + arb(asig)**2/(3*n)*I
    V = -n*g + Psi + EMb + Rem + ER
    ok = V.is_finite() and g.is_finite() and (V.upper() < 0) and (g.lower() > 1.5/N0)
    return ok, V, g, Rem

def certify(Llo, Lhi, ylo, yhi, depth, stats):
    nt = NT
    ok, V, g, R = box_bound(Llo, Lhi, ylo, yhi, nt)
    while (not ok) and nt < 1024 and not (V.is_finite() and g.is_finite()):
        nt *= 4
        ok, V, g, R = box_bound(Llo, Lhi, ylo, yhi, nt)
    if ok:
        stats['boxes'] += 1; stats['worstV'] = max(stats['worstV'], float(V.upper()))
        stats['maxR'] = max(stats['maxR'], float(R.upper()))
        return True
    if depth > 14:
        print('FAIL L=[%.5f,%.5f] y=[%.5f,%.5f] V=%s g=%s R=%s' % (Llo, Lhi, ylo, yhi, V, g, R), flush=True)
        return False
    if (Lhi-Llo)/0.02 >= (yhi-ylo)/0.1:
        mm = (Llo+Lhi)/2
        return certify(Llo, mm, ylo, yhi, depth+1, stats) and certify(mm, Lhi, ylo, yhi, depth+1, stats)
    mm = (ylo+yhi)/2
    return certify(Llo, Lhi, ylo, mm, depth+1, stats) and certify(Llo, Lhi, mm, yhi, depth+1, stats)

if __name__ == '__main__':
    t0 = time.time(); stats = {'boxes': 0, 'worstV': -1e300, 'maxR': 0.0}; ok = True
    Ls = [LLO + 0.05*i for i in range(41)]; Ls[0] = 0.0999   # AUDIT: cover the exact 0.1
    ys = [0.25*i for i in range(int(round(YE0/0.25))+1)]
    for i in range(len(Ls)-1):
        for j in range(len(ys)-1):
            ok = certify(Ls[i], Ls[i+1], ys[j], ys[j+1], 0, stats) and ok
    print('N0=%d YE0=%.1f: piece (E0) %s  boxes=%d worst V(N0)=%.2f max Rem=%.2f (%.0fs)' % (
        N0, YE0, 'CERTIFIED for all n >= N0' if ok else 'FAILED', stats['boxes'], stats['worstV'], stats['maxR'], time.time()-t0))
