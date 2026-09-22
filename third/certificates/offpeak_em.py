# offpeak_em.py -- piece (E) of the off-peak lemma (B): near a 5th root, 1.1 <= |y| <= YE, h = 1..4.
# Certified (ball arithmetic over full boxes) for all n >= N0, L in [0.1, 2.1]; L in [0, 0.1] is reduced to
# L = 0.1 by the monotonicity lemma (log|S(L)| <= log|S(L')| + n (L'-L), L <= L'), which costs 0.1 n
# and is charged on the first L-box.
#
# Bound (finite Euler-Maclaurin, Lemma A1/A2 of GAP_CENTRE2 with the integral form of the remainder):
#   log|S_n(q)| <= n Re p(sigma) + Re C_h(sigma) + K_h(sigma)/n,
#   K_h = sum_c [ |B2(c/5)|/2 |sigma| (|Phi_c'(sigma)| + |Phi_c'(0)|) + |sigma|^2/12 int_0^1 |Phi_c''(t sigma)| dt ].
# Phi_c is analytic on Re s > 0 and at s = 0 (zeta^{hc} != 1), so the EM identity holds for every y (no disc
# restriction is needed; only the size of K matters).  int_0^1 |Phi''| <= sum over NT t-pieces of (1/NT) sup.
# Requirement:  V(n) := -n [p(L) - Re p(sigma) - 0.1 chi] + Re C + K/n + 1.5 log n + 6.3 <= 0  at n = N0, and
#   g - 0.1 chi >= 1.5/N0 (then V is non-increasing for n >= N0).  chi = 1 on boxes touching L = 0.1.
# Symmetry: |S(conj q)| = |S(q)| maps (h, y) -> (5-h, -y); so h in {1,2} with both signs of y covers h = 1..4.
import sys, time
from offpeak_common import *

N0 = int(sys.argv[1]) if len(sys.argv) > 1 else 821
YE = float(sys.argv[2]) if len(sys.argv) > 2 else 6.0
NT = int(sys.argv[3]) if len(sys.argv) > 3 else 32
LLO, LHI = 0.1, 2.1
n = arb(N0)
ER = E_req(N0)
phi0p = {(h, c): abs(dPhi(h, c, acb(0))).upper() for h in (1, 2) for c in range(1, 5)}
C0 = {(h, c): Phi(h, c, acb(0)) for h in (1, 2) for c in range(1, 5)}

def box_bound(h, Llo, Lhi, ylo, yhi, nt=NT):
    sig = sigbox(Llo, Lhi, ylo, yhi)
    sc = acb(sig.real.mid(), sig.imag.mid())          # exact centre
    dsig = sig - sc                                    # centred box
    chi = 1 if Llo <= LLO else 0
    asig = abs(sig)
    # mean-value forms: p(sig) in p(sc) + p'(box) dsig, p'(s) = int_0^1 t Q'(t s) dt;
    # C(sig) in C(sc) + C'(box) dsig.  Integrals over t-pieces: piece integral lies in (weight) x (convex ball).
    pd = acb(0); Cd = acb(0); C = acb(0); K = arb(0)
    for c in range(1, 5):
        C += (arb(c)/5 - arb(1)/2)*(Phi(h, c, sc) - C0[(h, c)])
        d1 = dPhi(h, c, sig)
        Cd += (arb(c)/5 - arb(1)/2)*d1
        I = arb(0)
        for i in range(nt):
            t = ball(i/nt, (i+1)/nt)
            ts = acb(t)*sig
            w = Z[(h*c) % 5]*(-ts).exp()
            q1 = w/(1 - w)
            pd += (arb(i+1)**2 - arb(i)**2)/(2*nt*nt)*q1
            I += ddPhi_abs_ub(h, c, ts)/nt
        K += abs(B2(arb(c)/5))/2*asig*(dPhi_abs_ub(h, c, sig) + phi0p[(h, c)]) + asig**2/12*I
    P = p_h(sc) + pd*dsig
    Cb = C + Cd*dsig
    g = p_real_point(Lhi) - P.real - arb(0.1)*chi
    V = -n*g + Cb.real + K/n + ER
    ok = V.is_finite() and g.is_finite() and (V.upper() < 0) and (g.lower() > 1.5/N0)
    return ok, V, g, K

def certify(h, Llo, Lhi, ylo, yhi, depth=0, stats=None):
    nt = NT
    ok, V, g, K = box_bound(h, Llo, Lhi, ylo, yhi, nt)
    while (not ok) and nt < 1024 and not (V.is_finite() and g.is_finite()):
        nt *= 4
        ok, V, g, K = box_bound(h, Llo, Lhi, ylo, yhi, nt)
    if ok:
        stats['boxes'] += 1
        stats['worstV'] = max(stats['worstV'], float(V.upper()))
        stats['maxK'] = max(stats['maxK'], float(K.upper()))
        return True
    if depth > 14:
        print('FAIL h=%d L=[%.5f,%.5f] y=[%.5f,%.5f] V=%s g=%s K=%s' % (h, Llo, Lhi, ylo, yhi, V, g, K), flush=True)
        return False
    # split the dimension with the larger scaled width
    if (Lhi-Llo)/0.02 >= (yhi-ylo)/0.05:
        m = (Llo+Lhi)/2
        return certify(h, Llo, m, ylo, yhi, depth+1, stats) and certify(h, m, Lhi, ylo, yhi, depth+1, stats)
    m = (ylo+yhi)/2
    return certify(h, Llo, Lhi, ylo, m, depth+1, stats) and certify(h, Llo, Lhi, m, yhi, depth+1, stats)

if __name__ == '__main__':
    t0 = time.time()
    allok = True
    for h in (1, 2):
        for sgn in (1, -1):
            stats = {'boxes': 0, 'worstV': -1e300, 'maxK': 0.0}
            ok = True
            Ls = [LLO + 0.02*i for i in range(101)]; Ls[0] = 0.0999   # AUDIT: cover the exact 0.1
            ys = [1.1 + 0.1*i for i in range(int(round((YE-1.1)/0.1))+1)] + [YE]
            ys = sorted(set([min(v, YE) for v in ys]))
            for i in range(len(Ls)-1):
                for j in range(len(ys)-1):
                    a, b = ys[j], ys[j+1]
                    lo, hi = (a, b) if sgn > 0 else (-b, -a)
                    ok = certify(h, Ls[i], Ls[i+1], lo, hi, 0, stats) and ok
            print('h=%d sign=%+d: %s  boxes=%d worst V(N0)=%.2f max K=%.1f  (%.0fs)' % (
                h, sgn, 'CERTIFIED' if ok else 'FAILED', stats['boxes'], stats['worstV'], stats['maxK'], time.time()-t0), flush=True)
            allok = allok and ok
    print('N0=%d YE=%.2f: piece (E) %s' % (N0, YE, 'CERTIFIED for all n >= N0' if allok else 'FAILED'))
