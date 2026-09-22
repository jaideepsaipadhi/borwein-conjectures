# offpeak_tb.py -- piece (N) of the off-peak lemma (B): near a 5th root (h = 0..4), YE_h <= |y| <= UMAX n.
# n-uniform analytic bound on the triangle-inequality majorant of the j-series, certified by ball arithmetic.
#
#   log|S_n(q)| = -Re sum_j G_j/j <= TB := sum_j |G_j|/j,   G_j = sum_{k<=5n, 5 nmid k} q^{kj}
#   G_j = (1 - e^{-j sigma}) A_j^c / (1 - e^{-j omega}),  A_j^c = sum_{c=1}^4 zeta^{hcj} e^{-c j omega/5},  omega = sigma/n.
#   Abar_j = |sum_c zeta^{hcj}| (h != 0: 1 if 5 nmid j, 4 if 5 | j;  h = 0: 4);  |A_j| <= min(4, Abar_j + 2 j|omega|),
#   |A_j| <= alpha(jL/n),  alpha(a) = sum_c e^{-ca/5}.
# Part I (j <= J1 = floor(pi n/y); z = j omega has |Im z| <= pi, 0 <= Re z <= 20):
#   1/|1-e^{-z}| <= 1/|z| + c1,  c1 = 1 + 1/pi  (max modulus of 1/(1-e^{-z}) - 1/z on [0,20]x[-pi,pi]; see GAP_OFFPEAK.md)
#   => Part I <= n M(sigma) + R1,  M(sigma) = sum_{j>=1} |1-e^{-j sigma}| Abar_j/(j^2 |sigma|),
#      R1 <= sum_{j<=J1} (1+e^{-jL}) [2 + c1 min(4, Abar_j + 2 j |omega|)]/j.
# Part II (j > J1; period k >= 1: jy/n in ((2k-1)pi, (2k+1)pi]):
#   t_k = (y/((2k-1) pi n)) beta_k alpha(a_k^-) [2 G_k + 2 min(M G_k, D_k (1 + log+(2 M G_k/D_k)))],
#   a_k^{+-} = (2k+-1) pi L/y, beta_k = 1 + e^{-(2k-1) pi n L/y}, G_k = 1/(1-e^{-a_k^-}), M = pi n/y + 1,
#   D_k = e^{a_k^+/2} pi n/(2y).
# Requirement: TB <= n p(L) - (1.5 log n + 6.3) - 0.1 n chi  (chi = 1 on the L-box touching 0.1; L < 0.1 is
#   reduced to L = 0.1 by the monotonicity lemma).
# (N-a) YE_h <= y <= Y2: (L,y)-boxes, value at n = N0 plus monotonicity in n (II non-increasing in n, R1 = A + B log n).
# (N-b) Y2 <= y <= UMAX n: (L,u)-boxes, u = y/n, explicit a + b log n + c log^2 n bounds, value at N0 + monotonicity.
import sys, time
from offpeak_common import *

N0 = int(sys.argv[1]) if len(sys.argv) > 1 else 821
UMAX = arb(sys.argv[2]) if len(sys.argv) > 2 else arb('0.25')
Y2 = 200
YE = {1: 6.0, 0: 20.0}
ASW = arb(2)
C1 = 1 + 1/PI
n0 = arb(N0)
ER0 = E_req(N0)
JM = 400

def alpha(a): return sum(((-arb(c)*a/5).exp() for c in range(1, 5)), arb(0))
def Abar(h, j): return 4 if (h == 0 or j % 5 == 0) else 1
def Li2r(x): return Li2(acb(x)).real
def mh(h, L):   # sum_j (1+e^{-jL}) Abar_j/j^2
    if h == 0: return 4*ZETA2 + 4*Li2r((-L).exp())
    return arb(28)/25*ZETA2 + Li2r((-L).exp()) + arb(3)/25*Li2r((-5*L).exp())
def Hup(J): return 1 + arb(J).log() if J >= 1 else arb(0)   # H_J <= 1 + log J

def M_sigma(h, sig, Llo):
    s = arb(0); asig = abs(sig); cap_base = (-arb(Llo)).exp()
    for j in range(1, JM+1):
        v = abs(1 - (-j*sig).exp()).upper()
        cap = (1 + cap_base**j).upper()
        s += (v if v < cap else cap)*Abar(h, j)/arb(j)**2
    s += 8/arb(JM)      # tail: sum_{j>JM} (1+e^{-jL}) 4/j^2 <= 8/JM
    return s/absl(sig)

def tail_large(L, y_lo, n, u_up):
    # sum over k with a_k^- > ASW of t_k <= 2 (1 + 2u/pi) beta_1 [x phi(ASW)/ASW + (1/2) int_ASW^inf phi(a)/a da],
    # phi(a) = alpha(a)/(1-e^{-a}) (phi(a)/a decreasing), int <= (1/(ASW(1-e^{-ASW}))) sum_c (5/c) e^{-c ASW/5}.
    x = PI*L/y_lo
    phiA = alpha(ASW)/(1 - (-ASW).exp())
    integ = sum((arb(5)/c*(-arb(c)*ASW/5).exp() for c in range(1, 5)), arb(0))/(ASW*(1 - (-ASW).exp()))
    return 2*(1 + 2*u_up/PI)*arb(2)*(x*phiA/ASW + integ/2)   # beta_1 <= 2

def II_explicit(h, Lb, ylo, yhi, n):
    # sum_{k: a_k^- <= ASW} t_k, all quantities as balls over the box; n exact.
    tot = arb(0); k = 1
    ymid = ball(ylo, yhi)
    while True:
        am_lo = (2*k-1)*PI*arb(Lb.lower())/arb(yhi)          # smallest a_k^- over the box
        if am_lo > ASW: break
        am = (2*k-1)*PI*Lb/ymid; ap = (2*k+1)*PI*Lb/ymid
        am = arb(am.lower()) if am.lower() > 0 else am
        beta = 1 + (-(2*k-1)*PI*n*Lb/ymid).exp()
        G = 1/(1 - (-arb(am.lower())).exp())
        M = PI*n/ymid + 1
        D = arb(ap.upper()).exp().sqrt()*PI*n/(2*ymid)
        r = 2*M*G/D
        lp = r.log() if r.upper() > 1 else arb(0)
        lp = arb(0).union(lp) if False else (arb(lp.upper()) if lp.upper() > 0 else arb(0))
        mn = M*G
        alt = D*(1 + lp)
        m2 = mn if mn.upper() < alt.upper() else alt
        t = ymid/((2*k-1)*PI*n)*beta*alpha(arb(am.lower()))*(2*G + 2*m2)
        tot += arb(t.upper())
        k += 1
    return tot

def R1_bound(h, L, y_lo, sig_abs_up, n):
    # R1 <= B (1 + log(pi n/y)) + (2+4 c1)(-log(1-e^{-L})) + 4 pi c1 |sigma|/y,  B = 2 + 1.6 c1 (h!=0), 2 + 4 c1 (h=0)
    B = (2 + arb('1.6')*C1) if h else (2 + 4*C1)
    lg = (PI*n/y_lo).log()
    return B*(1 + lg) + (2 + 4*C1)*(-(1 - (-L).exp()).log()) + 4*PI*C1*sig_abs_up/y_lo, B

def box_Na(h, Llo, Lhi, ylo, yhi):
    Lb = ball(Llo, Lhi); chi = 1 if Llo <= 0.1 else 0
    sig = sigbox(Llo, Lhi, ylo, yhi)
    Mv = M_sigma(h, sig, Llo)
    g = p_real_point(Lhi) - Mv - arb('0.1')*chi
    R1, B = R1_bound(h, arb(Llo), arb(ylo), abs(sig).upper(), n0)
    II = II_explicit(h, Lb, ylo, yhi, n0) + tail_large(arb(Lhi), arb(ylo), n0, arb(yhi)/n0)
    V = -n0*g + R1 + II + ER0
    ok = V.is_finite() and V.upper() < 0 and g.lower() > (arb('1.5') + B)/n0
    return ok, V, g

def box_Nb(h, Llo, Lhi, ua, ub):
    r1 = box_Nb_opt(h, Llo, Lhi, ua, ub, True)
    if r1[0] or ua == 0: return r1
    r2 = box_Nb_opt(h, Llo, Lhi, ua, ub, False)
    return r2 if r2[0] else r1

def box_Nb_opt(h, Llo, Lhi, ua, ub, use_const):
    # y = u n in [max(Y2, ua n), ub n];  all n >= N0.  Bound written as  -n G + a + b log n + c log^2 n.
    L = ball(Llo, Lhi); Llw = arb(Llo); chi = 1 if Llo <= 0.1 else 0
    u_up = arb(ub); u_lo = arb(ua)
    beta1 = 1 + (-PI*Llw/u_up).exp()                          # beta_k <= 1 + e^{-pi L/u}
    m = mh(h, Llw)
    # n M(sigma) <= n m/|sigma| <= n m/y,  y >= max(Y2, u n)  ->  <= n m/Y2  (and <= m/u_lo if ua > 0)
    Gcoef = p_real_point(Lhi) - arb('0.1')*chi - beta1*u_up**2/Llw
    if ua > 0 and use_const:
        const_M = m/u_lo; linM = arb(0)
    else:
        const_M = arb(0); linM = m/Y2
    Gcoef = Gcoef - linM
    # R1: B (1 + log(pi/u)) + (2+4c1)(-log(1-e^{-L})) + 4 pi c1 (1 + L/Y2);  log(pi/u) <= log(pi n/Y2) if u small
    B = (2 + arb('1.6')*C1) if h else (2 + 4*C1)
    if ua > 0:
        R1a = B*(1 + (PI/u_lo).log()); R1b = arb(0)
    else:
        R1a = B*(1 + (PI/Y2).log()); R1b = B            # coefficient of log n
    R1a = R1a + (2 + 4*C1)*(-(1 - (-Llw).exp()).log()) + 4*PI*C1*(1 + arb(Lhi)/Y2)
    # II small-k, term 1: beta u^2 n/L  (in Gcoef) + (8 beta u/pi)(1 + (1/2) log(2K_sw)),  2K_sw <= ASW u n/(pi L) + 1
    #   log(ASW u n/(pi L) + 1) <= log n + l1,  l1 = log(ASW u/(pi L) + 1/N0)
    l1 = (ASW*u_up/(PI*Llw) + 1/n0).log()
    T1a = 8*beta1*u_up/PI*(1 + l1/2); T1b = 8*beta1*u_up/PI/2
    # term 2: beta psi e^{pi L/Y2} (1 + Lambda)(1 + (1/2) log(2 K_sw)),  Lambda <= lam0 + log n,
    #   lam0 = log(4(1+2u/pi)) + log(u/(pi L) + 1/N0)
    psi = (arb('0.3')*ASW).exp() + (arb('0.1')*ASW).exp() + 2
    pre = beta1*psi*(PI*arb(Lhi)/Y2).exp()
    lam0 = (4*(1 + 2*u_up/PI)).log() + (u_up/(PI*Llw) + 1/n0).log()
    lam0 = lam0 if lam0 > 0 else arb(0)
    # (1 + lam0 + log n)(1 + l1/2 + (1/2) log n)
    P0 = 1 + lam0; Q0 = 1 + (l1/2 if l1 > 0 else arb(0))
    T2a = pre*P0*Q0; T2b = pre*(P0/2 + Q0); T2c = pre/2
    TL = tail_large(arb(Lhi), arb(Y2), n0, u_up)
    a = const_M + R1a + T1a + T2a + TL + arb('6.3')
    b = R1b + T1b + T2b + arb('1.5')
    c = T2c
    lnN = n0.log()
    V = -n0*Gcoef + a + b*lnN + c*lnN**2
    ok = V.is_finite() and V.upper() < 0 and Gcoef.lower() > ((b + 2*c*lnN)/n0).upper()
    return ok, V, Gcoef

def cert(fn, h, a0, a1, b0, b1, wscale, depth, st):
    ok, V, g = fn(h, a0, a1, b0, b1)
    if ok:
        st['boxes'] += 1; st['worstV'] = max(st['worstV'], float(V.upper())); return True
    if depth > 12:
        print('FAIL %s h=%d L=[%.4f,%.4f] [%.5g,%.5g] V=%s g=%s' % (fn.__name__, h, a0, a1, b0, b1, V, g), flush=True)
        return False
    if (a1-a0)/0.05 >= (b1-b0)/wscale:
        m = (a0+a1)/2
        return cert(fn, h, a0, m, b0, b1, wscale, depth+1, st) and cert(fn, h, m, a1, b0, b1, wscale, depth+1, st)
    m = (b0+b1)/2
    return cert(fn, h, a0, a1, b0, m, wscale, depth+1, st) and cert(fn, h, a0, a1, m, b1, wscale, depth+1, st)

if __name__ == '__main__':
    t0 = time.time(); allok = True
    Ls = [0.0999] + [0.1 + 0.1*i for i in range(1, 21)]   # AUDIT: first box starts below the exact 0.1 (float 0.1 is 5.5e-18 too high)
    for h in (1, 0):
        st = {'boxes': 0, 'worstV': -1e300}; ok = True
        ys = [YE[h]]
        while ys[-1] < Y2:
            ys.append(min(Y2, ys[-1]*1.15 if ys[-1] > 20 else ys[-1] + 1.0))
        for i in range(20):
            for j in range(len(ys)-1):
                ok = cert(box_Na, h, Ls[i], Ls[i+1], ys[j], ys[j+1], 1.0, 0, st) and ok
        print('h=%d (N-a) y in [%.0f,%d]: %s boxes=%d worst V(N0)=%.1f (%.0fs)' % (h, YE[h], Y2, 'CERTIFIED' if ok else 'FAILED', st['boxes'], st['worstV'], time.time()-t0), flush=True)
        allok = allok and ok
        st = {'boxes': 0, 'worstV': -1e300}; ok = True
        us = [0.0] + [0.01*k for k in range(1, int(round(float(UMAX.mid())*100))+1)]
        for i in range(20):
            for j in range(len(us)-1):
                ok = cert(box_Nb, h, Ls[i], Ls[i+1], us[j], us[j+1], 0.01, 0, st) and ok
        print('h=%d (N-b) y in [%d, %.2f n]: %s boxes=%d worst V(N0)=%.1f (%.0fs)' % (h, Y2, float(UMAX.mid()), 'CERTIFIED' if ok else 'FAILED', st['boxes'], st['worstV'], time.time()-t0), flush=True)
        allok = allok and ok
    print('N0=%d: piece (N) %s' % (N0, 'CERTIFIED for all n >= N0' if allok else 'FAILED'))
