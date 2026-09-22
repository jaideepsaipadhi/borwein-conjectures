# centre3_peak.py -- certified near-peak piece (A) of the centre regime, n-uniform for n >= N0.
# Ball arithmetic (python-flint acb/arb) over FULL parameter boxes; no midpoint evaluation enters a bound.
#
# Notation: q = zeta^h e^{-sigma/(5n)}, sigma = L - i y, zeta = e^{2 pi i/5}, h = 1..4.
#   Phi_c(s) = log(1 - zeta^{hc} e^{-s})  (principal log; |zeta^{hc}e^{-s}|<=1 for Re s>=0),
#   Q = sum_c Phi_c,  p(sigma) = int_0^1 Q(t sigma) dt,
#   C_h(sigma) = sum_c (c/5 - 1/2)(Phi_c(sigma) - Phi_c(0)).
# Lemma A (finite EM, remainder (B2(c/5)/2)(f'(n)-f'(0)) + int B~2/2 f'', |B~2|<=1/6):
#   |log S_n(q) - n p(sigma) - C_h(sigma)| <= K_h(sigma)/n,
#   K_h = sum_c [ |B2(c/5)|/2 |sigma| (|Phi_c'(sigma)|+|Phi_c'(0)|) + |sigma|^2/12 sup_{0<=t<=1}|Phi_c''(t sigma)| ]
#   valid for |y| < 2pi/5 (no singularity on [0,sigma]).
# Laplace on |y|<=YIN with odd-term cancellation; strip YIN<=|y|<=Y2 by Re psi(y) = -int_0^y (y-u) Re p''(L-iu) du.
# Output per L-box and class s: margin M_s = sgn_s Re sum_h zeta^{-hs} e^{C_h(L)}, error bound
#   ERR = sum_h |e^{C_h}| (rho_h + strip_h), all evaluated at n = N0 and each term non-increasing in n >= N0.
# Normalisation: c(m) = N0f * [ sum_h zeta^{-hm} e^{C_h(L)}(1+rho_h) + strip + OFF ],
#   N0f = e^{n p(L) + 2 n mu L} sqrt(2 pi/(n a))/(10 pi n), a = p''(L);  OFF = region |y_h| > Y2 for all h (piece B, NOT here).
import sys, math
from flint import arb, acb, acb_series, ctx
ctx.prec = 80
N0 = int(sys.argv[1]) if len(sys.argv) > 1 else 1001
LMAX = 2.1; NL = int(sys.argv[2]) if len(sys.argv) > 2 else 210
L_FROM = int(sys.argv[4]) if len(sys.argv) > 4 else 0; L_TO = int(sys.argv[5]) if len(sys.argv) > 5 else NL
YIN, Y2 = 0.35, 1.10
NT = int(sys.argv[3]) if len(sys.argv) > 3 else 64
DEBUG = 0
PI = arb.pi(); I = acb(0, 1)
Z = [acb.exp_pi_i(acb(2*k)/5) for k in range(5)]
def B2(x): return x*x - x + arb(1)/6
def box(Llo, Lhi, ylo, yhi):   # sigma = L - i y
    return acb(arb((Llo+Lhi)/2, (Lhi-Llo)/2), -arb((ylo+yhi)/2, (yhi-ylo)/2))
def tball(i, N): return arb((i+0.5)/N, 0.5/N)
def phi_ser(h, c, s, K):
    ctx.cap = K+1
    x = acb_series([s, 1])
    f = (1 - Z[(h*c) % 5]*(-x).exp()).log()
    co = f.coeffs(); co += [acb(0)]*(K+1-len(co))
    return [co[k]*math.factorial(k) for k in range(K+1)]      # derivatives 0..K
def Q_der(s, K):
    d = [acb(0)]*(K+1)
    for c in range(1, 5):
        e = phi_ser(1, c, s, K); d = [d[k]+e[k] for k in range(K+1)]
    return d
def amax(x): return abs(x).upper() if isinstance(x, acb) else abs(x).upper()
def ub(x):  return arb(x.mid()+x.rad()) if False else x.upper()
def star_sup(fn, sig):  # sup over {t sig : t in [0,1]} of |fn(t sig)|, fn returns acb
    m = arb(0)
    for i in range(NT):
        v = abs(fn(acb(tball(i, NT))*sig)).upper()
        m = m if m > v else v
    return m
res = []
worst = {s: None for s in range(5)}
ybox = []
ny_in = 14; ny_st = 30
for j in range(ny_in):  ybox.append((-YIN + 2*YIN*j/ny_in, -YIN + 2*YIN*(j+1)/ny_in, 'in'))
for j in range(ny_st):
    a_, b_ = YIN + (Y2-YIN)*j/ny_st, YIN + (Y2-YIN)*(j+1)/ny_st
    ybox.append((a_, b_, 'st')); ybox.append((-b_, -a_, 'st'))
# fixed constants Phi_c'(0) per h
phi0p = {(h, c): abs(phi_ser(h, c, acb(0), 1)[1]).upper() for h in (1, 2, 3, 4) for c in range(1, 5)}
for iL in range(L_FROM, L_TO):
    Llo, Lhi = LMAX*iL/NL, LMAX*(iL+1)/NL
    Lb = arb((Llo+Lhi)/2, (Lhi-Llo)/2)
    # a = p''(L) enclosure (real), via int t^2 Q''(tL)
    a = arb(0)
    for i in range(NT):
        t0, t1 = arb(i)/NT, arb(i+1)/NT
        a += (t1**3 - t0**3)/3 * Q_der(acb(tball(i, NT)*Lb), 2)[2].real
    alo = a.lower()
    assert alo > 0
    # sup/inf over window boxes
    K = {h: arb(0) for h in (1, 2, 3, 4)}
    k1 = {h: arb(0) for h in (1, 2, 3, 4)}; k2 = dict(k1)
    k3 = arb(0); k4 = arb(0); bmin_in = None
    beta = {}   # lower bound of Re p'' on each strip/inner box (for |y| in box)
    for (ylo, yhi, kind) in ybox:
        sig = box(Llo, Lhi, ylo, yhi)
        # Re p''(sig) >= sum_i (t1^3-t0^3)/3 * min Re Q''(t sig)
        rp = arb(0); P3 = acb(0); P4 = acb(0)
        for i in range(NT):
            t0, t1 = arb(i)/NT, arb(i+1)/NT
            d = Q_der(acb(tball(i, NT))*sig, 4)
            rp += (t1**3 - t0**3)/3 * d[2].real.lower()
            # p^(k)(sig) = int_0^1 t^k Q^(k)(t sig) dt ; on each t-piece the mean of Q^(k) lies in the (convex) ball
            P3 += (t1**4 - t0**4)/4 * d[3]; P4 += (t1**5 - t0**5)/5 * d[4]
        beta[(ylo, yhi)] = rp.lower()
        if DEBUG and iL==0: print(ylo,yhi,kind,rp)
        if kind == 'in':
            bmin_in = rp.lower() if bmin_in is None else min(bmin_in, rp.lower())
            k3 = max(k3, (abs(P3)/6).upper()); k4 = max(k4, (abs(P4)/24).upper())
        for h in (1, 2, 3, 4):
            Kh = arb(0); C1 = acb(0); C2 = acb(0)
            for c in range(1, 5):
                d = phi_ser(h, c, sig, 2)
                s2 = star_sup(lambda s: phi_ser(h, c, s, 2)[2], sig)
                Kh += abs(B2(arb(c)/5))/2*abs(sig)*(abs(d[1]) + phi0p[(h, c)]) + abs(sig)**2/12*s2
                C1 += (arb(c)/5 - arb(1)/2)*d[1]; C2 += (arb(c)/5 - arb(1)/2)*d[2]
            K[h] = max(K[h], Kh.upper())
            if kind == 'in':
                k1[h] = max(k1[h], abs(C1).upper()); k2[h] = max(k2[h], (abs(C2)/2).upper())
            else:
                k1[h] = max(k1[h], abs(C1).upper())
    assert bmin_in > 0
    b = arb(bmin_in); n = arb(N0)
    # strip contribution (both sides), relative to sqrt(2 pi/(n a)): sqrt(n a/(2pi)) * int e^{n Re psi + Re dC + K/n}
    # Re psi(y) <= -sum over boxes [u0,u1] c [0,|y|] of beta_box * int_{u0}^{u1} (|y|-u) du (beta from boxes covering [-,+]).
    def beta_at(u0, u1):
        vals = [beta[(yl, yh)] for (yl, yh, k) in ybox if (yl <= u0 and yh >= u1) or (yl <= -u1 and yh >= -u0)]
        return min(vals)
    ub_pts = [0.0] + [YIN*(j+1)/(ny_in//2+0.5) for j in range(0)]  # simple: use box edges
    edges = sorted(set([0.0] + [abs(e) for (yl, yh, k) in ybox for e in (yl, yh)]))
    edges = [e for e in edges if e <= Y2 + 1e-12]
    def repsi_ub_box(lo, hi):
        # Re psi(y) = -int_0^|y| (|y|-u) Re p''(L-iu) du  (sign of y absorbed: beta_at takes min over both signs)
        # <= -sum_{pieces [u0,u1] in [0,lo]} beta_i (u1-u0)(Y-(u0+u1)/2) + max(0,-beta_box)(hi-lo)^2/2, Y in [lo,hi];
        # linear in Y -> sup at an endpoint.
        best = None
        for Y in (lo, hi):
            tot = arb(0); prev = 0.0
            for e in edges[1:]:
                if e > lo + 1e-12: break
                bb = arb(beta_at(prev, e))
                tot -= bb*(arb(e)-arb(prev))*(arb(Y)-(arb(prev)+arb(e))/2)
                prev = e
            bl = beta_at(lo, hi)
            tot += arb(max(0, -bl))*(arb(hi)-arb(lo))**2/2
            best = tot if best is None else (tot if tot.upper() > best.upper() else best)
        return best
    cR = {h: k1[h]*arb(Y2) for h in (1, 2, 3, 4)}
    strip = {}
    for h in (1, 2, 3, 4):
        st = arb(0)
        for (ylo, yhi, kind) in ybox:
            if kind != 'st': continue
            lo = min(abs(ylo), abs(yhi)); w = arb(yhi - ylo)
            D = arb((-repsi_ub_box(lo, max(abs(ylo), abs(yhi)))).lower())   # Re psi <= -D on this box
            assert D > 0
            # f(n)=sqrt(n) e^{-nD} non-increasing for n >= 1/(2D)
            assert n >= 1/(2*D)
            st += (n*a.upper()/(2*PI)).sqrt()*w*(-n*D + cR[h] + K[h]/n).exp()
        strip[h] = st.upper()
    # rho_h: inner Laplace error relative to sqrt(2pi/(n a)), a in ball -> use a (upper/lower) conservatively
    A = arb(alo); Aup = a.upper()
    rho = {}
    for h in (1, 2, 3, 4):
        cRi = k1[h]*arb(YIN)
        t_eps = ((K[h]/n).exp() - 1)*cRi.exp()*(arb(Aup)/b).sqrt()
        t_zt = 3*k4/(n*A*A) + k2[h]/(n*A)
        def zz(cc, pref):  # (1/2) int |Z|^2 e^{-n cc y^2/2} / sqrt(2pi/(n a)), |Z|^2 <= 2 n^2 k3^2 y^6 + 2 k1^2 y^2
            return pref*(arb(Aup)/cc).sqrt()*(15*k3*k3/(n*cc**3) + cRi*0 + k1[h]**2/(n*cc))
        t_z = zz(A, arb(1)) + zz(b, cRi.exp())
        t_tail = 2*(-n*A*YIN**2/2).exp()/(n*A*YIN)*(n*Aup/(2*PI)).sqrt()
        rho[h] = (t_eps + t_zt + t_z + t_tail).upper()
        if h == 1: parts = (float(t_eps.mid()), float(t_zt.mid()), float(t_z.mid()), float(t_tail.mid()), float(k3.mid()), float(k1[h].mid()))
    # margins
    eC = {}
    for h in (1, 2, 3, 4):
        Ch = acb(0)
        for c in range(1, 5):
            Ch += (arb(c)/5 - arb(1)/2)*(phi_ser(h, c, acb(Lb), 0)[0] - phi_ser(h, c, acb(0), 0)[0])
        eC[h] = Ch.exp()
    err = sum(abs(eC[h]).upper()*(rho[h] + strip[h]) for h in (1, 2, 3, 4))
    line = []
    for s in range(5):
        M = sum(Z[(-h*s) % 5]*eC[h] for h in (1, 2, 3, 4)).real*(1 if s == 0 else -1)
        ok = M.lower() > err.upper()
        r = (err/M).upper() if M > 0 else None
        line.append((s, M.lower(), ok, r))
        if worst[s] is None or not ok or (r is not None and r > worst[s][1]):
            worst[s] = (Llo, r, ok) if (worst[s] is None or not worst[s][2] is False) else worst[s]
    print('L=[%.3f,%.3f] a=%.4f b_in=%.4f K=%.3f rho1=%.4f rho2=%.4f strip=%.2e err=%.4f  M_4=%.4f %s %s' % (
        Llo, Lhi, float(a.mid()), float(b.mid()), float(max(K.values()).mid()), float(rho[1].mid()), float(rho[2].mid()),
        float(max(strip.values()).mid()), float(err.mid()), float(line[4][1].mid()),
        str(['%.3g'%x for x in parts]), 'ALL_OK' if all(x[2] for x in line) else 'FAIL ' + str([x[0] for x in line if not x[2]])), flush=True)
    res.append(line)
print('N0=%d  boxes=%d  all classes certified on every box: %s' % (N0, NL, all(all(x[2] for x in l) for l in res)))
for s in range(5):
    rs = [l[s][3] for l in res if l[s][3] is not None]
    print('class %d: worst err/margin = %.4f' % (s, max(float(x.mid()) for x in rs) if rs else float('nan')))
