# offpeak_chunk_arb.py -- piece (F) of the off-peak lemma (B), with every bound in arb ball arithmetic.
#
# Same certificate as offpeak_chunk.py (chunk lemma + Hadamard convexity in s + maximum principle in
# thetat = theta + i mu on rectangles, arb assembly), but the cell evaluator is arb:
#   for a cell (theta-interval, V-arc alpha-interval, mu-interval) every factor is bounded by
#   |1 - rho e^{i phi}|^2 = 1 - 2 rho cos(phi) + rho^2 <= max_{rho in {rho_lo, rho_hi}} (1 - 2 rho c_lo + rho^2),
#   phi = alpha + k' theta (arb ball over the cell), c_lo = lower end of the arb enclosure of cos(phi),
#   rho_lo/hi = ends of the arb enclosure of e^{-s - k' mu}; the product (arb) is compared to e^{2t} (arb).
# Cells are addressed by integer indices (dyadic subdivision of the rectangle, of [0, 2pi], of the mu range), and
# the arb balls are built from exact endpoints by union, so the leaves tile the parameter set exactly.
# The float64 checker of offpeak_chunk.py is used ONLY as a search heuristic (to decide where to split and which
# targets to try); a cell counts as certified only when the arb check passes. No float enters any bound.
# Other changes vs offpeak_chunk.py: the last rectangle ends at a float >= pi (float pi < pi); d_min in the
# assembly is chosen with exact arb comparisons.
#
# Usage:
#   python3 offpeak_chunk_arb.py N0 PIECE [RLO RHI] [-j JOBS]      PIECE in 0..7 (see PIECES), rects RLO..RHI-1
#   python3 offpeak_chunk_arb.py list N0                             list pieces and their rect counts
#   python3 offpeak_chunk_arb.py summarize LOG...                    check that the logs cover every rect with OK
import sys, time, math
import numpy as np
from flint import arb, acb, ctx
ctx.prec = 64
import offpeak_chunk as OC
from offpeak_chunk import Chunk, SG, Qa, dQa, LOG5, f_up

PI_A = arb.pi()
TP = 2*math.pi/5; D0 = 0.0499
PI_UP = math.nextafter(math.pi, 4.0)
assert arb(PI_UP) > PI_A
PIECES = [(40, 0.0499, 0.1), (40, TP-0.1, TP-D0), (40, TP+D0, TP+0.1), (40, 2*TP-0.1, 2*TP-D0), (40, 2*TP+D0, 2*TP+0.1),
          (20, 0.1, TP-0.1), (20, TP+0.1, 2*TP-0.1), (20, 2*TP+0.1, PI_UP)]

def rects_of(piece):
    T, lo, hi = PIECES[piece]; width = 0.01 if T == 40 else 0.02
    out = []; x = lo
    while x < hi - 1e-15:
        y = min(hi, x + width); out.append((x, y)); x = y
    return out

def cos_lo(ph):
    # rigorous lower bound for cos on the ball ph = m +- r (r = ph.rad(), exact mid m).  For |h| <= r:
    #   cos(m+h) = cos m cos h - sin m sin h >= min(cos m, cos m (1 - r^2/2)) - |sin m| r
    # (cos h in [1 - r^2/2, 1], |sin h| <= r).  Second order in r, unlike arb's generic cos(ball) >= cos m - r.
    # The result is max(this, arb's own bound): both are valid lower bounds.
    m = arb(ph.mid()); r = arb(ph.rad())
    cm = m.cos(); sm = m.sin()
    t1 = cm; t2 = cm*(1 - r*r/2)
    lo1 = min(t1.lower(), t2.lower())
    b = arb(lo1) - abs(sm)*r
    return arb(max(b.lower(), ph.cos().lower()))

class ArbChunk:
    def __init__(self, T, N0):
        self.F = Chunk(T, N0)            # float heuristic + estimate
        self.T = T; self.kp = [int(k) for k in self.F.kp]; self.K = len(self.kp)
        self.mu_up = self.F.mu_up        # float >= 2.1/(5 N0), exact dyadic
        assert arb(self.mu_up) >= self.F.MU
        self.rc = {}
        self.stats = {'arb': 0, 'arbfail': 0}

    def rho_arb(self, s, mlo, mhi, key):
        if key not in self.rc:
            sa = arb(s); out = []
            for k in self.kp:
                a = (-(sa + k*mlo)).exp(); b = (-(sa + k*mhi)).exp()   # a >= b
                out.append((arb(b.lower()), arb(a.upper())) if (mlo != mhi) else (arb(a.lower()), arb(a.upper())))
            self.rc[key] = out
        return self.rc[key]

    def arb_ok(self, thb, alb, rhos, thr):
        p = arb(1)
        for k, (rl, rh) in zip(self.kp, rhos):
            c = cos_lo(alb + k*thb)
            x1 = 1 - 2*rl*c + rl*rl
            if rl is rh or rl == rh:
                x = arb(x1.upper())
            else:
                x2 = 1 - 2*rh*c + rh*rh
                x = arb(max(x1.upper(), x2.upper()))
            p = p*x   # no early exit: later factors can be < 1, so a partial product > thr proves nothing
        return bool(p < thr)

    def bb(self, th_a, th_b, mu_lo, mu_hi, s, logthr, n_al0=64, maxcells=4_000_000):
        width = th_b - th_a
        n_th0 = max(1, int(math.ceil(width/(0.002/self.F.kmax*50)))) if width > 0 else 1
        thr = (2*arb(logthr)).exp()
        A_th, B_th = arb(th_a), arb(th_b); W_th = B_th - A_th
        A_mu, B_mu = arb(mu_lo), arb(mu_hi)
        # integer cell addresses
        IT = np.repeat(np.arange(n_th0), n_al0); LT = np.zeros_like(IT)
        IA = np.tile(np.arange(n_al0), n_th0); LA = np.zeros_like(IA)
        IM = np.zeros_like(IT); LM = np.zeros_like(IT)
        total = 0
        while len(IT):
            if total > maxcells: return False, total
            NT = n_th0*(2.0**LT); NA = n_al0*(2.0**LA)
            DL = width/(2*NT); TH = th_a + width*(IT + 0.5)/NT
            ET = math.pi/NA; AL = 2*math.pi*(IA + 0.5)/NA
            fmask = np.zeros(len(IT), dtype=bool)
            keys = set(zip(IM.tolist(), LM.tolist()))
            for (im, lm) in keys:
                idx = np.where((IM == im) & (LM == lm))[0]
                ml = mu_lo + (mu_hi - mu_lo)*im/2.0**lm; mh = mu_lo + (mu_hi - mu_lo)*(im+1)/2.0**lm
                rlo, rhi = self.F.rho_bounds(s, ml, mh) if mu_lo != mu_hi else self.F.rho_bounds(s, mu_lo, mu_lo)
                # heuristic prefilter (float) -- margins widened slightly so it is optimistic, never decisive
                okf = self.F.check(TH[idx], DL[idx], AL[idx], ET[idx], rlo, rhi, logthr)
                # rigorous arb check on the float-passing cells
                if mu_lo == mu_hi:
                    mla = mha = A_mu
                else:
                    den = 2**int(lm)
                    mla = A_mu + (B_mu - A_mu)*im/den; mha = A_mu + (B_mu - A_mu)*(im+1)/den
                rhos = self.rho_arb(s, mla, mha, (s, mu_lo, mu_hi, im, lm)) if mu_lo != mu_hi else self.rho_arb(s, A_mu, A_mu, (s, mu_lo, mu_hi))
                for jj in idx[okf]:
                    nt = n_th0*2**int(LT[jj]); na = n_al0*2**int(LA[jj])
                    if width > 0:
                        thb = (A_th + W_th*int(IT[jj])/nt).union(A_th + W_th*(int(IT[jj])+1)/nt)
                    else:
                        thb = A_th
                    alb = (2*PI_A*int(IA[jj])/na).union(2*PI_A*(int(IA[jj])+1)/na)
                    self.stats['arb'] += 1
                    if self.arb_ok(thb, alb, rhos, thr):
                        fmask[jj] = True
                    else:
                        self.stats['arbfail'] += 1
                total += len(idx)
            f = np.where(~fmask)[0]
            if len(f) == 0: return True, total
            IT, LT, IA, LA, IM, LM = IT[f], LT[f], IA[f], LA[f], IM[f], LM[f]
            DL = DL[f]; ET = ET[f] if np.ndim(ET) else ET
            if (np.asarray(ET) < 1e-7).any() and (DL < 1e-11).any(): return False, total
            if (LT > 40).any() or (LA > 40).any() or (LM > 40).any(): return False, total
            kmax = self.F.kmax
            cth = kmax*DL; cal = np.asarray(ET)*np.ones(len(IT)); cmu = (mu_hi - mu_lo)/2.0**LM*kmax*2
            s_th = (cth >= cal) & (cth >= cmu); s_mu = (~s_th) & (cmu > cal); s_al = ~(s_th | s_mu)
            nI = []
            for sel, which in ((s_th, 0), (s_al, 1), (s_mu, 2)):
                for b in (0, 1):
                    it, lt, ia, la, im, lm = IT[sel].copy(), LT[sel].copy(), IA[sel].copy(), LA[sel].copy(), IM[sel].copy(), LM[sel].copy()
                    if which == 0: it = 2*it + b; lt = lt + 1
                    elif which == 1: ia = 2*ia + b; la = la + 1
                    else: im = 2*im + b; lm = lm + 1
                    nI.append((it, lt, ia, la, im, lm))
            IT, LT, IA, LA, IM, LM = [np.concatenate([v[i] for v in nI]) for i in range(6)]
        return True, total

def certify_rect(C, th_a, th_b, side_cache, tol0=0.25, rel=None):
    t = []; cells = 0
    for s in SG:
        e = max(C.F.estimate(th_a, th_b, s, 0.0), C.F.estimate(th_a, th_b, s, C.mu_up))
        tol = tol0 if rel is None else max(tol0, rel*float((C.T*Qa(s)).mid()) - rel*e)
        while True:
            tgt = e + tol
            good = True
            for (m0, m1) in ((0.0, 0.0), (C.mu_up, C.mu_up)):
                okk, nc = C.bb(th_a, th_b, m0, m1, s, tgt); cells += nc
                if not okk: good = False; break
            if good:
                for key in (th_a, th_b):
                    if (key, s) in side_cache and side_cache[(key, s)] <= tgt: continue
                    okk, nc = C.bb(key, key, 0.0, C.mu_up, s, tgt); cells += nc
                    if not okk: good = False; break
                    side_cache[(key, s)] = tgt
            if good: break
            tol *= 2
            if tol > 8: return None, cells
        t.append(tgt)
    return t, cells

def assemble(T, N0, t):
    # = offpeak_chunk.assemble, with d_min chosen by exact arb comparison
    n0 = arb(N0); Tn = arb(T); d = []
    for i in range(len(SG)-1):
        D = arb(SG[i+1]) - arb(SG[i]); q1 = Qa(SG[i+1]); dq = abs(dQa(SG[i+1]))
        a1 = arb((Tn*q1 - arb(t[i+1])).lower()); a2 = arb((Tn*(q1 + dq*D) - arb(t[i])).lower())
        d.append(a1 if a1 < a2 else a2)
    E = 1.5*n0.log() + arb('6.3'); ok = True; worst = None
    Lg = [i/100 for i in range(211)]; shrink = 1 - Tn/n0
    for j in range(len(Lg)-1):
        La, Lb = arb(Lg[j]), arb(Lg[j+1])
        common = Tn*(Tn*Lb/n0 + (LOG5 - Qa(Lb))/2) + 4*(Tn-1)*(1 + (-(La*shrink)).exp()).log() + E
        best = None
        rel = [d[i] for i in range(len(d)) if SG[i] < Lg[j+1]]
        dmin = rel[0]
        for x in rel[1:]:
            if x < dmin: dmin = x
        if dmin > 0:
            v1 = (n0/Tn - 1)*dmin - common
            if (dmin/Tn - arb('1.5')/n0) > 0 and v1 > 0: best = v1
        if all(d[i] > 0 for i in range(len(d)) if SG[i] < Lg[j+1]):
            X = La*shrink; W = arb(0); run = None
            for i in range(len(d)):
                run = d[i] if run is None or d[i] < run else run
                lo = arb(SG[i]); hi = arb(SG[i+1])
                if lo >= X: break
                W += run*((hi if hi < X else X) - lo)
            if W > 0:
                v2 = n0/(Tn*Lb)*W - common
                if (W/(Tn*Lb) - arb('1.5')/n0) > 0 and v2 > 0:
                    if best is None or v2.lower() > best.lower(): best = v2
        if best is None or not (best > 0):
            ok = False; continue
        worst = float(best.lower()) if worst is None else min(worst, float(best.lower()))
    return ok, worst

def do_rect(args):
    N0, piece, idx, a, b = args
    T = PIECES[piece][0]; C = ArbChunk(T, N0); side_cache = {}
    t0 = time.time(); out = []; stack = [(a, b)]; allok = True; worst = 1e300
    while stack:
        x, y = stack.pop(); ok = False
        for rel in (0.3, None):
            t, cells = certify_rect(C, x, y, side_cache, rel=rel)
            if t is not None:
                ok, w = assemble(T, N0, t)
                if ok: break
        if ok:
            worst = min(worst, w)
            out.append('  sub [%.17g,%.17g] OK margin(N0)=%.2f t0=%.3f cells=%d' % (x, y, w, t[0], cells))
        elif y - x > 1e-4:
            m = (x + y)/2; stack.append((m, y)); stack.append((x, m))
        else:
            allok = False; out.append('  sub [%.17g,%.17g] FAILED' % (x, y))
    line = 'piece=%d rect=%d T=%d [%.17g,%.17g] %s min margin(N0)=%.2f arbcells=%d arbfail=%d (%.0fs)' % (
        piece, idx, T, a, b, 'OK' if allok else 'FAILED', worst, C.stats['arb'], C.stats['arbfail'], time.time()-t0)
    return line, out, allok

def summarize(logs):
    seen = {}
    for fn in logs:
        for l in open(fn):
            if l.startswith('piece='):
                f = l.split(); p = int(f[0][6:]); r = int(f[1][5:])
                seen[(p, r)] = seen.get((p, r), False) or (' OK ' in l)
    allok = True
    for p in range(len(PIECES)):
        n = len(rects_of(p)); missing = [r for r in range(n) if not seen.get((p, r), False)]
        print('piece %d (T=%d, [%.6f,%.6f]): %d rects, %s' % (p, PIECES[p][0], PIECES[p][1], PIECES[p][2], n,
              'all OK' if not missing else 'MISSING/FAILED %s' % missing))
        allok = allok and not missing
    print('(F) arb:', 'CERTIFIED for all n >= N0 (all rects OK)' if allok else 'INCOMPLETE')

if __name__ == '__main__':
    av = sys.argv[1:]
    if av[0] == 'summarize': summarize(av[1:]); sys.exit()
    if av[0] == 'list':
        for p in range(len(PIECES)): print(p, PIECES[p], len(rects_of(p)))
        sys.exit()
    jobs = 1
    if '-j' in av: i = av.index('-j'); jobs = int(av[i+1]); av = av[:i] + av[i+2:]
    N0 = int(av[0]); piece = int(av[1]); R = rects_of(piece)
    lo = int(av[2]) if len(av) > 2 else 0; hi = int(av[3]) if len(av) > 3 else len(R)
    tasks = [(N0, piece, i, R[i][0], R[i][1]) for i in range(lo, min(hi, len(R)))]
    t0 = time.time()
    if jobs > 1:
        from multiprocessing import Pool
        with Pool(jobs) as P:
            for line, out, ok in P.imap(do_rect, tasks):
                print(line, flush=True); [print(o, flush=True) for o in out]
    else:
        for tk in tasks:
            line, out, ok = do_rect(tk)
            print(line, flush=True); [print(o, flush=True) for o in out]
    print('done piece=%d rects %d..%d (%.0fs)' % (piece, lo, min(hi, len(R))-1, time.time()-t0), flush=True)
