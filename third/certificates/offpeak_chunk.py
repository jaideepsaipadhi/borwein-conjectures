# offpeak_chunk.py -- piece (F) of the off-peak lemma (B): theta at distance >= 0.05 from 0, 2pi/5, 4pi/5
# (theta in [0, pi]; [pi, 2pi] by conjugation).  n-uniform chunk / maximum-principle bound.
#
# Chunk lemma (GAP_OFFPEAK.md sec. 5).  T fixed, M = floor(n/T), T' = n - M T, mu = L/(5n) in (0, MU_MAX],
# MU_MAX = 2.1/(5 N0), lambda = 5 T mu = T L/n.  P(V) = prod_{k' <= 5T, 5 nmid k'} (1 - V e^{i k' thetat}),
# thetat = theta + i mu (the radial drift r^{k'} = e^{-k' mu} is an imaginary shift of theta).
#   log|S_n| <= sum_{m<M} B(m lambda) + 4(T-1) log(1 + e^{-(L-lambda)}),  B(s) = max_{|V|=e^{-s}} log|P(V)|.
#   n p(L) >= sum_{m<M} T Q(m lambda) - T[lambda + (log5 - Q(L))/2]     (Q convex decreasing, |Q'| <= 2).
# B(s; thetat) is (i) convex non-increasing in s (Hadamard), (ii) subharmonic in thetat (sup over V of log|entire|),
# so its sup over a rectangle R = [th_a, th_b] x [0, MU_MAX] is attained on the boundary: bottom edge (mu = 0),
# top edge (mu = MU_MAX) and the two sides.  We certify B(s_i) <= t_i(R) on the boundary for an s-grid by
# branch-and-bound over cells (theta-interval, V-arc, mu-interval) in float64 interval arithmetic that uses ONLY
# IEEE +,-,*,/ (correctly rounded) with explicit rounding margins; all transcendental constants (e^{i theta_c},
# e^{i alpha_c}, rho = e^{-s}, e^{-k' mu}, thresholds e^{2t}) come from arb.  No libm function enters a bound.
import sys, time, math
import numpy as np
from flint import arb, acb, ctx
ctx.prec = 80
U = 2.0**-53

def f_down(x):   # float <= arb x
    v = float(x.lower()); return v - abs(v)*4*U - 1e-300
def f_up(x):
    v = float(x.upper()); return v + abs(v)*4*U + 1e-300
def eix(x):      # float complex approx of e^{i x} for float x, |err| <= 2u
    z = acb(0, arb(x)).exp()
    return complex(float(z.real.mid()), float(z.imag.mid()))

class Chunk:
    def __init__(self, T, N0):
        self.T = T; self.N0 = N0
        self.kp = np.array([k for k in range(1, 5*T+1) if k % 5], dtype=np.int64)
        self.K = len(self.kp); self.kmax = 5*T
        self.MU = arb('2.1')/(5*N0)
        self.mu_up = f_up(self.MU)

    def rho_bounds(self, s, mu_lo, mu_hi):
        # per-factor radius interval [e^{-s-k' mu_hi}, e^{-s-k' mu_lo}], mu_* python floats (exact) or arb
        lo = np.empty(self.K); hi = np.empty(self.K)
        for i, k in enumerate(self.kp):
            lo[i] = f_down((-(arb(s) + k*arb(mu_hi))).exp())
            hi[i] = f_up((-(arb(s) + k*arb(mu_lo))).exp())
        return lo, np.minimum(hi, 1.0)

    def check(self, thc, dl, alc, et, rlo, rhi, logthr):
        """thc, dl, alc, et: float arrays (N,) of cell centres / half-widths (theta, alpha).
        rlo, rhi: (K,) radius bounds.  Returns boolean array: True if (1/2) log prod_k |1 - rho e^{i phi}|^2 <= logthr
        certified for every point of the cell.  Uses only IEEE +,-,*,/ plus arb-supplied constants."""
        N = len(thc)
        # base unit complex numbers from arb (|err| <= 2u each)
        zt = np.array([eix(x) for x in thc]); za = np.array([eix(x) for x in alc])
        # phases z_k = e^{i alpha} (e^{i theta})^{k'}; iterative powers, error <= (6k'+10) u  (<= 1.5e-13 for k'<=200)
        P = np.ones((N,), dtype=complex); zs = np.empty((N, self.K), dtype=complex)
        j = 0
        for k in range(1, self.kmax+1):
            P = P*zt
            if k % 5:
                zs[:, j] = za*P; j += 1
        c = zs.real; sn = np.abs(zs.imag)
        w = (et*(1 + 1e-12) + 1e-14)[:, None] + self.kp[None, :]*(dl*(1 + 1e-12) + 1e-15)[:, None] + 1e-12
        cm = np.where(c >= 0, c*(1 - w*w/2) - sn*w, c - sn*w)
        bad = (w >= 1.5707) | ((c <= 0) & (sn <= w))
        cm = np.where(bad, -1.0, cm) - 1e-12
        cm = np.maximum(cm, -1.0)
        # x = 1 - 2 rho cm + rho^2, convex in rho: max over endpoints; +32u absolute rounding margin
        x1 = 1 - 2*rlo[None, :]*cm + rlo[None, :]*rlo[None, :]
        x2 = 1 - 2*rhi[None, :]*cm + rhi[None, :]*rhi[None, :]
        x = np.maximum(x1, x2) + 32*U
        x = np.maximum(x, 1e-8)
        # product in groups of 16 with frexp bookkeeping
        G = (self.K + 15)//16
        mant = np.ones(N); expo = np.zeros(N, dtype=np.int64)
        for g in range(G):
            gp = np.prod(x[:, 16*g:16*(g+1)], axis=1)
            m_, e_ = np.frexp(gp)
            mant = mant*m_; expo += e_
            m2, e2 = np.frexp(mant); mant = m2; expo += e2
        nops = 2*(self.K + 3*G + 8)
        mant = mant*(1 + 2*U)**nops
        # threshold: prod <= e^{2 logthr}
        thr = (2*arb(logthr)).exp()
        tm, te = math.frexp(f_down(thr))
        rhs = np.ldexp(np.full(N, tm), (te - expo).clip(-1070, 1070))
        return mant <= rhs

    def estimate(self, th_a, th_b, s, mu, nth=None, nal=512):
        # non-rigorous float estimate of max over theta in [th_a,th_b], |V| = e^{-s}, drift mu
        nth = nth or max(3, int((th_b-th_a)/2e-4))
        th = np.linspace(th_a, th_b, nth); al = 2*np.pi*np.arange(nal)/nal
        rho = math.exp(-s)*np.exp(-self.kp*mu)
        best = -1e300
        for i0 in range(0, nth, 64):
            t = th[i0:i0+64]
            ph = al[None, :, None] + self.kp[None, None, :]*t[:, None, None]
            v = 0.5*np.log(1 - 2*rho*np.cos(ph) + rho*rho).sum(2).max()
            best = max(best, v)
        return best

    def bb(self, th_a, th_b, mu_lo, mu_hi, s, logthr, n_th0=None, n_al0=64, maxcells=4_000_000):
        """Branch and bound: certify B(s) <= logthr for theta in [th_a, th_b], mu in [mu_lo, mu_hi]
        (mu interval is split too).  Returns (ok, ncells_checked)."""
        width = th_b - th_a
        n_th0 = n_th0 or max(1, int(math.ceil(width/(0.002/self.kmax*50))))
        # cells: arrays
        dl0 = width/(2*n_th0)
        thc = th_a + dl0*(2*np.arange(n_th0) + 1)
        et0 = math.pi/n_al0
        alc = et0*(2*np.arange(n_al0) + 1)
        TH = np.repeat(thc, n_al0); AL = np.tile(alc, n_th0)
        DL = np.full(TH.shape, dl0 if width > 0 else 0.0); ET = np.full(TH.shape, et0)
        MUL = np.full(TH.shape, float(mu_lo)); MUH = np.full(TH.shape, float(mu_hi))
        total = 0
        rcache = {}
        while len(TH):
            if total > maxcells:
                return False, total
            # group by (mu_lo, mu_hi)
            keys = np.unique(np.stack([MUL, MUH], 1), axis=0)
            fails = []
            for (a, b) in keys:
                idx = np.where((MUL == a) & (MUH == b))[0]
                key = (a, b)
                if key not in rcache:
                    rcache[key] = self.rho_bounds(s, a, b)
                rlo, rhi = rcache[key]
                for i0 in range(0, len(idx), 4096):
                    ii = idx[i0:i0+4096]
                    ok = self.check(TH[ii], DL[ii], AL[ii], ET[ii], rlo, rhi, logthr)
                    fails.append(ii[~ok])
                total += len(idx)
            f = np.concatenate(fails) if fails else np.array([], dtype=np.int64)
            if len(f) == 0:
                return True, total
            TH, DL, AL, ET, MUL, MUH = TH[f], DL[f], AL[f], ET[f], MUL[f], MUH[f]
            if (ET < 1e-7).any() and (DL < 1e-11).any():
                return False, total
            # split rule: largest contribution to phase width / radius spread
            cth = self.kmax*DL; cal = ET; cmu = (MUH - MUL)*self.kmax*2
            s_th = (cth >= cal) & (cth >= cmu)
            s_mu = (~s_th) & (cmu > cal)
            s_al = ~(s_th | s_mu)
            nTH, nDL, nAL, nET, nMUL, nMUH = [], [], [], [], [], []
            # theta split
            i = s_th
            for sg in (-1, 1):
                nTH.append(TH[i] + sg*DL[i]/2); nDL.append(DL[i]/2); nAL.append(AL[i]); nET.append(ET[i]); nMUL.append(MUL[i]); nMUH.append(MUH[i])
            i = s_al
            for sg in (-1, 1):
                nTH.append(TH[i]); nDL.append(DL[i]); nAL.append(AL[i] + sg*ET[i]/2); nET.append(ET[i]/2); nMUL.append(MUL[i]); nMUH.append(MUH[i])
            i = s_mu
            mid = (MUL[i] + MUH[i])/2
            nTH += [TH[i], TH[i]]; nDL += [DL[i], DL[i]]; nAL += [AL[i], AL[i]]; nET += [ET[i], ET[i]]
            nMUL += [MUL[i], mid]; nMUH += [mid, MUH[i]]
            TH, DL, AL, ET, MUL, MUH = [np.concatenate(v) for v in (nTH, nDL, nAL, nET, nMUL, nMUH)]
        return True, total

# ------------------------------------------------------------------------------------------------------------
# Assembly (arb): from certified t_i >= sup_R B(s_i) to  n p(L) - log|S_n| >= 1.5 log n + 6.3  for all n >= N0,
# all L in (0, 2.1], all theta in R.
SG = [i/10 for i in range(22)]          # s-grid 0, 0.1, ..., 2.1
LOG5 = arb(5).log()
def Qa(s):
    s = arb(s)
    if s == 0: return LOG5
    return ((1 - (-5*s).exp())/(1 - (-s).exp())).log()
def dQa(s):
    s = arb(s)
    if s == 0: return arb(-2)
    return 5*(-5*s).exp()/(1 - (-5*s).exp()) - (-s).exp()/(1 - (-s).exp())

def assemble(T, N0, t):
    """t: list of certified upper bounds for B(s_i), i over SG.  Returns (ok, min margin at N0, details)."""
    n0 = arb(N0); Tn = arb(T)
    d = []
    for i in range(len(SG)-1):
        D = arb(SG[i+1]) - arb(SG[i])
        q1 = Qa(SG[i+1]); dq = abs(dQa(SG[i+1]))
        d.append(min((Tn*q1 - arb(t[i+1])).lower(), (Tn*(q1 + dq*D) - arb(t[i])).lower()))
    d = [arb(x) for x in d]
    E = 1.5*n0.log() + arb('6.3')
    worst = None; ok = True
    Lg = [i/100 for i in range(211)]
    shrink = 1 - Tn/n0
    for j in range(len(Lg)-1):
        La, Lb = arb(Lg[j]), arb(Lg[j+1])
        common = Tn*(Tn*Lb/n0 + (LOG5 - Qa(Lb))/2) + 4*(Tn-1)*(1 + (-(La*shrink)).exp()).log() + E
        best = None
        # bound 1: all chunks in cells with s_i < Lb
        rel = [d[i] for i in range(len(d)) if SG[i] < Lg[j+1]]
        dmin = min(rel, key=lambda x: float(x.lower()))
        if dmin > 0:
            v1 = (n0/Tn - 1)*dmin - common
            if (dmin/Tn - arb('1.5')/n0) > 0:
                best = v1
        # bound 2: running-min step function dt(s) = min_{s'<=s} d(s') (non-increasing), left sums >= integral:
        #   sum_{m<M} D0(m lam) >= (1/lam) int_0^{M lam} dt >= (n/(T Lb)) int_0^{La(1-T/N0)} dt
        if all(d[i] > 0 for i in range(len(d)) if SG[i] < Lg[j+1]):
            X = La*shrink; W = arb(0); run = None
            for i in range(len(d)):
                run = d[i] if run is None or d[i] < run else run
                lo = arb(SG[i]); hi = arb(SG[i+1])
                if lo >= X: break
                seg = (hi if hi < X else X) - lo
                W += run*seg
            if W > 0:
                v2 = n0/(Tn*Lb)*W - common
                if (W/(Tn*Lb) - arb('1.5')/n0) > 0:
                    if best is None or v2.lower() > best.lower(): best = v2
        if best is None or not (best > 0):
            ok = False
            worst = -1e9 if best is None else min(worst if worst is not None else 1e9, float(best.lower()))
            continue
        worst = float(best.lower()) if worst is None else min(worst, float(best.lower()))
    return ok, worst, [float(x.lower()) for x in d]

def certify_rect(C, th_a, th_b, side_cache, tol0=0.25, log=None, rel=None):
    t = []; cells = 0
    for s in SG:
        e = max(C.estimate(th_a, th_b, s, 0.0), C.estimate(th_a, th_b, s, C.mu_up))
        tol = tol0 if rel is None else max(tol0, rel*float((C.T*Qa(s)).mid()) - rel*e)
        while True:
            tgt = e + tol
            parts = [(th_a, th_b, 0.0, 0.0), (th_a, th_b, C.mu_up, C.mu_up)]
            good = True
            for (a, b, m0, m1) in parts:
                okk, nc = C.bb(a, b, m0, m1, s, tgt); cells += nc
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

def run_piece(T, N0, th_lo, th_hi, width, logf):
    C = Chunk(T, N0); side_cache = {}
    stack = []
    x = th_lo
    while x < th_hi - 1e-15:
        y = min(th_hi, x + width); stack.append((x, y)); x = y
    stack.reverse()
    allok = True; worst = 1e300; nrect = 0; t0 = time.time()
    while stack:
        a, b = stack.pop()
        ok = False
        for rel in (0.3, None):
            t, cells = certify_rect(C, a, b, side_cache, rel=rel)
            if t is not None:
                ok, w, d = assemble(T, N0, t)
                if ok: break
        if ok:
            nrect += 1; worst = min(worst, w)
            print('T=%d rect [%.6f,%.6f] OK  min margin(N0)=%.2f  t0=%.3f cells=%d (%.0fs)' % (T, a, b, w, t[0], cells, time.time()-t0), file=logf, flush=True)
        elif b - a > 1e-4:
            print('T=%d rect [%.6f,%.6f] split (t=%s)' % (T, a, b, None if t is None else '%.3f' % t[0]), file=logf, flush=True)
            m = (a + b)/2; stack.append((m, b)); stack.append((a, m))
        else:
            print('T=%d rect [%.6f,%.6f] FAILED' % (T, a, b), file=logf, flush=True)
            allok = False
    return allok, worst, nrect

if __name__ == '__main__':
    N0 = int(sys.argv[1]) if len(sys.argv) > 1 else 821
    which = sys.argv[2] if len(sys.argv) > 2 else 'all'
    TP = 2*math.pi/5; D0 = 0.0499
    pieces = {'T40': [(0.0499, 0.1), (TP-0.1, TP-D0), (TP+D0, TP+0.1), (2*TP-0.1, 2*TP-D0), (2*TP+D0, 2*TP+0.1)],
              'T20': [(0.1, TP-0.1), (TP+0.1, 2*TP-0.1), (2*TP+0.1, math.pi)]}
    logf = sys.stdout
    t0 = time.time(); allok = True
    for key in (['T40', 'T20'] if which == 'all' else [which]):
        T = 40 if key == 'T40' else 20
        for (a, b) in pieces[key]:
            ok, worst, nr = run_piece(T, N0, a, b, 0.01 if T == 40 else 0.02, logf)
            print('PIECE T=%d [%.6f,%.6f]: %s rects=%d min margin(N0)=%.2f (%.0fs)' % (T, a, b, 'CERTIFIED' if ok else 'FAILED', nr, worst, time.time()-t0), flush=True)
            allok = allok and ok
    print('N0=%d (F) %s: %s' % (N0, which, 'CERTIFIED for all n >= N0' if allok else 'FAILED'))
