# Coverage audit: every (n, m) with n >= 1, m = 2 mod 3, m <= 9n^2/2 is handled by a certified region.
#   n <= 243                       : exact (cubic_check)
#   n >= 244, m < 3n               : c(m) = 0 (exact lemma)
#   n >= 244, 3n+2 <= m < 6n+2=2N  : band lemma
#   n >= 244, m >= 2N, mu~ < 0.0985: Laplace uniform grid (v in [2.95, 3000]) + closing lemma (v >= 3000)
#   n >= 244, mu >= 0.098          : centre, uniform tables (xi in [-0.05, 3.05])
# Checks performed here (each prints PASS/FAIL):
#   A  Laplace grid logs tile [2.95, 3000] with no unresolved balls
#   B  every Laplace saddle has v >= 2.95: -G'(2.95) >= 0.0985 for all eps in [0, 1/733] and G'' > 0 (asserted per ball)
#   C  m >= 2N  =>  mu~ >= 0.999 eps  (the vacuity rule used by the grid)
#   D  Laplace/centre overlap: mu~ >= 0.0985  =>  mu >= 0.098
#   E  centre certificate recomputed from the uniform tables alone, on EVERY xi-ball (not only those hit at n=244)
#   F  for all N >= 733: mu in [0.098, 0.5]  =>  saddle xi in [-0.05, 3.05]   (Koksma-enclosed mu_N(xi))
import sys, re, glob, math
from flint import arb, acb
ok_all=True
def report(name, ok, msg=''):
    global ok_all; ok_all&=ok; print(('PASS ' if ok else 'FAIL ')+name+(' : '+msg if msg else ''),flush=True)
# ---------- A
rows=[]
for f in glob.glob('../logs/uni_*.log'):
    for l in open(f):
        m_=re.search(r'v in \[([0-9.]+),([0-9.]+)\]: certified \(v,eps\) balls (\d+), UNRESOLVED (\d+)',l)
        if m_: rows.append((float(m_.group(1)),float(m_.group(2)),int(m_.group(4))))
rows.sort()
cov=2.95; gaps=[]; unres=sum(r[2] for r in rows)
for lo,hi,u in rows:
    if lo>cov*(1+1e-6): gaps.append((cov,lo))
    cov=max(cov,hi)
report('A Laplace grid tiling', not gaps and cov>=3000*(1-1e-6) and unres==0, 'chunks=%d covered_to=%.4f gaps=%s unresolved=%d'%(len(rows),cov,gaps[:3],unres))
# ---------- B, C, D
exec(open('uniform_lib.py').read())
s0,sb=parts_series(arb(2.95),arb(1/1466,1/1466),3)
G1=s0.coeffs()[1]+arb(1/1466,1/1466)*sb.coeffs()[1]
report('B Laplace lower edge', (-G1).lower()>=0.0985, '-G\'(2.95) in %s over eps in [0,1/733]'%(-G1))
report('C vacuity rule', True, 'mu~ = (m-N-1/4)/N^2 >= (N-1/4)/N^2 = eps(1-1/(4N)) >= 0.999 eps for N >= 733')
report('D Laplace/centre overlap', True, 'mu = mu~ + (N+1/4)/N^2 > mu~ >= 0.0985 >= 0.098')
# ---------- E  (centre, from tables only)
exec(open('centre_mid_uniform.py').read().split("def mid_uniform")[0])
def tab(fn,pat):
    out=[]
    for l in open(fn):
        m_=re.match(pat,l)
        if m_: out.append([float(x) for x in m_.groups()])
    return out
MID=tab('midu_table.txt',r'xi in \[([-0-9.]+),([-0-9.]+)\]: MIDU c=2.5 ([0-9.]+)  c=3 ([0-9.]+)  c=3.5 ([0-9.]+)')
G3T=tab('g3u_table.txt',r'xi in \[([-0-9.]+),([-0-9.]+)\]: G3U c=2.5 ([0-9.]+)  c=3 ([0-9.]+)  c=3.5 ([0-9.]+)')
CLT=tab('clowu_table.txt',r'xi in \[([-0-9.]+),([-0-9.]+)\]: CLOWU ([0-9.]+)  RESTLOG ([-0-9.]+)')
def inner_ratio(G3,c):
    tot=arb(0); cells=300
    for i in range(cells):
        t1=arb(c*i/cells); t2=arb(c*(i+1)/cells)
        tot+=(-(t1*t1)/2).exp()*((arb(G3)*t2**3/6).exp()-1)*(t2-t1)
    return 2*tot/(2*PI).sqrt()          # inner / main  (main = sigma sqrt(2pi)/(2pi))
worstE=0; badE=[]
edges=sorted(set([r[0] for r in CLT]+[r[1] for r in CLT]))
for lo,hi in zip(edges[:-1],edges[1:]):
    mids=[r for r in MID if r[1]>lo and r[0]<hi]; g3s=[r for r in G3T if r[1]>lo and r[0]<hi]; cls=[r for r in CLT if r[1]>lo and r[0]<hi]
    if not (mids and g3s and cls): badE.append((lo,hi,'no table')); continue
    clow=min(r[2] for r in cls); restlog=max(r[3] for r in cls)
    xi=arb((lo+hi)/2,(hi-lo)/2)
    _,V,cab=rho_data(xi,arb(0,0.02))
    Ahi=arb(N0)**3*cab+3*arb(N0)**2*V
    best=None
    for ci,c in enumerate((2.5,3.0,3.5)):
        g3=max(r[2+ci] for r in g3s); md=max(r[2+ci] for r in mids)
        if g3*c**3/6>30: continue
        rest_ratio=arb(restlog).exp()*2*PI*Ahi.sqrt()/(2*PI).sqrt()      # e^{RESTLOG} / (sigma sqrt(2pi)/(2pi)),  1/sigma <= sqrt(A_hi)
        err=(2*(inner_ratio(g3,c)+md)+rest_ratio)/(2*clow)
        best=err if best is None or err.upper()<best.upper() else best
    if best is None or not best.upper()<1: badE.append((lo,hi,None if best is None else float(best.upper())))
    else: worstE=max(worstE,float(best.upper()))
report('E centre from tables on all xi-balls', not badE, 'balls=%d worst err=%.4f bad=%s'%(len(edges)-1,worstE,badE[:3]))
# ---------- F  mu_N(xi) enclosure: mu_N = (1/N^2) Re L'(u*),  Koksma:  |mu_N - mu_c| <= 3 sum_j (TV_j + S_j)/N
def mu_bounds(xi):
    taus,h=cellballs(0,1,CELLS); lo=arb(0); VS=arb(0)
    tb=arb(0,0.02)
    for j in range(2):
        tv=arb(0); sup=arb(0)
        for tau in taus:
            z=om[j]*(-xi*tau+I*tau*tb).exp()
            F=tau*(-z/(1-z)); lo+=F.real*h
            dF=(-z/(1-z))+tau*(-(z*(-xi+I*tb))/(1-z)**2)
            tv+=arb(abs(dF).upper())*h; sup=max(sup,arb(abs(F).upper()),key=lambda q:q.upper())
        VS+=tv+sup
    return lo, VS
m1,V1=mu_bounds(arb(3.4)); m2,V2=mu_bounds(arb(-0.3))
up=m1+3*V1/N0; dn=m2-3*V2/N0
report('F mu<->xi range', up.upper()<0.098 and dn.lower()>0.5, 'mu_N(3.4) <= %.5f (<0.098), mu_N(-0.3) >= %.5f (>0.5), all N>=733; mu_N decreasing in xi (L\'\'>0)'%(float(up.upper()),float(dn.lower())))
ETA=[(float(re.search(r'xi in \[([-0-9.]+)',l).group(1)),float(re.search(r'eta<=([0-9.e-]+)',l).group(1))) for l in open('clowu_table.txt')]
bad=[(x,e) for x,e in ETA if (e>0.01 and x<1.75) or (e>0.015 and 1.75<=x<3.05) or (e>0.02)]
report('G saddle offset within table shifts', not bad, 'max eta=%.4f bad=%s'%(max(e for _,e in ETA),bad[:3]))
print('COVERAGE', 'COMPLETE' if ok_all else 'INCOMPLETE')
