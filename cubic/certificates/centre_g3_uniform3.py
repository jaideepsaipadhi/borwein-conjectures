# Uniform-in-N bound for the centre inner-window constant G3 = H3(window) sigma^3, valid for N >= N0 = 733.
#   H3(window) <= N^4 [ h3_hi + 3 (TV+S)/N ],  sigma^3 <= (N^3 a2lo - 3N^2 V1)^{-3/2}  =>  G3 <= N^{-1/2} (h3_hi+3VS/N)/(a2lo-3V1/N)^{3/2}
#   window in scaled variables: xi' in xi +- N T, t in +-(0.01 + N T), N T = c N sigma <= c/sqrt(N0 a2lo - 3 V1): shrinks with N.
from flint import arb, acb, ctx
import sys, math
ctx.prec=80
exec(open('centre_mid_uniform.py').read().split("def mid_uniform")[0])
CELLS=100
def g3_uniform(xi,cs=(2.5,3.0,3.5)):
    N=arb(N0)
    lo0,V1,_=rho_data(xi,arb(0,0.02))
    # V1 over zone-1 t range is what mid_uniform used; here only |t|<=0.02+NT is needed, but reuse a safe V over |t|<=0.2
    pass
    a2lo=arb(lo0.lower())
    out={}
    for c in cs:
        NT=arb(c)/(N*a2lo-3*V1).sqrt()
        assert (0.02+NT).upper()<0.9
        xib=xi+arb(0,NT.upper()); tb=arb(0,(0.02+NT).upper())
        taus,h=cellballs(0,1,CELLS)
        SUB=5
        xlo,xhi=float(xib.lower()),float(xib.upper()); tlo,thi=float(tb.lower()),float(tb.upper())
        xs=[arb(xlo+(xhi-xlo)*(i+0.5)/SUB,(xhi-xlo)/(2*SUB)) for i in range(SUB)]
        ts=[arb(tlo+(thi-tlo)*(i+0.5)/SUB,(thi-tlo)/(2*SUB)) for i in range(SUB)]
        h3=arb(0); VS=arb(0)
        for j in range(2):
            tv=arb(0); sup=arb(0)
            for tau in taus:
                qmax=None; dqmax=None
                for xb in xs:
                    for tbb in ts:
                        aa=(-xb*tau).exp(); ps=2*PI*(j+1)/3+tau*tbb
                        d=(1-aa*ps.cos())**2+(aa*ps.sin())**2
                        assert d.lower()>0, ('d touches 0',xi,c,tau)
                        q=tau**3*aa*(1+aa)/(d*d.sqrt())
                        ap=-xb*aa; dp=-2*ap*ps.cos()+2*aa*ps.sin()*tbb+2*aa*ap
                        dq=3*tau*tau*aa*(1+aa)/(d*d.sqrt())+tau**3*(ap*(1+2*aa)/(d*d.sqrt())-arb(1.5)*aa*(1+aa)*dp/(d*d*d.sqrt()))
                        qu=q.upper(); du=abs(dq).upper()
                        qmax=qu if qmax is None or qu>qmax else qmax
                        dqmax=du if dqmax is None or du>dqmax else dqmax
                h3+=arb(qmax)*h; tv+=arb(dqmax)*h
                sup=max(sup,arb(qmax),key=lambda z_:z_.upper())
            VS+=tv+sup
        G3=N**arb(-0.5)*(h3+3*VS/N)/(a2lo-3*V1/N)**arb(1.5)
        out[c]=float(G3.upper())
    return out
if __name__=='__main__':
    lo,hi,K=float(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3])
    for i in range(K):
        a_=lo+(hi-lo)*i/K; b_=lo+(hi-lo)*(i+1)/K
        r=g3_uniform(arb((a_+b_)/2,(b_-a_)/2))
        print('xi in [%.4f,%.4f]: G3U c=2.5 %.4f  c=3 %.4f  c=3.5 %.4f'%(a_,b_,r[2.5],r[3.0],r[3.5]),flush=True)
