# Acceptance test for B_exact: exact c(m) vs the certificate's own window integral.
#   Certified: Im J in e^{Rh0}/(2pi) [main_lo - E, main_hi + E],  E = rest+eps+crude,  c = -2 Im J + R.
#   Here we check c/(-2 e^{Rh0} sigma sqrt(2pi)/(2pi)) lies within 1 +- (B + window slack) at single m.
import sys, math
sys.argv=['x','244']
src=open('B_exact.py').read(); exec(src[:src.index("if __name__")])
exec(open('exact_crt.py').read().split('if __name__')[0])
ms=[1466,2684,10745,26864,53726,107456]
ps,res=coeffs(244,max(ms))
for m in ms:
    c=crt_at(ps,res,m); mb=arb(m)
    w0=math.sqrt(1.0966/((m-N)/N**2))/N
    B,W,st=certify_ball(mb,w0)
    ser=Rh_series(W,mb,3).coeffs(); sig=1/(2*ser[2]).sqrt()
    S0,QD0,D0=hval(acb(W),float(W.lower()))
    logRh0=a/W+(mb-arb(1)/4)*W+S0.real+QD0.real.log()+(D0.real.sin()/D0.real).log()-N*W
    lg=float((logRh0+(sig*(2*PI).sqrt()/(2*PI)).log()+arb(2).log()).mid())
    lc=(abs(c).bit_length()-60)*math.log(2)+math.log(abs(c)>>(abs(c).bit_length()-60))
    print('m=%6d mu=%.5f  c/(-2 main_gauss) = %.5f   certified B = %.4g'%(m,m/N**2,math.exp(lc-lg)*(1 if c<0 else -1),(float(B.upper()) if B is not None else -1)),flush=True)
