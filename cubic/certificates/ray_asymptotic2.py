# Closing lemma, v >= V0 = 3000, all s = eps v^2/a in (0, S1], with the perturbation constants ENCLOSED (supersedes ray_asymptotic.py).
# Disc |zeta - v| <= v/2 (Re zeta >= v/2, |zeta| <= 3v/2).  x = r zeta eps.
#  (i)  |St(zeta)| <= MS := sum_r (3/r) e^{-r v/2} 2 (2/(3 r v) + eps_max)             [termwise magnitude bound; decreasing in v]
#  (ii) |h(x)| <= 4 for every r on the disc:  |x| < 0.9 -> Bernoulli enclosure (<= 1);  |x| >= 0.9 -> |h| <= 2(|x|/(3 Re x)+1) <= 2(1+1)
#       (|x|/Re x = |zeta|/Re zeta <= 3).  So |Qt - c1 h(zeta eps)| <= EQ := sum_{r>=2} c_r e^{-(r-1) v/2} 4,  |Delta| <= DM := (c1*4+EQ) e^{-v/2}.
#  (iii) r=1 term: log h(x) on the x-disc |x - x0| <= X, X = (3/2) a S1 / v: ball enclosure; Mb := diam(log h ball) + 2 EQ/Qlow + DM^2.
#  Derivative perturbations (Cauchy, relative to a k!/v^{k+1}):  d_k <= (MS + eps_max Mb) 2^k v / a.
#  Zone-1 quadratic terms: N Re[St diff] <= N t^2 (8 MS / v^2),  Re[b diff] <= t^2 (4 Mb / v^2)  =>  kappa loss <= 8 MS v/a + 4 Mb eps_max v/a.
#  Zone 2 (u >= 1/2) on the contour: log|Qt| <= log(2(1+3as/v)/(3as/v) * c1 + EQ), |sin Delta| <= |Delta| e^{|Delta|}, |Delta| <= Qtmax e^{-v},
#       N |St| <= (v^2/(a s)) 6 e^{-v}(1/(3v)+eps) :  exponent <= -(v/(5s))(1-d2) + log(4v/(3as)+6) + (12 + 8) v e^{-v}/(a s).
#  s-monotonicity (analytic): d/ds of the zone-2 exponent (with the N^{3/2} prefactor, 1.5 log(1/s)) is
#       v(1-d2)/(5 s^2) - 20 v e^{-v}/(a s^2) - 1/s - 1.5/s  >  0  for all s <= S1  since v/5 > 2.5 S1 + 20 v e^{-v}/a.   Same for R.
#  v-monotonicity: every enclosed constant is decreasing in v (exponential factors), the leading terms improve with v.
from flint import arb, acb
import math
exec(open('uniform_lib.py').read().split("def G_acb")[0])
V0=arb(3000); S1=arb('1.0011'); v=V0
eps_max=a*S1/(v*v)
MS=sum((arb(3)/r*(-r*v/2).exp()*2*(arb(2)/(3*r*v)+eps_max) for r in range(1,60)),arb(0))+(-60*v/2).exp()
c1=3*S3/2
EQ=sum((3*S3/(2*r)*(-(r-1)*v/2).exp()*4 for r in range(2,60) if r%3),arb(0))+(-59*v/2).exp()
DM=(c1*4+EQ)*(-v/2).exp()
X=arb(1.5)*a*S1/v
xball=acb(arb(0,X.upper()),arb(0,X.upper()))
hb=g_acb(3*xball)/(3*g_acb(xball)); lhb=hb.log()
Qlow=c1*arb(hb.real.lower())-EQ
Mb=2*abs(lhb.real).rad()+2*abs(lhb.imag).rad()+2*EQ/Qlow+DM*DM
Mb=arb(Mb.upper())
def dk(k): return (MS+eps_max*Mb)*arb(2)**k*v/a
kap_loss=8*MS*v/a+4*Mb*eps_max*v/a
print('MS <= %s  EQ <= %s  DM <= %s  Mb <= %.3e'%(MS,EQ,DM,float(Mb.upper())))
print('d2..d5 <= ', [float(dk(k).upper()) for k in (2,3,4,5)], ' kappa loss <= %.2e'%float(kap_loss.upper()))
G3n=(3/arb(2).sqrt())*(S1/V0).sqrt()*(1+dk(3))/(1-dk(2))**arb(1.5)
ph4_lo=-6*S1/V0*dk(4); ph4_hi=6*S1/V0*(1+dk(4))/(1-dk(2))**2
G5n=(arb(120)/arb(2)**arb(2.5))*(S1/V0)**arb(1.5)*(1+dk(5))*8/(1-dk(2))**arb(2.5)
c=3.0; cells=600; Lint=arb(0)
for i in range(cells):
    t1=arb(c*i/cells); t2=arb(c*(i+1)/cells)
    th3=G3n*t2**3/6+G5n*t2**5/120; cosl=1-th3*th3/2
    ex_lo=ph4_lo*t2**4/24-G5n*t2**5/120
    Lint+=arb(((-(t2*t2)/2).exp()*cosl*ex_lo.exp()).lower())*(t2-t1)
Lint=2*Lint
kappa=arb('0.8')*(1-dk(2))/(1+dk(2))-kap_loss
Z1=2*(PI/(2*kappa)).sqrt()*(arb(c)*(kappa/2).sqrt()).erfc()*(Mb+DM).exp()
s=S1; d2=dk(2)
Z2=4*PI*arb(2).sqrt()*v**arb(1.5)/(a*s**arb(1.5))*(-(v/(5*s))*(1-d2)+(4*v/(3*s*a)+6).log()+20*v*(-v).exp()/(a*s)).exp()
Rb=5*arb(2).sqrt()*PI*v**arb(1.5)/(a*s**arb(1.5))*(-(v/s)*(1-d2)+2*(v*v/(a*s)).log()+(4*PI*arb('0.2')).log()+arb(100).log()+20*v*(-v).exp()/(a*s)).exp()
weps=6*arb(c)*(ph4_hi*arb(c)**4/24+G5n*arb(c)**5/120).exp()*(-(v/s)*arb('0.9')).exp()
B=(Z1+Z2+Rb+weps)/Lint
mono=(v*(1-d2)/5-20*v*(-v).exp()/a-arb('2.5')*S1)
assert mono.lower()>0
print('Lint >= %.4f  Z1 <= %.4f  Z2 <= %s  R <= %s  =>  B <= %.4f   (s-monotonicity margin %.1f > 0)'%(
 float(Lint.lower()),float(Z1.upper()),Z2,Rb,float(B.upper()),float(mono.lower())))
assert B.upper()<1
print('CLOSING LEMMA (v >= 3000, all N): CERTIFIED')
