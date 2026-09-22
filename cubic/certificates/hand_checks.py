# Machine checks of the hand lemmas (HAND_LEMMAS.md).
import sympy as sp, mpmath as mp, math
# H2 midpoint Peano kernel: g(0) - int_{-1/2}^{1/2} g = -int K g'' with K = (1/2-|t|)^2/2, sup K = 1/8
t=sp.symbols('t'); g=sp.Function('g')
f=sp.exp(3*t)*sp.sin(t)+t**5   # test function
lhs=f.subs(t,0)-sp.integrate(f,(t,-sp.Rational(1,2),sp.Rational(1,2)))
K=(sp.Rational(1,2)-sp.Abs(t))**2/2
rhs=-(sp.integrate(((sp.Rational(1,2)+t)**2/2)*sp.diff(f,t,2),(t,-sp.Rational(1,2),0))+sp.integrate(((sp.Rational(1,2)-t)**2/2)*sp.diff(f,t,2),(t,0,sp.Rational(1,2))))
print('H2 Peano identity residual:', sp.N(lhs-rhs,30), ' sup K =', sp.Rational(1,8))
# H5 bookkeeping: L(m) = sum_j e(j) sqrt3 chi(m-j) P3(m-j)  equals  -2 Im [q^m] P(q) E_n(omega q) ... up to the stated sign convention
mp.mp.dps=50
n=20; N=3*n+1; M=900
e=[0]*(M+1); e[0]=1
for K_ in range(3*n+1,M+1):
    if K_%3==0: continue
    for _ in range(3):
        for i in range(K_,M+1): e[i]+=e[i-K_]
def P3(i):
    if i<1: return mp.mpf(0)
    B=mp.mpf(1)/2; C=2*i-B
    return (2*mp.pi/3)*mp.sqrt(B/C)*mp.besseli(1,(2*mp.pi/3)*mp.sqrt(B*C))
chi=lambda i:(1,-1,0)[i%3]
om=mp.expjpi(mp.mpf(2)/3)
worst=0
for m in (400,600,800,899):
    L1=sum(e[j]*mp.sqrt(3)*chi(m-j)*P3(m-j) for j in range(m+1))
    coef=sum(e[j]*om**j*P3(m-j) for j in range(m+1))        # [q^m] P(q) E_n(omega q)
    for form,val in (('-2 Im',-2*coef.imag),('+2 Im',2*coef.imag)):
        r=abs(val-L1)/abs(L1)
        if form=='-2 Im' or True: print('  m=%d %s: rel diff %s'%(m,form,mp.nstr(r,3)))
# H4/H5 Poisson form: P(e^{-w}) = sum_i P3(i) e^{-w i} = e^{-w/4}(e^{a/w} - 1 + eps(w)), |eps| <= 2
a=mp.pi**2/9; worst=0
for w in (mp.mpc(0.05,0), mp.mpc(0.05,1.0), mp.mpc(0.02,2.5), mp.mpc(0.1,-3.0), mp.mpc(0.03,0.5)):
    S=mp.nsum(lambda i: P3(int(i))*mp.e**(-w*i),[1,mp.inf]) if False else sum(P3(i)*mp.e**(-w*i) for i in range(1,int(60/mp.re(w)**2*0+ 4000/mp.re(w))))
    eps=S/mp.e**(-w/4)-(mp.e**(a/w)-1)
    worst=max(worst,abs(eps)); print('  w=%s |eps| = %s'%(mp.nstr(w,3),mp.nstr(abs(eps),4)))
print('H4 max |eps| at sample points:', mp.nstr(worst,4), '(<= 2 claimed; ball-certified <= 1.42 in poisson_check.py)')
