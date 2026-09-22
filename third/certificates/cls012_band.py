# Band region, classes 0,1,2, n > 2000 (N = 5n+1 >= 10001), all m < 4N.  c(m) = s(m) + sum_{j>=N} e(j) s(m-j).
# Weights: j in [N,2N): e(j) <= 1 (one part only); j < 4N: e(j) <= 1 + N + N^2/2 <= N^2 <= m^2.
# Hence |c(m) - s(m)| <= A(m-N) + m^2 A(m-2N) <= A(m-Nmin) + m^2 A(m-2Nmin),  A(I) = sum_{0<=i<=I}|s(i)| (0 if I<0),
# Nmin(m) = max(10001, floor(m/4)+1)   (m < 4N  <=>  N >= floor(m/4)+1).
# PART A (exact integers, m <= 200000): sign(s(m)) = sign(a(m)) and |s(m)| > A(m-Nmin) + m^2 A(m-2Nmin).
# PART B (arb, m > 200000): analytic majorant, each term C m^p e^{-kappa sqrt m} with sqrt(m) > 2p/kappa => decreasing.
import pickle, math
from flint import arb, ctx
ctx.prec=200
s=pickle.load(open('g5_coeffs.pkl','rb')); IM=len(s)-1; assert IM>=200000
IM=200000
A=[0]*(IM+1); acc=0
for i in range(IM+1): acc+=abs(s[i]); A[i]=acc
Af=lambda I: A[I] if I>=0 else 0
sgn={0:1,1:-1,2:-1}
zeros=[]; badsign=[]; badband=[]; worst=0.0
for m in range(IM+1):
    r=m%5
    if r>2: 
        assert s[m]==0; continue
    Nmin=max(10001,m//4+1)
    corr=Af(m-Nmin)+m*m*Af(m-2*Nmin)
    if s[m]==0 and corr==0: zeros.append(m); continue   # c(m)=s(m)=0: allowed (non-strict)
    if s[m]*sgn[r]<=0: badsign.append(m); continue
    if corr>0:
        if not abs(s[m])>corr: badband.append(m)
        else: worst=max(worst, corr/abs(s[m]) if corr<abs(s[m]) else 1.0)
# alpha = 1 check on the exact range
alpha_ok=all(abs(s[i])<=math.exp(2*math.sqrt(2*math.pi**2/75*i)) for i in range(0,IM+1,1) if s[i]!=0)
print('PART A: zero coefficients (c=0, allowed):',zeros,' sign failures',len(badsign),badsign[:5],' band failures',len(badband),badband[:5],
      ' worst corr/|s(m)| = %.3e'%worst, ' |s(i)|<=e^{2sqrt(a i)} on [0,2e5]:',alpha_ok)
# ---- PART B ----
PI=arb.pi(); a=2*PI**2/75; c=2*a.sqrt(); m0=arb(200000); ECON=arb('27.4636')
amin=(5-arb(5).sqrt())/2; amax=(5+arb(5).sqrt())/2
Y0=c*(m0-arb(1)/6).sqrt()
# P5(m) >= (2pi/5)(6m-1)^{-1/2} e^Y (1-1/Y)/sqrt(2 pi Y)  (I_1 >= I_{3/2})  >= h0 m^{-3/4} e^{c sqrt m}, m>=m0:
#   (6m-1)^{-1/2} >= (6m)^{-1/2};  Y <= c sqrt m;  Y - c sqrt m >= -c/(12 sqrt(m-1/6)) >= -c/(12 sqrt(m0-1/6));  1-1/Y >= 1-1/Y0
h0=(2*PI/5)/arb(6).sqrt()/(2*PI*c).sqrt()*(-c/(12*(m0-arb(1)/6).sqrt())).exp()*(1-1/Y0)
# P5(m) <= (2pi/5)(6m-1)^{-1/2} e^Y/sqrt(2 pi Y) (I_1 <= I_{1/2}) <= u0 e^{c sqrt m},  u0 at m0 (decreasing factors)
u0=(2*PI/5)/(6*m0-1).sqrt()/(2*PI*Y0).sqrt()
# R(m) <= Phi(K) P10(m) + E_con, Phi(K) <= K^2/10+K/2, K <= sqrt(4 pi m)+2;  P10(m) <= (pi/5)(6m-1)^{-1/2} e^{Y/2}/sqrt(pi Y) <= e^{c sqrt m/2}
#   => R(m) <= 1.3 m e^{c sqrt m/2} + E_con   for m >= m0 (checked below)
K0=(4*PI*m0).sqrt()+2; phi_ratio=(K0*K0/10+K0/2)/m0          # (K^2/10+K/2)/m is decreasing in m
P10fac=(PI/5)/(6*m0-1).sqrt()/(PI*Y0/2).sqrt()
assert (phi_ratio*P10fac).upper()<1.3
# alpha=1 for m >= m0: |s| <= amax*u0 e^{c sqrt m} + 1.3 m e^{c sqrt m/2} + E_con  <= e^{c sqrt m}  (terms decrease)
al=amax*u0+arb('1.3')*m0*(-c*m0.sqrt()/2).exp()+ECON*(-c*m0.sqrt()).exp()
print('PART B: alpha bound for i>=2e5: |s(i)|/e^{2 sqrt(a i)} <=',al.str(5))
# Q(m) = [B(3m/4) + m^2 B(m/2) + R(m)]/(amin P5(m)),  B(I) = e^{c sqrt I}(1+2 sqrt(I)/c)  (sum_{i<=I} e^{c sqrt i})
k1=c*(1-arb(3).sqrt()/2); k2=c*(1-1/arb(2).sqrt()); k3=c/2; k4=c
sm=m0.sqrt()
T1=m0**arb(0.75)*(1+(3*m0).sqrt()/c)*(-k1*sm).exp()
T2=m0**arb(2.75)*(1+(2*m0).sqrt()/c)*(-k2*sm).exp()
T3=arb('1.3')*m0**arb(1.75)*(-k3*sm).exp()
T4=ECON*m0**arb(0.75)*(-k4*sm).exp()
Q=(T1+T2+T3+T4)/(amin*h0)
# monotonicity: monomial m^p e^{-k sqrt m} decreasing for sqrt m > 2p/k; the largest p per term:
for p,k in ((1.25,k1),(3.25,k2),(1.75,k3),(0.75,k4)): assert (sm-2*p/k).lower()>0
print('PART B: sup_{m>=2e5} |c(m)-s(m)| + R  over  |a|min P5 :  Q <=',Q.str(5))
