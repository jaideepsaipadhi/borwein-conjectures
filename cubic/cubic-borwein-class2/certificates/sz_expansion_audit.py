# Audit of sz_joint Prop. (explicit expansion) specialised to G_3, delta = 3 (the input to R(m) in the cubic proof).
import math, mpmath as mp
mp.mp.dps=30
p=3; delta=3
B=mp.mpf(delta*(p-1))/12
print('1. growth exponent at p|k: pi*delta*(d^2-p)/(12 p z), d=3 ->', mp.mpf(delta*(9-p))/(12*p), '= B =', B, ' (a = pi^2 B^2 *4/... gives a_3 = pi^2/9 OK)')
print('2. correction term grows iff delta > 24/(p-1) =', 24/(p-1), '-> not extracted for delta=3')
# 3. counting: sum over growing cusps (3|k <= N, h coprime) of 1/k  <=  (p-1)/p^2 (N+1)
worst=0
phi=list(range(200001))
for i in range(2,200001):
    if phi[i]==i:
        for j in range(i,200001,i): phi[j]-=phi[j]//i
s=0.0
for N in range(1,200001):
    if N%3==0: s+=phi[N]/N
    if N%1000==0 or N<50: worst=max(worst, s/(N+1))
print('3. max_{N<=2e5} sum_{3|k<=N} phi(k)/k / (N+1) = %.5f  <=  (p-1)/p^2 = %.5f'%(worst,(p-1)/p**2))
print('   (proof: phi(3j)/(3j) <= 2/3 and there are N/3 multiples of 3)')
# 4. exp bound on the paths: exponent pi C z/k^2, C = 2n - 1/2, Re z <= 2k^2/(N+1)^2, N >= sqrt(4 pi n)+1  ->  <= e
print('4. on paths |e^{pi C z/k^2}| <= e^{2 pi C/(N+1)^2} <= e^{4 pi n /(4 pi n)} = e   (N+1 >= sqrt(4 pi n))')
# 5. remainder majorant monotone: terms e^{pi B w} e^{-2 pi j w}, j>=1, need B <= 2j
print('5. remainder terms non-increasing in w: B = %s <= 2'%B)
f=lambda x: 1/mp.qp(x)
t=mp.e**(-2*mp.pi)
Earc=mp.sqrt(2)*mp.pi*mp.mpf(p-1)/p**2*mp.e**(1+mp.pi*B)
Erem=2*mp.sqrt(2)*mp.mpf(p-1)/p**2*mp.e**(1+mp.pi*B)*(f(t)**delta*f(t**p)**delta-1)
Kp=mp.sqrt(p)*mp.e**(-mp.pi*(p-1)/(12*p))*f(mp.e**(-2*mp.pi/p))*f(t)
Eng=2*mp.sqrt(2)*mp.e*Kp**delta
print('6. E_arc=%s E_rem=%s K_3=%s E_ng=%s  E_con=%s'%(mp.nstr(Earc,8),mp.nstr(Erem,6),mp.nstr(Kp,8),mp.nstr(Eng,8),mp.nstr(Earc+Erem+Eng,8)))
# 7. non-growing chord sum: sum_{k<=N} phi(k) (1/k^2)(2 sqrt2 k/(N+1)) <= 2 sqrt2 N/(N+1) <= 2 sqrt2
print('7. non-growing chord count factor <= 2 sqrt2 (sum phi(k)/k <= N)')
# 8. K_3 bound: sup over Re(1/z) >= 1 of |G_3^delta| at p not| k:  sqrt3^delta * e^{-pi delta (p-1) w/(12p)} f(e^{-2 pi w/p})^delta f(e^{-2 pi w})^delta, max at w = 1
for w in (1,1.5,3):
    val=mp.sqrt(p)*mp.e**(-mp.pi*(p-1)*w/(12*p))*f(mp.e**(-2*mp.pi*w/p))*f(mp.e**(-2*mp.pi*w))
    print('   w=%s: bound factor %s'%(w,mp.nstr(val,8)))
