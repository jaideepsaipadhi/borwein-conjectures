# Band lemma certificate:  g(i) = p3(3i) + p3(3i+1) < 0  for all i >= 0.
#  (a) exact integers for 0 <= i <= I0 (p3 from P_inf^3 = (q;q)^3/(q^3;q^3)^3, Jacobi triple product + 3-colour partitions);
#  (b) i > I0 from the sz_joint expansion (G_3, delta = 3) with rigorous Bessel envelopes, in ball arithmetic:
#      g(i) = -sqrt3 (P3(3i+1) - P3(3i)) + rho(3i) + rho(3i+1),  |rho(j)| <= 2 P6(j) + K(j)^2 P9(j) + 50.1501,
#      P3(j+1) - P3(j) >= a I2(2 sqrt(a(j-1/4)))/(j+3/4)          [d/dx (sqrt(a/x) I1(2 sqrt(ax))) = a I2(2 sqrt(ax))/x]
#      I1(x) <= e^x/sqrt(2 pi x),  I2(x) >= kappa e^x/sqrt(2 pi x) for x >= x0,  kappa = sqrt(2 pi x0) e^{-x0} I2(x0)
#      (sqrt(2 pi x) e^{-x} I_nu(x) is increasing for nu >= 1/2).
#  The resulting bound is sum_t C_t j^{p_t} e^{-g_t sqrt(j - 1/4)} with p_t <= 3, g_t >= sqrt(a); it is decreasing for
#  sqrt(j) >= 2 p_t/g_t, and is < 1 at j = 3 I0.
import sys, math, time
from flint import arb
I0=int(sys.argv[1]) if len(sys.argv)>1 else 36000
t0=time.time(); J=I0+2
sig=[0]*(J+1)
for d in range(1,J+1):
    for mlt in range(d,J+1,d): sig[mlt]+=d
b=[0]*(J+1); b[0]=1
for j in range(1,J+1):
    s=0
    for t in range(1,j+1): s+=sig[t]*b[j-t]
    b[j]=3*s//j
def p3(i):
    s=0; k=0
    while k*(k+1)//2<=i:
        r=i-k*(k+1)//2
        if r%3==0 and r//3<=J: s+=(-1)**k*(2*k+1)*b[r//3]
        k+=1
    return s
# sanity: p3 against direct product for small indices
M=300; ser=[0]*(M+1); ser[0]=1
for k in range(1,M+1):
    if k%3==0: continue
    for _ in range(3):
        for i in range(M,k-1,-1): ser[i]-=ser[i-k]
assert all(p3(i)==ser[i] for i in range(M+1))
bad=[i for i in range(0,I0+1) if not p3(3*i)+p3(3*i+1)<0]
print('(a) exact: g(i) < 0 for all 0 <= i <= %d : %s  (%.0fs)'%(I0, 'PASS' if not bad else 'FAIL %s'%bad[:5], time.time()-t0),flush=True)
# (b)
a=arb.pi()**2/9; PI=arb.pi()
j=arb(3*I0)
def x_k(k,jj): return (2*PI/k)*(jj-arb(1)/4).sqrt()*(arb(1)/2).sqrt()*2**arb(0.5)   # (2pi/k) sqrt(B C), B=1/2, C=2j-1/2
def Pk_upper(k,jj):
    x=x_k(k,jj); return (2*PI/k)*((arb(1)/2)/(2*jj-arb(1)/2)).sqrt()*x.exp()/(2*PI*x).sqrt()
x0=x_k(3,j)
kappa=(2*PI*x0).sqrt()*(-x0).exp()*x0.bessel_i(2)
diff_lo=a*kappa*x0.exp()/(2*PI*x0).sqrt()/(j+arb(3)/4)
K2=((4*PI*(j+1)).sqrt()+2)**2
rest=2*(2*Pk_upper(6,j+1)+K2*Pk_upper(9,j+1)+arb('50.1501'))
ratio=rest/(arb(3).sqrt()*diff_lo)
print('(b) at i = %d: remainder/main <= %s'%(I0,ratio.str(5)))
gmin=float((2*PI/9).lower())      # smallest decay rate: the K^2 P9 term decays like e^{-x3/3}, x3/3 = (2 pi/9) sqrt(j-1/4)
jmono=(2*3/gmin)**2
print('    terms C j^p e^{-g sqrt(j)} with p <= 3, g >= %.3f are decreasing for j >= %.1f < 3 I0 = %d'%(gmin,jmono,3*I0))
ok=(not bad) and ratio.upper()<1 and jmono<3*I0
print('BAND LEMMA CERTIFICATE:', 'PASS' if ok else 'FAIL')
