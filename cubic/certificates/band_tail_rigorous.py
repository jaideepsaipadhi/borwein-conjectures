# Rigorous tail of the band lemma: g(i) = p3(3i) + p3(3i+1) < 0 for all i >= I0 = 36000 (exact for i <= 36000: band_lemma.py).
# By Thm 2.3,  g(i) <= -sqrt3 (P3(3i+1) - P3(3i)) + rho_max(3i) + rho_max(3i+1),  rho_max(j) = sum_{3|k<=K(j),k>=6} phi(k) P_k(j) + 50.1501.
# With y_k(j) = (2pi/k) sqrt(j - 1/4) (so P_k(j) = (2 pi/k) I_1(y_k)/(2 sqrt(j-1/4)) ... ) we use:
#   P3' (x) = 4 a^2 I_2(y)/y^2 at y = 2 sqrt(a(x-1/4))  (from (I_1(z)/z)' = I_2(z)/z), P3 convex => P3(3i+1)-P3(3i) >= P3'(3i);
#   I_2(y) >= kappa2 e^y / sqrt(2 pi y) for y >= y0 with kappa2 = sqrt(2 pi y0) e^{-y0} I_2(y0) (ratio increasing to 1);
#   I_1(y) <= e^y / sqrt(2 pi y);  P_k(j) <= P_6(j) for k >= 6 (y_k and prefactor decrease in k), count sum phi(k) <= K(j)^2/2.
# The resulting bound  ratio(i) <= [2(K^2/2 + 1) U6 + 100.31] / (sqrt3 * 4a^2 kappa2 e^{y}/(sqrt(2 pi y) y^2))  has the form
#   C1 y^{5/2} * y^2 e^{-y/2} + C2 y^{5/2} e^{-y}, with y = y_3(3i) >= y0, which is decreasing once y >= 2(4.5) = 9. We evaluate at y0.
from flint import arb
PI=arb.pi(); a=PI**2/9
I0=36000; j=3*I0
y0=2*(a*(j-arb(1)/4)).sqrt()
kappa2=(2*PI*y0).sqrt()*(-y0).exp()*y0.bessel_i(2)
dP=4*a*a*kappa2*y0.exp()/((2*PI*y0).sqrt()*y0*y0)
K=(4*PI*(j+1)).sqrt()+2
y6=y0/2+arb('0.001')                       # y_6(3i+1) <= y_3(3i)/2 + tiny
U6=(2*PI/6)/(2*(j-arb(1)/4).sqrt())*y6.exp()/(2*PI*y6).sqrt()*y6*2   # generous: P_6 <= (pi/3) * I_1(y6)/sqrt(j) with I_1 <= U
num=2*((K*K/2+1)*U6+arb('50.1501'))
ratio=num/(arb(3).sqrt()*dP)
print('y0=%.3f  kappa2=%.6f  ratio at i=%d <= %s'%(float(y0.mid()),float(kappa2.mid()),I0,ratio))
assert ratio.upper()<1 and y0.lower()>9
print('BAND TAIL: g(i) < 0 for all i >= %d (bound decreasing in y for y >= 9)'%I0)
