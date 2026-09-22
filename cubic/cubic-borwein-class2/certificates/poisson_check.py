# Check of Lemma C: eps(w) = -1 + sum_{k != 0} (-i)^k (e^{a/(w+2 pi i k)} - 1),  Re w > 0, |Im w| <= pi:  |eps| <= 2.
# (a) the proof's two ingredients:  sum |rho| <= (a^2/6) e^{a/pi} <= 0.285, and the blocks-of-four bound on the linear part;
# (b) ball-arithmetic enclosure of |eps| on a grid of boxes covering Re w in (0, 10], |Im w| <= pi (tail k>K bounded).
from flint import acb, arb, ctx
import math
ctx.prec=60
PI=arb.pi(); a=PI**2/9; I=acb(0,1)
rho_bd=a*a/6*(a/PI).exp()
print('sum|rho| <= %.4f'%float(rho_bd.upper()))
K=400
def eps_box(W):
    s=acb(-1)
    for k in range(-K,K+1):
        if k==0: continue
        s+=(-I)**(k%4)*((a/(W+2*PI*I*k)).exp()-1)
    # tail |k|>K: terms (-i)^k (a/(w+2 pi i k)) + rho; linear part grouped in 4-blocks is O(sum 1/k^2) -> <= a*8*pi/(pi^2 K) ; rho <= a^2 e^{a/pi}/(pi^2 k^2)
    tail=arb(8)*a/(PI*K)+a*a*(a/PI).exp()/(PI*PI*K)
    return abs(s)+tail
worst=arb(0)
for i in range(40):
    for j in range(40):
        re=arb(10.0*(i+0.5)/40, 10.0/80); im=arb(-math.pi+2*math.pi*(j+0.5)/40, math.pi/40)
        if float(re.lower())<=0: re=arb((1e-9+10.0/40)/2,(10.0/40-1e-9)/2)
        e=eps_box(acb(re,im)); worst=max(worst,arb(e.upper()),key=lambda q:q.upper())
print('max |eps| over boxes Re w in (0,10], |Im w|<=pi: <= %.3f'%float(worst.upper()))
