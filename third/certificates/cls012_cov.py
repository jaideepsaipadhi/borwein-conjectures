# Coverage for the classes-0,1,2 Laplace certificate: every (n,m) with mu = m/(10n^2) < 0.105 has saddle V >= 2.
# Saddle equation: A(V) := alpha(V)/V^2 = mu' := (m - 1/6)/(25 n^2) < 0.4 mu < 0.042.
# A is continuous on (0,inf), A(V) -> 1/5 (V -> 0), A(V) -> 0 (V -> inf), so a root exists (IVT); it suffices that
# A(V) >= 0.042 on (0, 2]: then EVERY root has V > 2.
#  (i)  V in [0.001, 2]: alpha is increasing (alpha' = V q(V)/5 > 0, q > 0 because y/(e^y-1) decreases), so on
#       a tile [v1,v2], A(V) >= alpha(v1)/v2^2 with alpha(v1) evaluated at the exact point v1.
#  (ii) V in (0, 0.001]: alpha(V) = (1/5) int_0^V s q(s) ds, q = -Q' = 1/(e^s-1) - 5/(e^{5s}-1) >= 2 - 25 s/12
#       (from y/(e^y-1) >= 1 - y/2 and y/(e^y-1) <= 1 - y/2 + y^2/12), so A(V) >= (2 - 25V/12)/10 >= 0.199.
exec(open('cls012_lib.py').read())
ctx.prec=200
v=arb('0.001'); worst=None; cnt=0
while v.lower()<2:
    w=v*arb('0.002') if float(v.mid())<0.05 else arb('0.0005')
    A=alpha(v)/(v+w)**2
    assert A.lower()>0.042, (v,A)
    worst=A if worst is None or A.lower()<worst.lower() else worst
    v=v+w; cnt+=1
print('A(V) >= 0.042 on [0.001, 2]:',cnt,'tiles, min A =',worst.lower())
print('A(V) on (0,0.001] >= ',(2-arb(25)/12*arb('0.001'))/10)
# independent check of the saddle map against direct root finding
import mpmath as mp
for n,m in ((291,1.05*291**2),(1000,1.049*1000**2),(10**6,1.0499e12)):
    mu=(mp.mpf(m)-mp.mpf(1)/6)/(25*mp.mpf(n)**2)
    f=lambda V: float((alpha(arb(V))/arb(V)**2).mid())-float(mu)
    print(n, 'root V ~', mp.findroot(f,2.2))
