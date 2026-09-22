# Randomised checks of H1 (Cauchy derivative estimate) and H3 (Koksma, one node per cell).
import mpmath as mp, random
mp.mp.dps=20; random.seed(1)
w1=0
for trial in range(40):
    c0=mp.mpc(random.uniform(-1,1),random.uniform(-1,1)); rho=mp.mpf(random.uniform(0.3,1)); d=rho*random.uniform(0,0.8)
    p=c0+d*mp.e**(1j*random.uniform(0,6.28)); k=random.randint(1,5)
    f=lambda z: mp.e**(z*z)/(z-3)+mp.log(z+4)
    M=max(abs(f(c0+rho*mp.e**(1j*t))-f(c0)) for t in mp.linspace(0,2*mp.pi,400))*1.02
    w1=max(w1,abs(mp.diff(f,p,k))/(mp.factorial(k)*M*rho/(rho-d)**(k+1)))
print('H1 Cauchy: max lhs/rhs = %s over 40 random cases (<= 1 required)'%mp.nstr(w1,4),flush=True)
w3=0
for trial in range(12):
    Mn=random.randint(5,200); h=mp.mpf(1)/Mn; A=random.uniform(1,30); B=random.uniform(0,5)
    f=lambda x: mp.sin(A*x)*mp.e**(-B*x)+x**2
    fp=lambda x: A*mp.cos(A*x)*mp.e**(-B*x)-B*mp.sin(A*x)*mp.e**(-B*x)+2*x
    xs=[(i+random.random())*h for i in range(Mn)]
    lhs=abs(h*sum(f(x) for x in xs)-mp.quad(f,[0,1]))
    TV=mp.quad(lambda x: abs(fp(x)),mp.linspace(0,1,60))
    w3=max(w3,lhs/(h*TV))
print('H3 Koksma: max lhs/(h TV) = %s over 12 random cases (<= 1 required)'%mp.nstr(w3,4),flush=True)
