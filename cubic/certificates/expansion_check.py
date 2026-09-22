# Extends the exact check of the sz_joint expansion bound |p3(i) - sqrt3 chi(i) P3(i)| <= sum_{3|k<=K(i),k>=6} phi(k) P_k(i) + 50.1501
# from i <= 3000 to i <= 3*J using  P_infty^3 = (q;q)^3/(q^3;q^3)^3,  (q;q)^3 = sum_k (-1)^k (2k+1) q^{k(k+1)/2}  (Jacobi),
# and 1/(q;q)^3 coefficients b(j) from  j b(j) = 3 sum_{t=1}^{j} sigma(t) b(j-t).
import sys, math, time, mpmath as mp
J=int(sys.argv[1]); t0=time.time()
sig=[0]*(J+1)
for d in range(1,J+1):
    for mlt in range(d,J+1,d): sig[mlt]+=d
b=[0]*(J+1); b[0]=1
for j in range(1,J+1):
    s=0
    for t in range(1,j+1): s+=sig[t]*b[j-t]
    b[j]=3*s//j
print('b done %.0fs'%(time.time()-t0),flush=True)
M=3*J
def p3(i):
    s=0; k=0
    while k*(k+1)//2<=i:
        r=i-k*(k+1)//2
        if r%3==0: s+=(-1)**k*(2*k+1)*b[r//3]
        k+=1
    return s
mp.mp.dps=int(2*math.sqrt(1.1*3*J)/2.302)+40; a=mp.pi**2/9
def Pk(k,i):
    B=mp.mpf(1)/2; C=2*i-B
    return (2*mp.pi/k)*mp.sqrt(B/C)*mp.besseli(1,(2*mp.pi/k)*mp.sqrt(B*C))
chi=lambda i:(1,-1,0)[i%3]
worst=0; checked=0
for i in list(range(2990,M+1,31)):
    main=mp.sqrt(3)*chi(i)*Pk(3,i); rest=p3(i)-main
    K=int(math.sqrt(4*math.pi*i))+1
    bound=sum(sum(1 for h in range(1,k) if math.gcd(h,k)==1)*Pk(k,i) for k in range(6,K+1,3))+mp.mpf('50.1501')
    r=float(abs(rest)/bound); worst=max(worst,r); checked+=1
print('checked %d values of i in [2990, %d]: max |rest|/bound = %.4f  (%.0fs)'%(checked,M,worst,time.time()-t0),flush=True)
