# Adaptive driver for B_exact: bisect failing m-balls down to single m.  Usage: python3 run_exact.py n mu_a mu_b K0
import sys, math
args=sys.argv[:]
n=int(args[1]); sys.argv=['x',str(n)]
src=open('B_exact.py').read(); exec(src[:src.index("if __name__")])
mu_a,mu_b,K0=float(args[2]),float(args[3]),int(args[4])
edges=[mu_a*(mu_b/mu_a)**(j/K0) for j in range(K0+1)]
stack=[]
for j in range(K0):
    lo=max(math.ceil(edges[j]*N*N), 2*N if j==0 else 0); hi=math.floor(edges[j+1]*N*N)
    if hi>=lo: stack.append((lo,hi))
worst=0.0; nballs=0; hard=[]
while stack:
    lo,hi=stack.pop()
    mb=arb((lo+hi)/2,(hi-lo)/2)
    w0=math.sqrt(1.0966/max((float(mb.mid())-N)/N**2,1e-9))/N
    try:
        B,W,st=certify_ball(mb,w0)
    except Exception as e:
        B,W,st=None,None,'exc:'+str(e)[:40]
    ok = B is not None and B.upper()<1
    if ok:
        worst=max(worst,float(B.upper())); nballs+=1
    elif hi>lo:
        mid=(lo+hi)//2; stack.append((lo,mid)); stack.append((mid+1,hi))
    else:
        hard.append((lo,st if B is None else float(B.upper())))
print('EXACT-LAPLACE n=%d m in [%d,%d]: certified balls %d, UNRESOLVED single m %d %s ; worst B = %.4f'%(n,math.ceil(mu_a*N*N),math.floor(mu_b*N*N),nballs,len(hard),hard[:5],worst))
