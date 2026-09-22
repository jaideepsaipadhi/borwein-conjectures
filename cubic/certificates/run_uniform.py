# Adaptive driver for B_uniform over v in [v_a, v_b], all eps in [0, 1/733].  Bisects v-balls and eps-balls on failure.
# Usage: python3 run_uniform.py v_a v_b  [max_depth]      (independent chunks in v can be run in parallel)
import sys, math, time
args=sys.argv[:]
exec(open('B_uniform.py').read().split("if __name__")[0])
va,vb=float(args[1]),float(args[2]); maxdepth=int(args[3]) if len(args)>3 else 14
def eps_cap(vlo):
    # pairs with eps above this are vacuous: m >= 2N gives eps <= mu~/0.999 and mu~ = -G'(v) <= a/v^2 + |St'(v)| + eps|b'(v)|,
    # |St'| <= 30 e^{-v}/v (termwise), eps|b'| <= 1e-3 a/v^2 for v >= 2.9;  margin 1.05.
    return min(1.0/733, 1.05*(1.0966/vlo**2+30*math.exp(-vlo)/vlo)/0.999)
def eps_grid(vhi, vlo=None):
    cap=eps_cap(vlo if vlo else vhi)
    ec=min(1.0/733, 1.0966/(3*vhi*vhi))
    if ec<cap:
        g=[0.0]+[ec*(cap/ec)**(i/3) for i in range(4)]
    else:
        g=[0.0, cap/8, cap/4, cap/2, cap]
    return g
rel0=3e-3 if va<4 else (4e-3 if va<20 else 1e-2)
stack=[]; v=va
while v<vb:
    h=min(vb, v*(1+rel0)); stack.append((v,h,0)); v=h
done=0; worst=0.0; unresolved=[]; t0=time.time()
while stack:
    lo,hi,dep=stack.pop()
    V=arb((lo+hi)/2,(hi-lo)/2)
    eg=eps_grid(hi,lo); epairs=[(eg[i],eg[i+1],0) for i in range(len(eg)-1)]
    # region above the cap: must be vacuous (m >= 2N forces eps <= mu~/0.999); verified rigorously, bisecting eps if needed
    capv=eg[-1]; vac=[(capv,1.0/733,0)] if capv<(1.0/733)*(1-1e-9) else []
    while vac:
        e0,e1,ed=vac.pop()
        Eb=arb((e0+e1)/2,(e1-e0)/2)
        if 1.0966/hi**2 > 0.999*e1:      # a/v^2 <= mu~ + |h'| bound cannot be below eps here: certify directly
            epairs.append((e0,e1,4)); continue
        try: ok=mu_upper(V,Eb).upper()<0.999*e0
        except Exception: ok=False
        if ok: continue
        if ed<12: vac+=[(e0,(e0+e1)/2,ed+1),((e0+e1)/2,e1,ed+1)]
        else: epairs.append((e0,e1,4))
    while epairs:
        e0,e1,ed=epairs.pop()
        E=arb((e0+e1)/2,(e1-e0)/2)
        try: B,st=certify(V,E)
        except Exception as ex: B,st=None,'exc:'+str(ex)[:40]
        if B is not None and B.upper()<1:
            worst=max(worst,float(B.upper())); done+=1; continue
        if ed<4 and e1>0:
            m_=(e0+e1)/2 if e0>0 else e1/8
            epairs+= [(e0,m_,ed+1),(m_,e1,ed+1)]
        elif dep<maxdepth:
            stack+= [(lo,(lo+hi)/2,dep+1),((lo+hi)/2,hi,dep+1)]; break
        else:
            unresolved.append((lo,hi,e0,e1,st if B is None else float(B.upper()))); break
print('UNIFORM-LAPLACE v in [%.4f,%.4f]: certified (v,eps) balls %d, UNRESOLVED %d %s ; worst B = %.4f ; %.0fs'%(va,vb,done,len(unresolved),unresolved[:3],worst,time.time()-t0),flush=True)
