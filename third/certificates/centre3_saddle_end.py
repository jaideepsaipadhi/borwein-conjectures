# Certify mu(L=2.1) = -p'(2.1)/2 < 0.105 (so the centre mu in [0.105,1/2] has saddle L in [0,2.1]; p'' > 0 there by centre3_peak.py).
# p'(L) = int_0^1 t Q'(tL) dt; on each t-piece the integral is (t1^2-t0^2)/2 * Q'(xi), xi in piece (real MVT).
from flint import arb, acb, acb_series, ctx
ctx.prec = 80
Z = [acb.exp_pi_i(acb(2*k)/5) for k in range(5)]
def Qd1(s):
    ctx.cap = 2; x = acb_series([s, 1]); f = acb_series([0])
    for c in range(1, 5): f = f + (1 - Z[c]*(-x).exp()).log()
    return f.coeffs()[1]
NT = 4096; L = arb(2.1); pd = arb(0)
for i in range(NT):
    t0, t1 = arb(i)/NT, arb(i+1)/NT
    pd += (t1**2 - t0**2)/2 * Qd1(acb(arb((i+0.5)/NT, 0.5/NT)*L)).real
print("p'(2.1) in", pd, " mu(2.1) = -p'/2 in", -pd/2, " < 0.105:", (-pd/2) < arb(0.105))
