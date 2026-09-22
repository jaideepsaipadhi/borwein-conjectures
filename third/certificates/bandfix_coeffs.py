# Independent recomputation of s(i) = [q^i] (q;q)_inf/(q^5;q^5)_inf for i <= 200000 and comparison with
# g5_coeffs.pkl; then the exact one-part-block facts:  g(t) = s(5t)+s(5t+1)+s(5t+2),  G(T) = sum_{t<=T} g(t).
import pickle, time
from flint import fmpz_poly
M = 200000
s = pickle.load(open('g5_coeffs.pkl', 'rb'))
assert len(s) == M + 1
# partitions p(k), k <= M//5, by Euler's recurrence
K = M // 5
p = [0] * (K + 1); p[0] = 1
pent = []
j = 1
while True:
    g1 = j * (3 * j - 1) // 2; g2 = j * (3 * j + 1) // 2
    if g1 > K: break
    pent.append((g1, 1 if j % 2 else -1))
    if g2 <= K: pent.append((g2, 1 if j % 2 else -1))
    j += 1
for k in range(1, K + 1):
    tot = 0
    for g, sg in pent:
        if g > k: break
        tot += sg * p[k - g]
    p[k] = tot
# (q;q)_inf = sum_j (-1)^j q^{j(3j-1)/2}, j in Z
eul = [0] * (M + 1)
j = 0
while True:
    done = True
    for jj in (j, -j) if j else (0,):
        e = jj * (3 * jj - 1) // 2
        if e <= M:
            eul[e] += -1 if jj % 2 else 1; done = False
    if done and j > 0: break
    j += 1
P5 = [0] * (M + 1)
for k in range(K + 1): P5[5 * k] = p[k]
prod = fmpz_poly(eul).mul_low(fmpz_poly(P5), M + 1)
c = [int(x) for x in prod.coeffs()] + [0] * (M + 1 - prod.degree() - 1)
assert c == s, 'pkl disagrees with independent recomputation'
print('g5_coeffs.pkl agrees with independent recomputation for i <= %d' % M)
assert all(s[i] == 0 for i in range(M + 1) if i % 5 in (3, 4))
G = 0; nonneg = []; Gpos = []; Gzero = []
for t in range((M - 2) // 5 + 1):
    g = s[5 * t] + s[5 * t + 1] + s[5 * t + 2]
    if g >= 0: nonneg.append((t, g))
    G += g
    if G > 0: Gpos.append(t)
    if G == 0: Gzero.append(t)
print('t range 0..%d:  g(t) >= 0 only at %s;  G(T) > 0 at %s;  G(T) = 0 at %s' % (t, nonneg, Gpos, Gzero))
assert nonneg == [(1, 1), (5, 0)] and Gpos == [] and Gzero == [1]
print('ONE-PART BLOCK (exact): g(t) < 0 for all t <= %d except t = 1 (g=1), t = 5 (g=0); G(T) <= 0 for all T <= %d' % (t, t))
