# Exact band certification: for every n in [n_lo, n_hi] compute c_n(m) = [q^m] (q;q)_{5n}/(q^5;q^5)_n EXACTLY
# (integer arithmetic) for all 0 <= m < 4N, N = 5n+1, and check the Borwein sign pattern in ALL FIVE classes.
# Method: c_n = G_5 * E_n with E_n = prod_{k>5n, 5 not| k} (1-q^k)^{-1}.  Below q^{4N} only partitions into <= 3
# parts >= N contribute, so with p1 = sum_{N<=k<4N, 5 not| k} q^k (cycle index of S_2, S_3):
#   6 E_n = 6 + 6 p1 + 3 (p1^2 + p1(q^2)) + (p1^3 + 3 p1(q^2) p1 + 2 p1(q^3))   (mod q^{4N}).
# Also records, for classes 3,4 with 2N <= m < 4N, the exact ratio (b2+b3)/|b1| where b1 = G(T) (one-part block).
import pickle, sys, time
from multiprocessing import Pool
from flint import fmpz_poly
s = pickle.load(open('g5_coeffs.pkl', 'rb'))

def E6(n):
    N = 5 * n + 1; L = 4 * N
    c1 = [0] * L
    for k in range(N, L):
        if k % 5: c1[k] = 1
    def sub(d):
        out = [0] * L
        for k in range(N, L):
            if k % 5 and k * d < L: out[k * d] = 1
        return fmpz_poly(out)
    p1 = fmpz_poly(c1); p1_2 = sub(2); p1_3 = sub(3)
    sq = p1.mul_low(p1, L)
    return 6 + 6 * p1 + 3 * (sq + p1_2) + sq.mul_low(p1, L) + 3 * p1_2.mul_low(p1, L) + 2 * p1_3, L

def coeffs(n):
    E, L = E6(n)
    c6 = E.mul_low(fmpz_poly(s[:L]), L)
    cf = [int(x) for x in c6.coeffs()]; cf += [0] * (L - len(cf))
    assert all(x % 6 == 0 for x in cf)
    return [x // 6 for x in cf]

def validate(nmax=14):
    for n in range(1, nmax + 1):
        N = 5 * n + 1; L = 4 * N
        P = fmpz_poly([1])
        for k in range(1, 5 * n + 1):
            if k % 5: P = P.mul_low(fmpz_poly([1] + [0] * (k - 1) + [-1]), L)
        d = [int(x) for x in P.coeffs()]; d += [0] * (L - len(d))
        assert d == coeffs(n), n
    print('construction validated against the direct product S_n for n = 1..%d' % nmax)

# prefix sums of g for b1
Gp = []; acc = 0
for t in range((len(s) - 3) // 5 + 1):
    acc += s[5 * t] + s[5 * t + 1] + s[5 * t + 2]; Gp.append(acc)

def check(n):
    N = 5 * n + 1; c = coeffs(n)
    bad = [m for m in range(4 * N) if (c[m] < 0 if m % 5 == 0 else c[m] > 0)]
    worst = (0.0, None); zeros34 = 0
    for m in range(4 * N):
        r = m % 5
        if r in (3, 4):
            T = (m - r) // 5 - n
            b1 = Gp[T] if T >= 0 else 0
            if m < 2 * N:
                assert c[m] == b1          # one-part block only
                if c[m] == 0: zeros34 += 1
            else:
                b23 = c[m] - b1
                rat = b23 / -b1           # b1 < 0 here (T >= n >= 291)
                if abs(rat) > abs(worst[0]): worst = (rat, m)
    return n, bad, worst, zeros34

if __name__ == '__main__':
    validate()
    lo, hi = int(sys.argv[1]), int(sys.argv[2])
    t0 = time.time(); nbad = 0; W = (0.0, None, None)
    with Pool(2) as P:
        for n, bad, (rat, m), z in P.imap_unordered(check, range(lo, hi + 1), chunksize=4):
            if bad: nbad += 1; print('FAIL n=%d m=%s' % (n, bad[:10]), flush=True)
            if abs(rat) > abs(W[0]): W = (rat, n, m)
            if n % 250 == 0: print('  n=%d done (%.0fs)' % (n, time.time() - t0), flush=True)
    print('n in [%d,%d]: sign failures (any class, m < 4N): %d;  max |b2+b3|/|b1| (classes 3,4, 2N<=m<4N) = %.4e at n=%s, m=%s (t=m/N=%.3f)'
          % (lo, hi, nbad, W[0], W[1], W[2], W[2] / (5 * W[1] + 1)))
