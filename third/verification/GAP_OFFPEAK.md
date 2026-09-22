# GAP_OFFPEAK: the off-peak lemma (B) of the centre regime

*19 Sept 2026. Status: **(B) is certified for every n >= 821.** With piece (A) (GAP_CENTRE2, n >= 821) and the exact
computation (n <= 820), this closes the centre regime mu in [0.105, 1/2] for all n >= 291, subject to the caveat
on piece (F) in section 7.*

## 0. Statement (as needed by GAP_CENTRE2 section 3)

For n >= 821, L in [0, 2.1] and every theta with |5n(theta - 2 pi h/5)| > 1.1 for h = 1..4:

    log|S_n(e^{-L/(5n)} e^{i theta})| - n p(L) <= -(1.5 log n + 6.3),     p(L) = (1/L) int_0^L Q.      (B)

Notation: q = zeta^h e^{-sigma/(5n)}, sigma = L - i y, y = 5n(theta - 2 pi h/5), where h in {0,..,4} is the nearest
5th root; Phi_c(s) = log(1 - zeta^{hc} e^{-s}); omega = sigma/n. E(n) := 1.5 log n + 6.3 (E(821) = 16.36).
By |S(conj q)| = |S(q)|, it suffices to take theta in [0, pi]. Near h = 1, 2 both signs of y are used, which also
covers h = 3, 4.

## 1. Result

| piece | theta-range | L | method | status | worst certified value at N0 = 821 (the bound, minus the requirement) |
|---|---|---|---|---|---|
| (A) | \|y\| <= 1.1 | [0, 2.1] | GAP_CENTRE2 | certified earlier | not part of (B) |
| (E) | h = 1..4, 1.1 <= \|y\| <= 6 | [0.1, 2.1] (+ lemma M) | finite EM, sigma-boxes | **certified**, n >= 821 | V <= -12.0 |
| (E0) | h = 0, \|y\| <= 20 | [0.1, 2.1] (+ lemma M) | finite EM from t = 1, sigma-boxes | **certified**, n >= 821 | V <= -0.36 |
| (N) | h = 0..4, 6 (20 for h = 0) <= \|y\| <= 0.25 n | [0.1, 2.1] (+ lemma M) | analytic triangle bound on the j-series | **certified**, n >= 821 | V <= -1.8 |
| (F) | dist(theta, {0, 2pi/5, 4pi/5}) >= 0.0499 | (0, 2.1] | chunks + maximum principle | **certified**, n >= 821 (float64 interval, see section 7) | margin >= 47.9 |

The pieces overlap: (N) ends at |theta - 2 pi h/5| = 0.05, while (F) starts at 0.0499. Every certificate is uniform
in n >= N0 = 821, so **n_0 = 821**, and no exact computation beyond n = 820 is needed for the centre. The margins are
listed as V <= (negative number), meaning the certified bound for log|S_n| - n p(L) + E(N0) is below 0 by at least
that much. For (E0) and (N), V is taken after adaptive bisection, which stops at the first passing box, so it
understates the true slack: at the same points, floating point gives slack in the hundreds.

Commands and runtimes (2-core box, shared):
- `python3 offpeak_em.py 821 6.0 32` produces `offpeak_em_821.log` in 85 s.
- `python3 offpeak_em0.py 821 20 32` produces `offpeak_em0_821.log` in 8 s.
- `python3 offpeak_tb.py 821 0.25` produces `offpeak_tb_821.log` in 2 s.
- `python3 offpeak_chunk.py 821 all` produces `offpeak_chunk_821.log` in 470 s: 163 rectangles, 5 T=40 pieces and 3 T=20 pieces.
- `python3 offpeak_crosscheck.py` produces `offpeak_crosscheck.log`. It is a sanity check only and not part of the certificate.

Shared code is in `offpeak_common.py`: closed forms, box constructors and rigorous lower bounds for |z| over a ball.

## 2. Two elementary lemmas used throughout

**Lemma M (monotonicity in L). [proved]** For 0 <= L <= L', log|S_n(r e^{i theta})| <= log|S_n(r' e^{i theta})| + n(L' - L),
where r = e^{-L/5n} and r' = e^{-L'/5n}.

Proof:
- |1 - rho e^{ia}|^2 = (1-rho)^2 + 4 rho sin^2(a/2).
- For rho' <= rho <= 1, both (1-rho)^2 <= (1-rho')^2 and rho = (rho/rho') rho' hold, so |1-rho e^{ia}|^2 <= (rho/rho') |1-rho' e^{ia}|^2.
- Take rho = r^k and rho' = r'^k, and use sum_{k<=5n, 5 nmid k} k = 10 n^2.

Use: for L in [0, 0.1), bound at L' = 0.1 and pay 0.1 n. Since theta is fixed, y is too, and p is decreasing, so
n p(L) >= n p(0.1). Pieces (E), (E0) and (N) therefore require the stronger inequality on every box touching L = 0.1.
Piece (F) needs no such reduction.

**Closed forms. [standard]** With the principal Li2, which is valid because |e^{-s}| < 1 along the whole path:
- For h != 0: p(sigma) = (1/sigma) int_0^sigma Q = (Li2(e^{-5 sigma})/5 - Li2(e^{-sigma}) + 4 zeta(2)/5)/sigma. This uses sum_{c=0}^4 Li2(zeta^c x) = Li2(x^5)/5.
- For h = 0: p0(sigma) = 4(Li2(e^{-sigma}) - zeta(2))/sigma.

p is decreasing on L > 0, because p' = (Q - p)/L and Q is decreasing. Q is convex with |Q'| <= 2:
- Q'' = [g(s/2) - g(5s/2)]/s^2, where g(x) = x^2/sinh^2 x is decreasing.
- Q'(0+) = -2, and Q' rises to 0.

## 3. Piece (E): Euler-Maclaurin near the peaks, 1.1 <= |y| <= 6

This is GAP_CENTRE2 Lemma A2, with one observation added. **Phi_c is analytic on Re s > 0 and at s = 0.** Its
singularities lie on the imaginary axis away from 0. So the finite EM identity on [0, n] holds for every y. The disc
|sigma| < 2pi/5 was only needed for the Taylor/Laplace step, not for the upper bound.

    log|S_n| <= n Re p(sigma) + Re C_h(sigma) + K_h(sigma)/n,
    K_h = sum_c [ |B2(c/5)|/2 |sigma| (|Phi_c'(sigma)| + |Phi_c'(0)|) + |sigma|^2/12 int_0^1 |Phi_c''(t sigma)| dt ].

The remainder uses the integral form, (1/12) int |f''|, with the integral bounded by sums over t-pieces.

**Certification.** 20 000 (L, y)-boxes, adaptively bisected, with L in [0.1, 2.1], 1.1 <= |y| <= 6 and h in {1,2}.
- p(sigma) and C_h(sigma) use mean-value forms, p(sc) + p'(box)(sigma - sc) with p'(s) = int_0^1 t Q'(t s) dt, enclosed on t-pieces.
- |Phi'| and |Phi''| use rigorous lower bounds of |1 - w| over the rectangle.
- Each box satisfies V(N0) = -N0 g + Re C + K/N0 + E(N0) < 0, with g = p(L_hi) - Re p(sigma) - 0.1 chi.
- It also satisfies g > 1.5/N0, which makes V(n) non-increasing for n >= N0.
- Worst V is -12.0, at L near 2.1 and |y| near 1.1. The largest K is 5.6e3, at L = 0.1 near a singular crossing, where K/n is about 7.

## 4. Piece (E0): theta near 0, |y| <= 20

Here Phi(s) = log(1 - e^{-s}) is singular at s = 0, so the t = 0 terms are split off. The finite EM is applied on
t = 1..n-1, and sum_c B1(c/5) = 0:

    log|S_n| <= n Re p0(sigma) + Psi(omega) + 0.08 (|omega Phi'(omega)| + |omega||Phi'(sigma)|) + (|omega||sigma|/3) int_{1/n}^1 |Phi''(t sigma)| dt,
    Psi(omega) = 4 + log(24/625) + sum_c Re l(c omega/5) - 4 Re[(1/omega) int_0^omega l],  l(s) = log((1-e^{-s})/s).

The log|omega| terms cancel exactly.

Analytic bounds [proved], for Re s >= 0 and |s| <= 0.1:
- |(1-e^{-s})/s - 1| <= |s|e^{|s|}/2.
- Hence |l(s)| <= |s|.
- |s Phi'(s)| <= 1/(1 - |s|e^{|s|}/2).
- |Phi''(s)| <= 1/(|s|^2 m^2), with m = 1 - r e^r/2.

Consequences:
- Psi <= 0.7403 + 6|sigma|/N0.
- The part of the integral on [1/n, 1/64] contributes at most 1/(3m^2).
- The part on [1/64, 1] is enclosed on t-pieces and multiplied by |sigma|^2/(3 N0).

The first t-piece of p0' is handled through x/(e^x - 1).

Certified on 3200 initial (L, y)-boxes (4739 after bisection), with y in [0, 20] and L in [0.1, 2.1]. The worst V is
-0.36, which is a bisection artefact.

## 5. Piece (N): the j-series triangle bound, 6 <= |y| <= 0.25 n (h = 0: 20 <= |y|)

Since |q| < 1, log|S_n| = -Re sum_j G_j/j <= sum_j |G_j|/j, where

    G_j = (1 - e^{-j sigma}) A_j/(1 - e^{-j omega}),  A_j = sum_c zeta^{hcj} e^{-cj omega/5},

with the following bounds on A_j:
- |A_j| <= min(4, Abar_j + 2j|omega|) and |A_j| <= alpha(jL/n), where alpha(a) = sum_c e^{-ca/5}.
- Abar_j = |sum_c zeta^{hcj}| is 1 or 4 (h != 0), and 4 (h = 0).

**Part I (j <= pi n/y).** Here 1/|1-e^{-z}| <= 1/|z| + c1 with c1 = 1 + 1/pi. This is the maximum modulus of
1/(1-e^{-z}) - 1/z on [0,20] x [-pi, pi]:
- On the left edge the value is at most 0.593.
- On the top and bottom edges it is at most 1 + 1/pi.
- On the right edge it is at most 1.06.

So Part I <= n M(sigma) + R1, with:
- M(sigma) = sum_j |1-e^{-j sigma}| Abar_j/(j^2|sigma|).
- R1 <= B(1 + log(pi n/y)) + O(1), with B = 2 + 1.6 c1.

**Part II (j > pi n/y).** Split into periods k >= 1, where jy/n is in ((2k-1)pi, (2k+1)pi]. Then:
- 1/j <= y/((2k-1)pi n).
- |1 - e^{-j omega}| >= max(1 - e^{-a}, 2e^{-a/2}|sin(b/2)|).
- |sin(b/2)| >= |b - 2 pi k|/pi.
- The two j nearest the pole are bounded by G_k = 1/(1-e^{-a_k^-}).
- The others use sum_{i<=M} min(G, D/i) <= min(MG, D(1 + log+(2MG/D))).

Three consequences:
- The pole terms sum to at most beta u^2 n/L, with u = y/n and beta = 1 + e^{-pi L/u}.
- The rest is O(log^2 n).
- Periods with a_k^- > 2 give a constant tail, with int_2^inf phi(a)/a da in closed form.

**Certification.**
- (N-a): y in [YE, 200], on (L, y)-boxes. It uses the exact M(sigma), 400 terms plus a tail, evaluated at n = N0. Monotonicity holds because Part II is non-increasing in n, R1 = A + B log n, and g >= (B + 1.5)/N0 is checked.
- (N-b): y in [200, 0.25n], on (L, u)-boxes with u = y/n. The bound has the explicit form -nG + a + b log n + c log^2 n, and G >= (b + 2c log N0)/N0 is checked.

Worst certified values (V at N0):

| h | (N-a) | (N-b) |
|---|---|---|
| 1 | -147.9 | -23.3 |
| 0 | -71.3 | -1.8 |

Floating point says the true slack of the triangle bound is larger than this. For example, (TB - n p)/n is about
-0.6 at u = 0.25, L = 0.1, n = 821.

Why (N) stops at u = 0.25: the pole term beta u^2/L must stay below p(L) - 0.1 at L = 0.1. The TB majorant blows up
exactly at low-denominator rationals such as theta = pi and 2pi/3, which is where the chunk method takes over.

## 6. Piece (F): chunks and the maximum principle, dist >= 0.0499

T = 40 is used on the bands 0.05 <= dist < 0.1, where a chunk must rotate by 5T x dist >= 2pi. T = 20 is used elsewhere.

**Chunk lemma. [proved]** Let M = floor(n/T), lambda = TL/n and mu = L/(5n) <= 2.1/(5N0). Define
P(V) = prod_{k'<=5T, 5 nmid k'}(1 - V e^{i k' thetat}), with thetat = theta + i mu, so that the drift r^{k'} is an
imaginary shift of theta. Then:

    log|S_n| <= sum_{m<M} B(m lambda) + 4(T-1) log(1 + e^{-(L-lambda)}),   B(s) = max_{|V|=e^{-s}} log|P(V)|,
    n p(L) >= sum_{m<M} T Q(m lambda) - T[lambda + (log 5 - Q(L))/2]   (tangent line + convexity telescoping).

**Reductions. [proved]**
- B is convex and non-increasing in s (Hadamard three circles). So on an s-grid of step 0.1, B <= chord and TQ >= tangent. This gives cell lower bounds d_i for D0 = TQ - B.
- sup over V of log|P| is subharmonic in thetat. So on R = [th_a, th_b] x [0, mu_max], the sup is taken on the boundary: the bottom edge (mu = 0), the top edge (mu = mu_max) and the two vertical sides.
- **The radial drift therefore needs no treatment in the interior.**

**Assembly. [arb]** Let dt be the running minimum of d_i, which is non-increasing. Left Riemann sums then give
sum_m D0(m lambda) >= (n/(T L_b)) int_0^{L_a(1-T/N0)} dt. For small L, M d_min is used instead. Both are checked on
L-intervals of width 0.01 over (0, 2.1], at n = N0, together with the slope condition, so the bound holds for all
n >= N0.

**Certification.** On each boundary part, B(s_i) <= t_i is shown by branch-and-bound over cells (theta-interval,
V-arc, mu-interval). There are 163 rectangles and 22 s-values. The minimum assembled margin at N0 is 47.9, near
theta = pi with T = 20. The tightest T = 40 band, near theta = 0.05, has margin 51.8.

## 7. What is and is not rigorous here (honest list)

- **(E), (E0), (N)** are arb ball arithmetic over full parameter boxes. Their only analytic inputs are the lemmas quoted above, all proved in the text.
- **(F): not arb, but float64 interval arithmetic.** The inner loop uses only IEEE +, -, *, / (correctly rounded) and numpy frexp/ldexp, which are exact.
  - All transcendental inputs come from arb: e^{i theta_c}, e^{i alpha_c}, e^{-s-k'mu} and the thresholds e^{2t}.
  - Explicit a-priori rounding margins are added:
    - 1e-12 on each phase half-width. The worst-case error of the iterated powers is at most 4.3 k' u, which is at most 1e-13.
    - 1e-12 on cos.
    - 32u absolute on each |1 - rho e^{i phi}|^2.
    - A (1+2u)^{2(K+3G+8)} factor on the product.
  - No libm function enters a bound.
  - `offpeak_crosscheck.py` checks the cell evaluator against arb on 200 random point cells: 0 unsound cases. It also checks the (E), (E0) and (N-a) box bounds against exact log|S_821| at 80 random points: 0 violations.
  - If a pure-arb certificate is required, re-running `Chunk.check` with acb balls is a mechanical change. The cost is about 50-100x the 470 s.
- The chunk argument covers every L in (0, 2.1], including L -> 0, directly. L = 0 follows as the limit case, since both sides are continuous and the left side may be -infinity.
- **Not re-checked here:** piece (A) (GAP_CENTRE2) and the exact runs for n <= 820. The claim that (B) with 6.3 is what (A)'s normalization needs is taken from GAP_CENTRE2 section 3.

## 8. Consequence

(B) holds for all n >= 821 with constant 6.3, and a fortiori with 5.6 for n >= 1001. Together with (A) for n >= 821
and the exact signs for n <= 820, the centre regime has no remaining open piece, subject only to the float-interval
caveat for (F). The running 820..920 exact job is no longer needed for the centre.
