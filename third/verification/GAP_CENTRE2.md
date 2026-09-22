# GAP_CENTRE2: the centre regime (mu in [0.105, 1/2]), second pass

*18 Sept 2026. This pass certifies two things and leaves one piece open. The centre is now reduced to a
single explicit off-peak inequality (B) for n >= 821.*

## Summary

| piece | range | status | evidence |
|---|---|---|---|
| exact signs, **all** m, all classes | 291 <= n <= 820 | **[exact]**, 0 wrong signs (820 < n <= 1000: job running or left to the user, see section 1) | `centre3_exact.py` + logs |
| (A) near-peak, \|y\| <= 1.1 | n >= 821 (also run at 1001), L in [0, 2.1] (all mu in [0.104, 1/2]), all 5 classes | **[certified]**, ball arithmetic over full boxes | `centre3_peak.py`, `centre3_peak_821.log`, `centre3_peak_1001.log`, `centre3_saddle_end.py` |
| (B) off-peak, \|y_h\| > 1.1 for all h | n >= 821 | **[open]**: the only remaining gap in the centre | — |

Worst err/margin in (A), counting peak plus strip only (each term is non-increasing in n):
- N0 = 821: class 0: 0.057, class 1: 0.223, class 2: 0.224, class 3: 0.798, class 4: 0.847 (L ~ 2.1).
- N0 = 1001: 0.046, 0.182, 0.184, 0.620, 0.658.

## 1. Exact computation (`centre3_exact.py NA NB`)

- S_n is palindromic of degree 10n^2, with c(D-m) = c(m) and D = 0 mod 5 (so classes s and -s swap). Checking m <= 5n^2 is therefore enough.
- The script works mod x^{5 NB^2 + 1}. It builds S_NA incrementally from S_0 = 1, then gets each next n by four exact shift-subtracts, S <- S - x^k S for k = 5n+1..5n+4.
- It checks every coefficient: a strict wrong sign counts as a violation. It also lists zeros, which occur only at m <~ 5n (mu ~ 1e-3). Those lie in the band m < 4N, far from the centre, and matter only to the band/Laplace pieces, which treat sign non-strictly just as `verify_small.py` does.
- Done, with 0 wrong signs: 291..410 (`centre3_exact_291_410.log`, 114 s), 410..580 (419 s) and 580..820 (1357 s).
- 820..920 was launched in background (`centre3_exact_820_920.log`). It is followed by 920..1000, which may exceed 8 GB RAM, so check that log for `DONE`.
- The build is incremental from n = 0 because a mul_low product tree ran out of memory at n ~ 820.
- Commands for the user, run on a machine with 16 GB or more; each is about 1-1.5 h:
  `python3 centre3_exact.py 820 920 > centre3_exact_820_920.log`
  `python3 centre3_exact.py 920 1000 > centre3_exact_920_1000.log`
- Until those print `DONE ... total violations 0`, the exact range is n <= 820. Piece (A) passes at N0 = 821 (`centre3_peak_821.log`), so the handover at 820/821 needs no further exact runs. The 820-1000 runs only enlarge the slack available for (B).
- Memory is about 4 x 5 NB^2 x (2.3 NB bits), measured at 3.4 GB RSS for NB = 820. So NB = 1000 needs about 6 GB, NB = 1400 about 16 GB and NB = 2000 about 45 GB. Run it in segments with NB <= 1.4 NA:
  `python3 centre3_exact.py 1000 1400 > centre3_exact_1000_1400.log` (about 2-3 h and 12 GB RAM, estimated from the n^3 step cost).

## 2. Piece (A): certified near-peak lemma, n-uniform for n >= N0 (run at N0 = 821 and at 1001)

Notation:
- q = zeta^h e^{-sigma/(5n)} with sigma = L - i y, so y = 5n(theta - 2 pi h/5) and r = e^{-L/(5n)}.
- Phi_c(s) = log(1 - zeta^{hc} e^{-s}), Q = sum_{c=1..4} Phi_c, and p(sigma) = int_0^1 Q(t sigma) dt.
- C_h(sigma) = sum_c (c/5 - 1/2)(Phi_c(sigma) - Phi_c(0)).

**Lemma A1 (finite Euler-Maclaurin). [proved; identity checked numerically to 1e-15]** For 0 <= a < 1:

  sum_{t=0}^{n-1} f(t+a) = int_0^n f + B1(a)(f(n)-f(0)) + (B2(a)/2)(f'(n)-f'(0)) - int_0^n (B~2(x-a)/2) f''.

- Constants: B2(1/5) = B2(4/5) = 1/150 and B2(2/5) = B2(3/5) = -11/150, so the coefficient is at most 11/300. Since B~2 lies in [-1/12, 1/6], the last term is at most (1/12) int |f''|.
- The 11/300 constant in GAP_CENTRE is therefore correct.
- The script uses the per-c values |B2(c/5)|/2 rather than the uniform 11/300.

**Lemma A2. [proved + certified constants]** Apply A1 to f(x) = Phi_c(x sigma/n) with a = c/5, and use int_0^n f = (n/sigma) int_0^sigma Phi_c. Then

  |log S_n(q) - n p(sigma) - C_h(sigma)| <= K_h(sigma)/n,
  K_h = sum_c [ (|B2(c/5)|/2)|sigma|(|Phi_c'(sigma)| + |Phi_c'(0)|) + (|sigma|^2/12) sup_{t in [0,1]} |Phi_c''(t sigma)| ].

This holds for |y| < 2pi/5 and L >= 0, uniformly down to L = 0: no singularity of Phi_c lies on [0, sigma]. (The main-term constant is B1(c/5) = c/5 - 1/2, which gives C_h.) A float check at n = 50 and 200 confirms the O(1/n) behaviour. Certified sup of K over the window: 7.3 at L ~ 0, and 0.75 to 1.4 elsewhere.

Facts used: |e^{C_h(L)}| = 1 for real L (enclosed, not assumed), and the four peaks have equal modulus exactly. Class margins are M_s(L) = sgn_s Re sum_h zeta^{-hs} e^{C_h(L)}. This replaces the arctan law by an exact enclosure, and it matches the old binding value 0.126 (class 4, L ~ 2.08).

**Laplace step. [certified]** With the saddle condition p'(L) = -2mu:

  c(m) = N(n,L) [ sum_h zeta^{-hm} e^{C_h(L)} (1 + rho_h) + strip + OFF ],
  N = e^{n p(L) + 2n mu L} sqrt(2pi/(n a))/(10 pi n),  a = p''(L).

- **Inner window |y| <= 0.35.** Write psi(y) = p(L-iy) - p(L) + i p'(L) y and Z = n psi + n a y^2/2 + DeltaC. The odd part T = n psi_3 y^3 + C_1 y integrates to exactly 0 against the Gaussian. The remaining error uses |e^Z - 1 - Z| <= |Z|^2/2 max(1, e^{Re Z}) and Re psi <= -b y^2/2, where b is a certified lower bound for Re p'' on the window. Taylor constants come from p^{(k)}(sigma) = int t^k Q^{(k)}(t sigma) dt, enclosed as complex balls through the convex-hull mean value on t-pieces. The bound also carries e^{K/n} - 1 and the Gaussian tail. Every term has the form c/n, (e^{K/n}-1) or a decreasing exponential, so a bound at n = N0 holds for all n >= N0.
- **Strip 0.35 <= |y| <= 1.1.** Lemma A2 plus Re psi(y) = -int_0^{|y|} (|y|-u) Re p''(L-iu) du, with Re p'' bounded below on boxes, give sqrt(n)e^{-nD} terms. These are non-increasing for n >= 1/(2D), which is asserted in the script.
- **Grid.** 210 L-boxes of width 0.01 on [0, 2.1]. Each has 28 inner and 60 strip y-boxes and 64 t-pieces; every quantity is enclosed over its full box. `centre3_saddle_end.py` certifies mu(2.1) = -p'(2.1)/2 < 0.105, so every centre mu has its saddle in [0, 2.1]. Positivity a > 0 is asserted on each box, so the saddle is unique.
- **Result.** M_s > sum_h |e^{C_h}|(rho_h + strip_h) on every box for every class at n = N0, hence for all n >= N0. The minimum slack M_4 - err is about 0.019 at N0 = 821 and about 0.042 at N0 = 1001 (L ~ 2.1); elsewhere it is larger.

## 3. Piece (B): what remains, stated exactly

Covering OFF needs the following. For n >= 821, L in [0, 2.1] and every theta with |5n(theta - 2pi h/5)| > 1.1 for h = 1..4 (mod 2pi):

  log|S_n(e^{-L/(5n)} e^{i theta})| - n p(L) <= -(1.5 log n + 6.3).                      (B)

This makes |OFF| <= 10 pi n sqrt(n a/2pi) e^{-(...)} <= 0.0187 <= min(M_s - err) at N0 = 821 (worst at L ~ 2.1, class 4). If the exact range reaches 1000, (B) is needed only for n >= 1001, where 1.5 log n + 5.6 suffices. An L-dependent right-hand side read from the logs is weaker and would also do.

Float evidence (`centre2_landscape.py`, not a proof): the true deficit for |y| >= 1 is at least 0.06n, which is >= 49 at n = 821 against the required 16.4.

Why this pass did not prove (B): I estimated a rational-approximation / residue-class Euler-Maclaurin argument (major/minor arcs with O(1) cost per class and per singularity crossing). With explicit constants it closes only for n of order 1e4 to 1e5, which is beyond the exact range reachable here (n ~ 1000 at 7 GB; about 2000 at 35 GB). The chunk + maximum-principle route has the defect recorded in GAP_CENTRE near the 5th roots. Neither route was certified.

Two directions look viable:
- (i) Raise the exact range to ~2000 on a large machine, and prove (B) with the residue-class EM argument using sharp constants.
- (ii) Prove (B) only for n >= some n_B, by the same boxing technique as (A) applied to the continuum at rationals a/b. Its error is O(b + #crossings), and n_B must come out at most the exact limit.

## Status of old audit items
- C.2.4: constant verified (11/300; the script uses per-c values).
- C.2.5 and C.2.6: fixed. Derivatives are enclosed over whole windows and the arctan/equal-modulus laws are replaced by exact enclosures.
- C.2.7 (L -> 0): resolved; the y-variable is used and L = 0 is inside the certified range.
- C.2.1-3 (off-peak tail): open. This is (B).

**Exact: the centre (indeed every coefficient) is proved for n <= 820. (A) is certified for all n >= 821. The centre remains OPEN for n >= 821, and exactly (B) is missing.**
