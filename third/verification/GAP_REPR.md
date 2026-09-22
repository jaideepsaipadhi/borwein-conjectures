# GAP_REPR: the representation error L(m) - (one-period integral)

Scripts: `repr_D0.py` (certified kernel constant), `repr_check.py` (+ `repr_check.log`, brute force),
`repr_tile.py`, `repr_tile2.py`, `repr_tile3.py` (+ `.log`, tiling with the new term).

## 1. Exact identity: fold the kernel, not the integrand

The kernel has a closed-form fold. Let w = W + it, w_k = w + 2 pi i k, om = e^{2 pi i (m-1/6)} = e^{-i pi/3},
r_1(u) = e^{a5/u} - 1 - a5/u. Then

    Phi(w) := sum_{i>=1} P_5(i) e^{-(i-1/6) w}
            = a5 e^{-5w/6}/(1 - e^{-w}) + sum_{k in Z} r_1(w_k) om^k        (absolutely convergent)
            = e^{a5/w} - 1 + eps(w),
    eps(w)  = a5 g(w) + sum_{k != 0} r_1(w_k) om^k,     g(w) = e^{-5w/6}/(1-e^{-w}) - 1/w,

and for every W > 0, **exactly**,

    L(m) = (1/2pi) int_{-pi}^{pi} Phi(w) F_m(e^{-w}) e^{(m-1/6) w} dt.                         (*)

*Proof.* Start from section 12 with K = 1. |F_m| <= 4 E_n(e^{-W}) on the line (see section 3 below), so
Fubini lets the whole-line r_1 integral fold onto one period, which gives sum_k r_1(w_k) om^k. The head
a5 sum_{j: m-j>1/6} e(j) equals (1/2pi) int_{-pi}^{pi} a5 e^{-5w/6}/(1-e^{-w}) F_m e^{(m-1/6)w} dt,
because e^{-5w/6} e^{(m-1/6)w}/(1-e^{-w}) = sum_{n>=1} e^{(m-n)w}, and the constant term picks out
exactly the j <= m-1. So **the head sum is absorbed into the folded kernel**, and no separate head term is
left. The closed form sum_k om^k/w_k = e^{-5w/6}/(1-e^{-w}) (the Mittag-Leffler / Poisson form) gives the
second expression. (*) is just Cauchy's formula for [q^m] of (sum_i P_5(i) q^i) F_m(q), so it is a
consistency check, not a new claim.

*Checks.* The kernel alone reproduces P_5(M) for M >= 1, and 0 for M <= 0, to about 1e-11 (the residual is
the |k| <= 3000 truncation). The full (*) matches L(m), summed directly from Bessel values, to within
3e-8 to 5e-13 relative error at 9 points (n, m), n = 4..60 (`repr_check.log`, column relerr).

## 2. The error term and its bound

The certificate integrates e^{a5/w} over one period. Its error is therefore

    E'(m) = L(m) - (1/2pi) int_{-pi}^{pi} e^{a5/w} F_m e^{(m-1/6)w} dt
          = (1/2pi) int_{-pi}^{pi} (-1 + eps(w)) F_m(e^{-w}) e^{(m-1/6)w} dt.

**Certified** (`repr_D0.py`, arb, maximum principle on [-0.25, 0.25] x [-pi, pi], valid for all
0 < W <= 0.02, which contains the whole Laplace region since W = a5/(tv) <= a5/20):

    sup|g| <= 0.44525,   sum_{k!=0}|r_1| <= (a5/pi)^2 e^{a5/pi} pi^2/8 = 0.009416,
    sup|eps| <= 0.12660,   D0 := sup|-1+eps| <= 1.12601.

This gives

    |E'(m)| <= D0 * 4 E_n(e^{-W}) e^{(m-1/6)W}.

The certificate's Z3 = (4/(K0 e^{-v})) e^{-a5/W} 2pi/win in `B_ball` is exactly this bound with D0 = 1.
It covered the -1 and dropped eps.

## 3. Audit item 9 is wrong: sup|F_m| <= 4 E_n(e^{-W}) holds

E_n(q) = prod_{k>5n, 5 nmid k} (1-q^k)^{-1} (NOTES section 1) has **non-negative coefficients**. So
|E_n(z^h q)| <= E_n(|q|), and |F_m| <= 4 E_n(e^{-W}) on the whole line. The audit's "4 prod(1+e^{-kW})" and
"E_n(e^{-W}) < 1" describe the *finite* product prod_{k<=5n}(1-q^k). `subtract_Fm.py` uses that finite
product as a stand-in polynomial, which is fine for a bookkeeping test but is not E_n. Measured
sup_t |F_m| / (4 E_n(e^{-W})) = 0.04 to 0.25 at every test point.

## 4. The real problem is the normalisation, and it is a new finding

Divided by the **true** main term, the bound is not small in the way `B_ball` assumes. The thin-factor
lemma bounds |F_m(e^{-W})| >= K0 e^{-v} **|E_n(z5 e^{-W})|**. The measured ratio
|F|/(|E_n(z5 q)| e^{-v}) is 1.003 to 1.027 at four points. The lemma does not bound F_m relative to
E_n(e^{-W}); that ratio is 4e-5 to 2.6e-2 at the same points. The two quantities differ by e^{Dz}, and

    Dz = log E_n(e^{-W}) - log|E_n(z5 e^{-W})|  <=  Li2(x)/W + 5 Li1(x),   x = e^{-v}

(proved by the block argument: sum_{c=1..4}(1 - cos 2 pi c r/5) = 5, and 1/(1-e^{-y}) <= 1/y + 1). The
measured Dz is 6.28 against a bound of 6.6 at n = 40. **`B_ball`'s logmain uses log E_n(e^{-W})**, so its
main term is too large. Measured: `B_ball`'s main is 130x the true |L(m)| at (40, 1624) and 7700x at
(60, 3804). Every term with an absolute positive-coefficient majorant in its numerator must therefore
carry e^{Dz}: Z_E, and also the existing Z3 and Z_R. I have not checked whether Z2's C(v) ratio is
affected in the same way. It must be checked.

Corrected term (decreasing in t at fixed v, since it is W^{-3/2} e^{-(a5 - Li2(x))/W} and a5 > Li2(e^{-2.5})):

    Z_E = D0 e^{-W/6} e^{Dz} (4/(K0 e^{-v})) e^{-a5/W} 2pi/win.

## 5. Numbers (ball arithmetic, 16 sub-cells, min over c in {2, 2.5, 3, 3.5})

At the PROOF.md table thresholds t_min (`repr_tile2.log`):

| v | t | B | Z_E (corrected) | B + Z_E | also correcting Z3, Z_R |
|---|---|---|---|---|---|
| [2.5, 3.0] | 8 | 0.832 | 0.555 | **1.388 FAIL** | 380 |
| [3.0, 3.6] | 8.25 | 0.357 | 1.8e-3 | 0.359 | **3.37 FAIL** |
| [3.6, 4.5] | 6 | 0.688 | 4.6e-3 | 0.692 | **2.94 FAIL** |
| [4.5, 6.0] | 5 | 0.750 | 1.2e-3 | 0.751 | **1.065 FAIL** |
| v >= 6 (all rows) | t_min | <= 0.43 | <= 1.1e-2 | <= 0.455 | <= 0.455 |
| tail cell, v >= 20 | 2 | 0.0660 / 0.2969 | 1.3e-4 / 1.2e-3 | 0.0662 / 0.2981 | same (e^{Dz} = 1.000) |

The thresholds are much lower than the proof needs. On the proof's range N >= 1456 and t = a5 N/v^2, so
t >= t_cov = a5 * 1456 / v_hi^2 (42.6 on [2.5, 3], 29.6, 18.9, 10.6, 7.8, 5.3, 3.8 on the next rows). At
t_eff = max(t_min, t_cov) (`repr_tile3.log`):

| v | [2.5,3] | [3,3.6] | [3.6,4.5] | [4.5,6] | [6,7] | [7,8.5] | [8.5,10] | [10,13] | [13,16] | [16,20] |
|---|---|---|---|---|---|---|---|---|---|---|
| t_eff | 42.57 | 29.56 | 18.92 | 10.64 | 7.82 | 5.30 | 3.83 | 3.00 | 2.25 | 2.00 |
| Z_E (corr.) | 3e-25 | 5e-25 | 2e-20 | 1e-13 | 2e-13 | 6e-9 | 2e-6 | 1e-4 | 4e-3 | 6.5e-3 |
| B, all corrected | 0.081 | 0.093 | 0.139 | 0.252 | 0.148 | 0.226 | 0.225 | 0.377 | 0.326 | 0.327 |

Z_E and the corrected Z_R keep decreasing for larger t. At v = 2.5, e^{Dz} Z_R = 8e-14, 3e-32 and 1e-144
at t = 40, 80 and 320.

**Brute force** (`repr_check.log`): the bound D0 * 4 E_n(e^{-W}) e^{(m-1/6)W} exceeds the true |E'(m)| by
a factor of 39 to 4600 at all 9 test points. In the Laplace region, |E'|/|L| = 1.5e-5 at (40, 1624; v = 2.56,
t = 8.1) and 1.5e-8 at (60, 3804; v = 2.50, t = 12.6).

## 6. Head sums

With (*) there is no head term: it is folded into Phi. If you keep the K = 1 form anyway, the head is
|a5 sum_{j<m} e(j)| <= 4 a5 sum_{j<m} [q^j]E_n <= 4 a5 E_n(e^{-W}) e^{(m-1)W} (Rankin). Relative to the true
main term this is a5 e^{W/6}/(2pi D0) times Z_E, which is negligible wherever Z_E is. For K = 2 the extra
factor (m-j-1/6) <= m adds a factor a5 m/2 = a5^2/(2W^2) to the bound. That is still exponentially
dominated by e^{-(a5-Li2)/W}, but it is pointless given (*).

## Verdict

**Closed with caveats.**
- The identity (*) is exact, and the error term is rigorously bounded by D0 = 1.12601 (certified) times the
  valid majorant 4 E_n(e^{-W}).
- Relative to the true main term it needs the factor e^{Dz}. With that factor, B stays below 1 in every
  row once each row is evaluated at the t forced by N >= 1456. Worst case: B + Z_E = 0.377 at [10, 13];
  the tail cell gives 0.0662 or 0.2981.
- It does **not** stay below 1 at the PROOF.md t_min for v in [2.5, 3.0].

Caveats:
1. `B_ball` normalises its main term by E_n(e^{-W}) instead of |E_n(z5 e^{-W})|, which is wrong by up to
   ~1e3 at the table thresholds and 1e15 at t = 42. That is outside this gap, but it invalidates the Z3/Z_R
   values at the table t_min for v < 6. It is repaired here only by using t_eff, which has to be stated in
   PROOF.md section 5, with the monotonicity-in-t argument.
2. Whether Z2's C(v) carries the same missing e^{Dz} factor is unchecked.
3. B_ball's Z_R omits a 2pi relative to the true main term (e^{Rh} win/(2pi)). This does not matter at t_eff.
