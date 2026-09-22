# §3 (full proof). The explicit expansion of G_5

*Standalone write-up of the transplanted circle-method argument, 18 September 2026. It replaces the citation to `GAP_CITATION.md` in PROOF.md §3.*

**Notation.** In this section the coefficient index is written **n**, not i, so that i always means √−1. The statement is about the coefficient index, so here n plays the role of the i used in PROOF.md.

- e(x) = e^{2πix}.
- F(x) = f(x) = 1/(x;x)_∞ = ∏_{m≥1}(1−x^m)^{−1} for |x| < 1.
- s(h,k) = Σ_{r=1}^{k−1} (r/k)((hr/k)) is the Dedekind sum, with ((x)) = x − ⌊x⌋ − 1/2 for x ∉ ℤ and ((x)) = 0 for x ∈ ℤ.
- ω(h,k) = e^{πi s(h,k)}.
- φ is Euler's function.

Standard facts that are used without proof are tagged **[standard: reference]**.

---

## 3.1 Statement

Let G_5(q) = (q;q)_∞/(q^5;q^5)_∞ = Σ_{n≥0} s(n) qⁿ. Put

- B = 1/3 and C = C(n) = 2n − 1/3;
- P_k(n) = (2π/k)·√(B/C)·I_1((2π/k)√(BC));
- K(n) = ⌊√(4πn)⌋ + 2;
- a(n) = Σ_{h=1}^{4} w_h ζ^{−nh}, where w_h = e^{−πi s(h,5)} and ζ = e(1/5).

**Theorem 3.1.** For every n ≥ 1,

  |s(n) − a(n)P_5(n)| ≤ Σ_{5|k, 10≤k≤K(n)} φ(k)P_k(n) + E_arc + E_rem + E_ng,

where

- E_arc = √2·π·(4/25)·e^{1+π/3} ≤ 5.506447,
- E_rem = 2√2·(4/25)·e^{1+π/3}·(f(t)f(t⁵) − 1) ≤ 0.006571, with t = e^{−2π},
- E_ng = 2√2·e·K_5 ≤ 21.950524, with K_5 = √5·e^{−π/15}·f(e^{−2π/5})·f(e^{−2π}) ≤ 2.8549955.

Their sum is E_con ≤ 27.463542 ≤ **27.4636**.

The general-p formulas specialise correctly at p = 5: (p−1)/p² = 4/25, π(p−1)/(12p) = π/15, and B = (p−1)/12 = 1/3.

Two remarks on the statement:

- **Every P_k(n) is positive.** For n ≥ 1 we have C ≥ 5/3 > 0, and I_1(x) > 0 for x > 0.
- **Closed form of a(n).** Section 3.9 verifies s(1,5) = 1/5, s(2,5) = s(3,5) = 0 and s(4,5) = −1/5. Hence w_1 = e^{−πi/5}, w_2 = w_3 = 1 and w_4 = e^{πi/5}, so

  a(n) = 2cos(2πn/5 + π/5) + 2cos(4πn/5).

  For n ≡ 0, 1, 2, 3, 4 (mod 5) this equals (5+√5)/2, −√5, −(5−√5)/2, 0, 0.

---

## 3.2 Choice of the Farey order

Take **N = K(n) = ⌊√(4πn)⌋ + 2**. Two properties of this choice are used.

**(N1) N ≥ 5 for every n ≥ 1.** Since √(4π) = 3.5449…, we have N ≥ 3 + 2 = 5. So the denominator k = 5 always occurs in the dissection, which is needed to produce the main term a(n)P_5(n).

The choice N = K(n) − 1, suggested as a possibility in the brief, **fails** at n = 1: it gives N = 4, and then no fraction with k = 5 occurs. With N = K(n), the sum over 5 | k, 10 ≤ k ≤ N in the error term is exactly the sum in the statement.

**(N2) N + 1 > √(4πn), so 2πC/(N+1)² < 1.** This holds because N + 1 = ⌊√(4πn)⌋ + 3 > √(4πn) + 2. Then

  2πC/(N+1)² = (4πn − 2π/3)/(N+1)² < 1.

This also satisfies the source's stronger requirement N ≥ √(4πn) + 1.

---

## 3.3 The Rademacher path

**Lemma 3.2 (Farey dissection)** [standard: Apostol, *Modular Functions and Dirichlet Series in Number Theory*, 2nd ed., Ch. 5, §§5.4–5.7; Rademacher, *Topics in Analytic Number Theory*, Ch. 14]. Let N ≥ 1, and let h_1/k_1 < h/k < h_2/k_2 be consecutive fractions in the Farey sequence of order N, taken cyclically on [0,1). Let C(h/k) be the Ford circle with centre h/k + i/(2k²) and radius 1/(2k²). Let γ_{h,k} be its upper arc between its points of tangency with C(h_1/k_1) and C(h_2/k_2). Then:

1. The arcs γ_{h,k} join up into a path from some τ_0 to τ_0 + 1 in the upper half-plane ℍ.
2. Since G_5(e(τ))e(−nτ) is holomorphic on ℍ and 1-periodic, the coefficient formula s(n) = ∫_{τ_0}^{τ_0+1} G_5(e(τ))e(−nτ) dτ, taken along a horizontal segment, may be deformed onto this path by Cauchy's theorem.
3. The substitution τ = h/k + iz/k² maps C(h/k) onto the circle K: |z − 1/2| = 1/2.

It follows that

  s(n) = Σ_{k≤N} Σ_{0≤h<k, (h,k)=1} (i/k²)·e(−nh/k) ∫_{z'}^{z''} G_5(e(h/k + iz/k²))·e^{2πnz/k²} dz,   (3.1)

where:

- z' = (k² + ikk_1)/(k² + k_1²) and z'' = (k² − ikk_2)/(k² + k_2²);
- the integral runs along K clockwise from z' (where Im > 0) through z = 1 to z'' (where Im < 0).

This uses dτ = (i/k²)dz and e(−nτ) = e(−nh/k)·e^{2πnz/k²}. The tangency point of C(h/k) with C(h_1/k_1) is h/k − k_1/(k(k²+k_1²)) + i/(k²+k_1²), and it maps to z'. The same holds for z''.

**Lemma 3.3 (geometry).** With the notation of Lemma 3.2:

- **(b)** |z'|, |z''| ≤ √2·k/(N+1), and Re z', Re z'' ≤ 2k²/(N+1)².
- **(c)** On the chord s_{h,k} = [z', z''] we have Re(1/z) ≥ 1 and 0 < Re z ≤ 2k²/(N+1)². Also |s_{h,k}| ≤ 2√2·k/(N+1).
- **(d)** On the arcs of K from z'' to 0 and from 0 to z' (the complement of the Rademacher arc), Re(1/z) = 1 and Re z ≤ 2k²/(N+1)². Each of these arcs has length at most (π/2)·√2·k/(N+1).

*Proof.*

(b) Consecutive fractions in the Farey sequence of order N satisfy k + k_j ≥ N + 1 [standard: Hardy–Wright, *An Introduction to the Theory of Numbers*, Thm 30]. Hence k² + k_j² ≥ (k+k_j)²/2 ≥ (N+1)²/2. Therefore

- |z'| = k/√(k² + k_1²) ≤ √2·k/(N+1), and
- Re z' = k²/(k² + k_1²) ≤ 2k²/(N+1)².

The same argument applies to z''.

(c) The map z ↦ 1/z sends K∖{0} onto the line Re w = 1, and sends the open disc bounded by K onto the half-plane Re w > 1. The chord lies in the closed disc and avoids 0, so Re(1/z) ≥ 1 on it. Re z is affine along a segment, so it is at most its value at an endpoint, and it is positive at both endpoints. Finally |s_{h,k}| ≤ |z'| + |z''|.

(d) Re(1/z) = 1 on K∖{0}. The point z' lies on the upper semicircle of K, which runs from 1 to 0. Let θ ∈ (0, π] be the central angle of the arc from z' to 0 along that semicircle. The diameter of K is 1, so the chord from 0 to z' has length |z'| = sin(θ/2), and the arc has length θ/2. By Jordan's inequality, θ/2 ≤ (π/2)·sin(θ/2) for θ/2 ∈ [0, π/2]. Hence the arc has length at most (π/2)|z'|. On the semicircle, Re z = (1 + cos ψ)/2 is monotone, so on this arc Re z ≤ Re z'. The lower arc from z'' to 0 is handled the same way. ∎

**Corollary 3.4 (the exponential factor).** On every chord s_{h,k} and on every arc in (d),

  |e^{πCz/k²}| = e^{πC·Re z/k²} ≤ e^{2πC/(N+1)²} < e.

This follows from Lemma 3.3 and (N2).

---

## 3.4 The modular transformation at a cusp h/k

**Fact 3.5 (transformation of the partition generating function)** [standard: Apostol, op. cit., Thm 5.1; this is the Dedekind η transformation in Rademacher's normalisation]. Let k ≥ 1, (h,k) = 1, Re z > 0, and let H be an integer with hH ≡ −1 (mod k). Put

- x = exp(2πih/k − 2πz/k²),
- x' = exp(2πiH/k − 2π/z).

Then

  F(x) = ω(h,k)·(z/k)^{1/2}·exp(π/(12z) − πz/(12k²))·F(x').

Here the square root is principal. Equivalently,

  (x;x)_∞ = ω(h,k)^{−1}·(z/k)^{−1/2}·exp(−π/(12z) + πz/(12k²))·(x';x')_∞.

Note that this is the congruence hH ≡ **−1**, the Hardy–Ramanujan convention. The source records the same correction.

Take τ = h/k + iz/k², so that x = e(τ) and G_5(x) = (x;x)_∞·F(x⁵). There are two cases, according to d = gcd(5,k).

### Case d = 5 (growing cusps), k = 5k_1

Since gcd(h,k) = 1, also gcd(h,k_1) = 1. Write

  x⁵ = exp(2πih/k_1 − 2πz_1/k_1²),  with z_1 = 5z·k_1²/k² = z/5,

so Re z_1 > 0. Apply Fact 3.5 to (h, k_1, z_1), with H_1 chosen so that hH_1 ≡ −1 (mod k_1). In that application:

- the square-root factor is (z_1/k_1)^{1/2} = (z/k)^{1/2};
- π/(12z_1) = 5π/(12z);
- πz_1/(12k_1²) = 5πz/(12k²).

Multiplying the result by the formula for (x;x)_∞ gives

  G_5(x) = Φ_{h,k}·exp(π/(3z) − πz/(3k²))·Ĝ_{h,k}(z),   (3.2)

where

- Φ_{h,k} = ω(h,k_1)/ω(h,k) = exp(πi(s(h,k/5) − s(h,k))), so |Φ_{h,k}| = 1;
- Ĝ_{h,k}(z) = (x';x')_∞·F(x_5'), with x_5' = exp(2πiH_1/k_1 − 10π/z).

The square roots cancel exactly: (z/k)^{−1/2}·(z/k)^{1/2} = 1, because both are principal powers of the same number. So **no branch choice enters** (the source's "δ = 1 is an integer" remark). The exponents combine as

  −π/(12z) + 5π/(12z) = **+π/(3z)**  and  πz/(12k²) − 5πz/(12k²) = −πz/(3k²).

For k = 5 we have k_1 = 1, s(h,1) = 0 and Φ_{h,5} = e^{−πi s(h,5)} = w_h.

### Case d = 1 (non-growing cusps), 5 ∤ k

Now gcd(5h,k) = 1. Write

  x⁵ = exp(2πi·5h/k − 2π·(5z)/k²).

Apply Fact 3.5 to (5h, k, 5z), with 5hH_5 ≡ −1 (mod k). This gives

  G_5(x) = Φ'_{h,k}·√5·exp(−π/(15z) − πz/(3k²))·(x';x')_∞·F(x_5''),   (3.3)

where

- Φ'_{h,k} = ω(5h,k)/ω(h,k), so |Φ'_{h,k}| = 1;
- x_5'' = exp(2πiH_5/k − 2π/(5z)).

Here (z/k)^{−1/2}·(5z/k)^{1/2} = √5 exactly, because z and 5z have the same argument in (−π/2, π/2). The exponents combine as −π/(12z) + π/(60z) = −π/(15z), and πz/(12k²) − 5πz/(12k²) = −πz/(3k²).

### Relation to η and the shift in C

Since η(τ) = e(τ/24)(x;x)_∞, we have G_5 = e(−τ/6)·η(τ)/η(5τ). This is the "q^{1/6}" shift. In (3.2) and (3.3) the shift appears as the factor e^{−πz/(3k²)}, which combines with the circle-method factor e^{2πnz/k²} in (3.1) to give

  e^{2πnz/k²}·e^{−πz/(3k²)} = e^{πCz/k²},  C = 2n − 1/3.

The growth exponent π/(3z) in (3.2) defines B = 1/3. At a general p, (d² − p)/(12p) gives (25 − 5)/60 = 1/3 when d = 5, and −(p−1)/(12p) = −1/15 when d = 1, in agreement with the source.

**Numerical confirmation.** Formulas (3.2) and (3.3) were checked against direct q-series evaluation, including all phases, at z = 1.3 + 0.4i for (h,k) = (1,5), (2,5), (3,10), (7,15), (1,1), (1,3), (2,7), (3,4). The relative error was ≤ 4·10^{−60} at 60 digits.

### Expansion and majorant of Ĝ

Put u = e^{−2π/z} and w = Re(1/z), so |u| = e^{−2πw}. Then

- (x';x')_∞ = ∏_m (1 − e(mH/k)u^m), and
- F(x_5') = ∏_m (1 − e(mH_1/k_1)u^{5m})^{−1}.

Expanding both gives Ĝ_{h,k} = 1 + Σ_{j≥1} g_j·u^j. Each g_j is a sum of products of roots of unity, and the sum has the same number of terms as the j-th coefficient of the product of majorant series. Coefficientwise,

  |g_j| ≤ [y^j] ∏_m (1 + y^m)·f(y⁵) ≤ [y^j] f(y)f(y⁵) =: φ_j,

using ∏(1 + y^m) ≼ ∏(1 − y^m)^{−1} coefficientwise (distinct partitions are partitions). Hence

  |Ĝ − 1| ≤ Σ_{j≥1} φ_j·e^{−2πjw}.   (3.4)

---

## 3.5 Growing cusps: the main term

By (3.1) and (3.2), the contribution of a growing cusp h/k (5 | k, k ≤ N) is

  (i/k²)·e(−nh/k)·Φ_{h,k}·[M_{h,k} + R_{h,k}],

where

- M_{h,k} = ∫_{z'}^{z''} exp(π/(3z) + πCz/k²) dz, and
- R_{h,k} = ∫_{z'}^{z''} exp(π/(3z) + πCz/k²)·(Ĝ_{h,k}(z) − 1) dz.

Both integrals run along the Rademacher arc of K.

**Fact 3.6 (Bessel integral)** [standard: DLMF 10.9.19; Watson, *Treatise on the Theory of Bessel Functions*, §6.22]. For x > 0 and c > 0,

  I_1(x) = ((x/2)/(2πi))·∫_{c−i∞}^{c+i∞} e^{s + x²/(4s)}·s^{−2} ds.

DLMF gives this with a Hankel loop around (−∞, 0]. For integer order the loop may be opened to the vertical line, because the integrand is O(|s|^{−2}) there and e^s decays on the left.

**Lemma 3.7.** Let B, C > 0 and let K be oriented negatively (clockwise). Then

  (i/k²)·∮_K exp(πB/z + πCz/k²) dz = P_k(B,C) := (2π/k)·√(B/C)·I_1((2π/k)√(BC)).

*Proof.*

1. Put w = 1/z = 1 + iτ with τ ∈ ℝ. The clockwise traversal of K∖{0} corresponds to τ increasing from −∞ to ∞. Indeed, Im(1/(1+iτ)) = −τ/(1+τ²), so the circle is left from z = 1 downward.
2. Since dz = −dw/w², the integral becomes −∫_{1−i∞}^{1+i∞} exp(πBw + πC/(k²w))·w^{−2} dw.
3. Substitute s = πBw. This gives −πB·∫_{πB−i∞}^{πB+i∞} e^{s + x²/(4s)}·s^{−2} ds, with x = 2π√(BC)/k.
4. By Fact 3.6 this equals −πB·(2πi)·(2/x)·I_1(x).
5. Multiplying by i/k² gives 4π²B·I_1(x)/(k²x) = (2π/k)√(B/C)·I_1(x). ∎

**Numerical confirmation.** The left side of Lemma 3.7 was evaluated at k = 10, n = 7 by oscillatory quadrature on Re w = 1. It gives 0.08174002496770177, against P_10 = 0.08174002496770177.

**Completing the arc.** The integrand of M_{h,k} has modulus e^{π/3}·|e^{πCz/k²}| on K∖{0}, which is bounded. So ∮_K converges absolutely. Going around K clockwise, the complement of the Rademacher arc z' → 1 → z'' is the arc z'' → 0 → z'. Hence

  M_{h,k} = ∮_K(…) − ∫_{z''→0}(…) − ∫_{0→z'}(…).

On the two removed arcs, Lemma 3.3(d) and Corollary 3.4 give

  |integrand| ≤ e^{π/3}·e = e^{1+π/3}.

The two arcs have total length at most π√2·k/(N+1). Therefore

  (i/k²)·M_{h,k} = P_k(n) + ρ_{h,k},  with |ρ_{h,k}| ≤ (1/k²)·π√2·(k/(N+1))·e^{1+π/3} = π√2·e^{1+π/3}/(k(N+1)).   (3.5)

**Amplitudes.** Define

  a_k(n) = Σ_{0≤h<k, (h,k)=1} Φ_{h,k}·e(−nh/k)  for 5 | k.

It is a sum of φ(k) unimodular terms, so **|a_k(n)| ≤ φ(k)**. For k = 5 it equals Σ_{h=1}^{4} w_h ζ^{−nh} = a(n). The main terms of all growing cusps with k ≤ N therefore add up to

  a(n)P_5(n) + Σ_{5|k, 10≤k≤N} a_k(n)P_k(n).

---

## 3.6 Growing cusps: the remainder R_{h,k}

The integrand of R_{h,k} is holomorphic on Re z > 0, and the region between the Rademacher arc and the chord s_{h,k} lies in Re z > 0. So by Cauchy's theorem R_{h,k} may be taken along the chord. There w = Re(1/z) ≥ 1 by Lemma 3.3(c). By (3.4),

  |e^{π/(3z)}(Ĝ − 1)| ≤ Σ_{j≥1} φ_j·e^{(π/3 − 2πj)w}.

Each exponent π/3 − 2πj is negative for j ≥ 1, so each term is decreasing in w. This is where the condition B = 1/3 < 2 is used; it is the analogue of the source's δ ≤ 4. The sum is therefore at most its value at w = 1:

  e^{π/3}·Σ_{j≥1} φ_j·t^j = e^{π/3}·(f(t)f(t⁵) − 1),  with t = e^{−2π}.

With Corollary 3.4 and the chord length from Lemma 3.3(c),

  |(i/k²)·R_{h,k}| ≤ (1/k²)·(2√2k/(N+1))·e^{1+π/3}·(f(t)f(t⁵) − 1).   (3.6)

**The correction term is not extracted.** The first correction term, j = 1, grows like e^{(π/3 − 2π)/z}, which decays. So, unlike Q_8, where δ > 2, no correction term is extracted. At general p the correction grows only when δ(p−1)/12 > 2, that is when δ > 24/(p−1) = 6. The bracketed [·] terms of the source's constants are therefore absent.

---

## 3.7 Non-growing cusps

Here 5 ∤ k, and there are φ(k) values of h for each such k ≤ N. Move each integral in (3.1) onto the chord. The integrand G_5·e^{2πnz/k²} is holomorphic on Re z > 0. By (3.3), on the chord

  |G_5(x)·e^{2πnz/k²}| = √5·e^{−πw/15}·|(x';x')_∞|·|F(x_5'')|·|e^{πCz/k²}|.

The factors are bounded as follows, using w ≥ 1:

- |(x';x')_∞| ≤ ∏(1 + e^{−2πwm}) ≤ f(e^{−2πw}) ≤ f(e^{−2π});
- |F(x_5'')| ≤ f(|x_5''|) = f(e^{−2πw/5}) ≤ f(e^{−2π/5});
- e^{−πw/15} ≤ e^{−π/15}.

Each bound is a non-increasing function of w, which justifies evaluating at w = 1. Together with Corollary 3.4,

  |integrand| ≤ √5·e^{−π/15}·f(e^{−2π/5})·f(e^{−2π})·e = K_5·e.

Hence each non-growing cusp contributes at most (1/k²)·(2√2k/(N+1))·e·K_5.   (3.7)

---

## 3.8 Assembly and the three constants

Collecting (3.1) and (3.5)–(3.7):

  s(n) − a(n)P_5(n) − Σ_{5|k, 10≤k≤N} a_k(n)P_k(n) = E(n),

where

  |E(n)| ≤ Σ_{5|k≤N} φ(k)·[π√2·e^{1+π/3} + 2√2·e^{1+π/3}·(f(t)f(t⁵) − 1)]/(k(N+1)) + Σ_{5∤k≤N} φ(k)·2√2·e·K_5/(k(N+1)).

**Counting.**

- For 5 | k we have φ(k)/k ≤ 4/5, and there are ⌊N/5⌋ ≤ N/5 multiples of 5 up to N. So Σ_{5|k≤N} φ(k)/(k(N+1)) ≤ (4/25)·N/(N+1) < 4/25 = (p−1)/p².
- For the non-growing cusps, φ(k)/k ≤ 1, and there are at most N values of k. So the sum is at most N/(N+1) < 1.

This bound is wasteful, because it ignores that 5 ∤ k and that φ(k) < k. It is kept to match the source's E_ng, and it is still valid.

This gives

- |E(n)| ≤ E_arc + E_rem + E_ng, where
- E_arc = √2π·(4/25)·e^{1+π/3} (from the arc completions, (3.5)),
- E_rem = 2√2·(4/25)·e^{1+π/3}·(f(t)f(t⁵) − 1) (from the remainders, (3.6)),
- E_ng = 2√2·e·K_5 (from the non-growing cusps, (3.7)).

Finally, |a_k(n)| ≤ φ(k) and P_k(n) > 0, and N = K(n). Hence

  |s(n) − a(n)P_5(n)| ≤ Σ_{5|k, 10≤k≤K(n)} φ(k)P_k(n) + E_arc + E_rem + E_ng.

This proves Theorem 3.1, modulo the numerical values in §3.9. ∎

---

## 3.9 Numerical values

The values below were computed with mpmath at 60 significant digits (f(x) = 1/qp(x)). Each "rounded up" value is at least the true value.

| constant | value (20 digits) | rounded up |
|---|---|---|
| K_5 = √5·e^{−π/15}·f(e^{−2π/5})·f(e^{−2π}) | 2.8549954049908738961 | 2.8549955 |
| E_arc = √2π·(4/25)·e^{1+π/3} | 5.5064468677889376771 | 5.506447 |
| E_rem = 2√2·(4/25)·e^{1+π/3}·(f(t)f(t⁵) − 1) | 0.0065708632316640301 | 0.006571 |
| E_ng = 2√2·e·K_5 | 21.950523842235242584 | 21.950524 |
| E_con | 27.463541573255844291 | **27.4636** |

The value "K_5 = 2.8549954" quoted in PROOF.md is a truncation, 4.99·10^{−8} below the true value. It is harmless, because E_ng is computed from the exact K_5, but the rounded-up value is 2.8549955.

The following were also checked at 60 digits:

- s(1,5) = 1/5, s(2,5) = s(3,5) = 0 and s(4,5) = −1/5 (the middle two to 10^{−61});
- the five values of a(n) against the closed form in §3.1;
- the transformation formulas (3.2) and (3.3) (§3.4);
- Lemma 3.7 (§3.5).

As an end-to-end sanity check, the bound was evaluated for every n ≤ 400 from exact coefficients. The largest ratio of left side to right side was 0.766. The earlier ball-arithmetic survey in PROOF.md covers n ≤ 200 000.

---

## 3.10 Status of the steps

- **Complete.** Every step above is either proved here or is one of the four cited standard facts:
  - Fact 3.5, the η / partition-function transformation (Apostol Thm 5.1);
  - Lemma 3.2, the Rademacher/Ford path (Apostol Ch. 5; Rademacher Ch. 14);
  - the Farey-neighbour inequality k + k' ≥ N + 1 (Hardy–Wright Thm 30);
  - Fact 3.6, the Bessel contour integral (DLMF 10.9.19 / Watson §6.22).
- **No hypothesis fails at p = 5, δ = 1:**
  - there is no branch issue, since the square roots cancel exactly;
  - no correction term is extracted, since B = 1/3 < 2;
  - there is no large-n hypothesis, since only n ≥ 1 is used, for C > 0;
  - the Farey order N = K(n) satisfies both (N1) and (N2).
- **A deviation from the brief.** The Farey order must be N = K(n), not K(n) − 1. With K(n) − 1 the cusp k = 5 is missing at n = 1, since the order is then 4.
- **Deliberately wasteful but valid bounds.**
  - E_rem uses f(t) in place of (−t;t)_∞.
  - E_ng sums over all k ≤ N with φ(k) ≤ k.
  - The proof actually gives the error with factor N/(N+1) < 1.
  - Sharper constants are available but not needed.
