# The Third Borwein Conjecture (p = 5): a computer-assisted proof

*Integrated write-up, 19 September 2026. This document replaces the working file `PROOF.md` (§§1–16) and the gap and audit reports in `/home/claude/third-borwein`. It is meant to be read on its own. Wherever a number comes from a computation, the script and log that produce it are named. Appendix A lists every certificate together with its command and the output it must print.*

---

## Summary

Let $S_n(q)=(q;q)_{5n}/(q^5;q^5)_n=\prod_{k\le 5n,\,5\nmid k}(1-q^k)$. The Third Borwein conjecture says that the coefficient $c_n(m)=[q^m]S_n$ is $\ge 0$ when $5\mid m$ and $\le 0$ otherwise.

We give a **computer-assisted proof**. Every certificate is in Arb ball arithmetic or exact integer arithmetic. Piece (F) of the off-peak lemma (§7.2) was first certified in float64 interval arithmetic and has since been re-run in full in Arb, with identical margins (minimum 47.91; `arb_F_summary_full.log`).

By palindromy it suffices to treat the half-range $0\le m\le 5n^2$. The proof has an exact part and three analytic regions.

1. **Exact computation.** Every coefficient in every class, for $n\le 979$ (§4).
2. **The band** $t=m/N<4$, with $N=5n+1$ (§5). The finite range is covered exactly and the infinite range by an analytic majorant.
3. **The Laplace region** $t\ge 4$, $\mu=m/(10n^2)<0.105$ (§6). Here a saddle-point analysis is applied to an exact Bessel–Laplace representation. It is certified by Arb tilings in one variable $V\in[2,12]$ plus an analytic tail cell for $V\ge 11.99$.
4. **The centre** $\mu\ge 0.105$ (§7). This part uses a near-peak Laplace analysis together with an off-peak lemma bounding $|S_n|$ on the rest of the circle.

The three regions of items 2–4 partition the half-range for each $n$ (§8); item 1 covers the small $n$.

The only analytic input beyond elementary lemmas proved here is an explicit circle-method expansion of $G_5=(q;q)_\infty/(q^5;q^5)_\infty$ with error constant $E_{\rm con}\le 27.4636$ (§3). It is proved in full from four textbook facts and was refereed independently.

**Status.** Every region is certified, and no item is pending. The audit record is in Appendix B:

- The certificates for the band (classes 3–4), the Laplace region (classes 3–4), the centre near-peak and combination step, and the off-peak lemma have each passed two independent adversarial audits.
- The classes 0–2 certificates (band and Laplace) have two independent adversarial audits (GAP_AUDIT012.md, AUDIT2_CLASS012.md), the second with its own exact method (CRT over 82 primes) at 76 hard points: all signs correct, true error 15–260× below the certified bound.
- §3 was refereed independently.
- An outside referee read the whole document end to end. They found one coverage gap, a float-drift sliver at $V=12$ in the classes 0–2 Laplace tiling, and it is now closed (§6.3). Their other points were presentational and are fixed in this version.

---

## Notation

| symbol | meaning |
|---|---|
| $S_n(q)$, $c_n(m)=c(m)$ | the polynomial above and its coefficients; $\deg S_n = 10n^2$ |
| $G_5(q)=\sum_i s(i)q^i$ | $(q;q)_\infty/(q^5;q^5)_\infty$ |
| $E_n(q)=\sum_j e(j)q^j$ | $\prod_{k>5n,\,5\nmid k}(1-q^k)^{-1}$ |
| $N$ | $5n+1$ (in §3 only, $N$ is the Farey order instead; see there) |
| $t,\ \mu$ | $t=m/N$, $\mu=m/(10n^2)$ |
| $a_5$ | $2\pi^2/75$ |
| $\zeta$ | $e^{2\pi i/5}$, except inside §3, where $\zeta$ is also $e(1/5)$ |
| $a(i)$ | $2\cos(2\pi i/5+\pi/5)+2\cos(4\pi i/5)$, the $k=5$ amplitude |
| $P_k(i)$ | $(2\pi/k)\sqrt{B/C}\,I_1\!\big((2\pi/k)\sqrt{BC}\big)$, with $B=1/3$ and $C=2i-1/3$ |
| $\mathrm{Li}_1(x)$ | $-\log(1-x)$ |
| $B_1(x)=x-\tfrac12$, $B_2(x)=x^2-x+\tfrac16$ | Bernoulli polynomials; $\tilde B_k$ are their 1-periodic extensions |
| "class" | $m \bmod 5$ |

The labels used throughout:

- **[proved]**: a complete argument is written out here.
- **[certified]**: an inequality verified in ball arithmetic (Arb, via python-flint) over its whole parameter range.
- **[exact]**: verified by exact integer arithmetic.
- **[float-interval]**: verified in float64 interval arithmetic with directed-rounding margins, but not in Arb.
- **[standard]**: a textbook fact, listed in Appendix C.

---

## 1. Statement and reduction

**Conjecture (Andrews; P. Borwein).** For every $n\ge1$ and every $m$,
$$c_n(m)\ \ge 0 \text{ if } m\equiv 0 \pmod 5,\qquad c_n(m)\ \le 0 \text{ if } m\not\equiv 0\pmod 5 .$$

**Lemma 1.1 (palindromy). [proved]** $c_n(10n^2-m)=c_n(m)$. The map $m\mapsto 10n^2-m$ preserves class 0 and swaps classes $1\leftrightarrow4$ and $2\leftrightarrow3$. Hence it suffices to prove the sign pattern for $0\le m\le 5n^2$, that is, for $0\le\mu\le\tfrac12$.

*Proof.* $S_n$ is a product of $4n$ factors $1-q^k$. Each factor satisfies $q^k(1-q^{-k})=-(1-q^k)$. Therefore
$$q^{D}S_n(1/q)=(-1)^{4n}S_n(q)=S_n(q),\qquad D=\sum_{k\le 5n,\,5\nmid k}k=\tfrac{5n(5n+1)}2-\tfrac{5n(n+1)}2=10n^2 .$$
Since $D\equiv0\pmod 5$, we have $10n^2-m\equiv -m\pmod 5$, and the sign pattern is invariant under $s\mapsto -s$ on classes. ∎

**Scaled variables.** Write $N=5n+1$, $t=m/N$ and $\mu=m/(10n^2)$. The three regions are:

- **band:** $t<4$;
- **Laplace region:** $t\ge 4$ and $\mu<0.105$;
- **centre:** $\mu\ge 0.105$.

Their union is $\{0\le m\le 5n^2\}$ trivially. §8 checks that each region's certificate really covers the whole of its region for $n\ge 980$.

---

## 2. Factorisation and the vanishing classes

**Lemma 2.1 (factorisation). [proved]** $S_n=G_5\cdot E_n$. The series $E_n$ has non-negative coefficients, with $e(0)=1$ and $e(j)=0$ for $0<j<N$. Consequently
$$c_n(m)=\sum_{j\ge0}e(j)\,s(m-j)=s(m)+\sum_{j\ge N}e(j)s(m-j).$$

*Proof.* We have $(q;q)_\infty/(q^5;q^5)_\infty=\prod_{5\nmid k}(1-q^k)$. Split this product at $k=5n$. The factor over $k>5n$ is $E_n^{-1}$, and $E_n$ is a product of geometric series in $q^k$ with $k\ge 5n+1=N$. ∎

In particular $e(j)$ counts the partitions of $j$ into parts $>5n$ that are not divisible by 5. For $j<4N$ at most three parts fit.

**Lemma 2.2 (vanishing classes). [proved]** $s(i)=0$ whenever $i\equiv 3,4\pmod 5$.

*Proof.* By Euler's pentagonal theorem, $(q;q)_\infty=\sum_{k\in\mathbb Z}(-1)^kq^{k(3k-1)/2}$. Modulo 5, the exponent $k(3k-1)/2$ takes the values $0,1,0,2,2$ for $k\equiv0,1,2,3,4$. So $(q;q)_\infty$ is supported on the residues $\{0,1,2\}$. Multiplying by $1/(q^5;q^5)_\infty$, a series in $q^5$, preserves this. ∎

**Consequence.** In classes 3 and 4 the main term of $c(m)$ cancels, and the sign comes from the $E_n$-convolution; these are the "thin" classes. In classes 0, 1, 2, $s(m)$ has a non-vanishing main term $a(m)P_5(m)$ (§3), where
$$a(m)=\tfrac{5+\sqrt5}2,\ -\sqrt5,\ -\tfrac{5-\sqrt5}2,\ 0,\ 0\qquad(m\equiv0,1,2,3,4).$$

**Remark (cusp by cusp). [certified numerically]** In the expansion of §3, the $k=10$ and $k=15$ amplitudes also vanish exactly on $i\equiv3,4\pmod 5$, to 390 or more digits (`cusp10.py`, `cusp15b.py`; PROOF.md §2). This is not used in the proof.

---
## 3. The expansion of $G_5$

*Source: `SEC3_PROOF.md`, refereed in `REFEREE_SEC3.md` (verdict: correct, no error, no missing term). The referee listed three presentational gaps, (a) to (c), and all three are filled below; they are marked **[gap (a)]**, **[gap (b)]** and **[gap (c)]**.*

**Notation for this section only.** The coefficient index is written $\nu$, and $i=\sqrt{-1}$. $N$ is the Farey order, not $5n+1$. The remaining symbols:

- $e(x)=e^{2\pi ix}$;
- $f(x)=F(x)=1/(x;x)_\infty$;
- $s(h,k)=\sum_{r=1}^{k-1}\frac rk\big(\big(\frac{hr}k\big)\big)$ is the Dedekind sum, and $\omega(h,k)=e^{\pi i s(h,k)}$.

### 3.1 Statement

Put $B=1/3$, $C=C(\nu)=2\nu-1/3$, $K(\nu)=\lfloor\sqrt{4\pi\nu}\rfloor+2$, $P_k(\nu)=\frac{2\pi}k\sqrt{B/C}\,I_1\big(\frac{2\pi}k\sqrt{BC}\big)$ and $a(\nu)=\sum_{h=1}^4w_h\zeta^{-\nu h}$, where $w_h=e^{-\pi is(h,5)}$ and $\zeta=e(1/5)$.

**Theorem 3.1. [proved]** For every $\nu\ge1$,
$$\big|s(\nu)-a(\nu)P_5(\nu)\big|\le\sum_{5\mid k,\ 10\le k\le K(\nu)}\varphi(k)P_k(\nu)+E_{\rm arc}+E_{\rm rem}+E_{\rm ng},$$
where

- $E_{\rm arc}=\sqrt2\pi\cdot\frac4{25}e^{1+\pi/3}\le 5.506447$;
- $E_{\rm rem}=2\sqrt2\cdot\frac4{25}e^{1+\pi/3}\big(f(\tau)f(\tau^5)-1\big)\le0.006571$, with $\tau=e^{-2\pi}$;
- $E_{\rm ng}=2\sqrt2\,e\,K_5\le 21.950524$, with $K_5=\sqrt5e^{-\pi/15}f(e^{-2\pi/5})f(e^{-2\pi})\le2.8549955$.

Their sum is $E_{\rm con}\le 27.463542\le\mathbf{27.4636}$.

Every $P_k(\nu)$ is positive: $C\ge5/3$, and $I_1>0$ on $(0,\infty)$. §3.9 shows $s(1,5)=1/5$, $s(2,5)=s(3,5)=0$ and $s(4,5)=-1/5$. So $w_1=e^{-\pi i/5}$, $w_2=w_3=1$, $w_4=e^{\pi i/5}$, and $a(\nu)=2\cos(2\pi\nu/5+\pi/5)+2\cos(4\pi\nu/5)$, as in §2.

### 3.2 Farey order

Take $N=K(\nu)$. Two properties are used.

- **(N1)** $N\ge5$ for all $\nu\ge1$, since $\sqrt{4\pi}>3$. So the denominator $k=5$ always occurs. (The choice $N=K(\nu)-1$ fails at $\nu=1$: it gives $N=4$.)
- **(N2)** $N+1=\lfloor\sqrt{4\pi\nu}\rfloor+3>\sqrt{4\pi\nu}$. Hence $2\pi C/(N+1)^2=(4\pi\nu-2\pi/3)/(N+1)^2<1$.

### 3.3 The Rademacher path

**Lemma 3.2 (Farey dissection) [standard, Fact C.2].** Take the Farey fractions of order $N$ cyclically on $[0,1)$, and let $h_1/k_1<h/k<h_2/k_2$ be neighbours. Let $C(h/k)$ be the Ford circle with centre $h/k+i/(2k^2)$ and radius $1/(2k^2)$, and let $\gamma_{h,k}$ be its upper arc between its tangency points with $C(h_1/k_1)$ and $C(h_2/k_2)$. Then:

- these arcs join into a path from $\tau_0$ to $\tau_0+1$ in $\mathbb H$;
- $s(\nu)=\int_{\tau_0}^{\tau_0+1}G_5(e(\tau))e(-\nu\tau)\,d\tau$, taken over a horizontal segment, may be deformed onto this path by Cauchy's theorem. The integrand is holomorphic on $\mathbb H$ and 1-periodic.

The substitution $\tau=h/k+iz/k^2$ maps $C(h/k)$ onto $K:\ |z-\tfrac12|=\tfrac12$. Also $d\tau=(i/k^2)dz$ and $e(-\nu\tau)=e(-\nu h/k)e^{2\pi\nu z/k^2}$. This gives
$$s(\nu)=\sum_{k\le N}\ \sum_{\substack{0\le h<k\\(h,k)=1}}\frac i{k^2}e(-\nu h/k)\int_{z'}^{z''}G_5\big(e(h/k+iz/k^2)\big)e^{2\pi\nu z/k^2}dz, \tag{3.1}$$
with $z'=(k^2+ikk_1)/(k^2+k_1^2)$ and $z''=(k^2-ikk_2)/(k^2+k_2^2)$. The integral runs clockwise along $K$ from $z'$ ($\operatorname{Im}>0$) through $z=1$ to $z''$ ($\operatorname{Im}<0$). The tangency point of $C(h/k)$ with $C(h_1/k_1)$ is $h/k-k_1/(k(k^2+k_1^2))+i/(k^2+k_1^2)$, and it maps to $z'$; the same holds for $z''$.

**[gap (b)] The arc of $0/1$.** In the cyclic dissection of $[0,1)$, the circle $C(0/1)$ contributes two pieces:

- the arc from $\tau_0=i$ to its tangency with $C(1/N)$;
- the arc of $C(1/1)=C(0/1)+1$ from its tangency with $C((N-1)/N)$ to $\tau_0+1=i+1$.

The integrand $G_5(e(\tau))e(-\nu\tau)$ is 1-periodic. Translate the second piece by $-1$: it becomes the arc of $C(0/1)$ from the tangency with $C(-1/N)$ to $i$. The two pieces then join into a single arc of $C(0/1)$, from the tangency with $C(-1/N)$ through $i$ to the tangency with $C(1/N)$.

In the $z$-variable ($h=0$, $k=1$, $\tau=iz$), this is the arc of $K$ from $z'$ (with $k_1=N$) through $z=1$ (that is, $\tau=i$) to $z''$ (with $k_2=N$). So the term $h/k=0/1$ in (3.1) is one Rademacher arc like every other, and Lemma 3.3 applies to it with $k_1=k_2=N$. (These neighbours satisfy $k+k_j=N+1$.)

**Lemma 3.3 (geometry). [proved]**

- **(b)** $|z'|,|z''|\le\sqrt2k/(N+1)$, and $\operatorname{Re}z',\operatorname{Re}z''\le2k^2/(N+1)^2$.
- **(c)** On the chord $s_{h,k}=[z',z'']$: $\operatorname{Re}(1/z)\ge1$, and $0<\operatorname{Re}z\le 2k^2/(N+1)^2$. Also $|s_{h,k}|\le2\sqrt2k/(N+1)$.
- **(d)** On the arcs of $K$ from $z''$ to $0$ and from $0$ to $z'$ (the complement of the Rademacher arc): $\operatorname{Re}(1/z)=1$ and $\operatorname{Re}z\le2k^2/(N+1)^2$. Each of these arcs has length at most $(\pi/2)\sqrt2k/(N+1)$.

*Proof.*

(b) Farey neighbours satisfy $k+k_j\ge N+1$ [standard, Fact C.3]. So $k^2+k_j^2\ge(N+1)^2/2$. Therefore $|z'|=k/\sqrt{k^2+k_1^2}\le\sqrt2k/(N+1)$ and $\operatorname{Re}z'=k^2/(k^2+k_1^2)\le2k^2/(N+1)^2$. The same holds for $z''$.

(c) The map $z\mapsto1/z$ sends $K\setminus\{0\}$ onto the line $\operatorname{Re}w=1$, and the open disc onto $\operatorname{Re}w>1$. The chord lies in the closed disc and avoids 0. $\operatorname{Re}z$ is affine on the chord and positive at both ends. Finally, $|s_{h,k}|\le|z'|+|z''|$.

(d) $\operatorname{Re}(1/z)=1$ on $K\setminus\{0\}$. Let $\theta\in(0,\pi]$ be the central angle of the arc from $z'$ to 0. Its chord has length $|z'|=\sin(\theta/2)$ and its arc length is $\theta/2\le(\pi/2)\sin(\theta/2)$, by Jordan's inequality. On a semicircle, $\operatorname{Re}z=(1+\cos\psi)/2$ is monotone, so $\operatorname{Re}z\le\operatorname{Re}z'$ on this arc. The arc from $z''$ is handled the same way. ∎

**Corollary 3.4.** On every chord and every arc in (d), $|e^{\pi Cz/k^2}|\le e^{2\pi C/(N+1)^2}<e$. This follows from Lemma 3.3 and (N2).

**Remark (the regions of deformation). [gap (c)]** Two facts about the region $\Delta_{h,k}$ bounded by the Rademacher arc $z'\to1\to z''$ and the chord $s_{h,k}$ are used below.

- **$\Delta_{h,k}$ lies in the half-plane $\operatorname{Re}z\ge\min(\operatorname{Re}z',\operatorname{Re}z'')>0$.** On each of the two semicircular pieces of the arc, $\operatorname{Re}z$ is monotone, with its maximum at $z=1$. So on the arc, $\operatorname{Re}z\ge\min(\operatorname{Re}z',\operatorname{Re}z'')$. On the chord this holds by affinity. The region $\Delta_{h,k}$ is a circular segment, which is the convex hull of the arc, so it lies in that half-plane as well.
- **The integrands are holomorphic on $\Delta_{h,k}$.** In particular $\Delta_{h,k}$ avoids $z=0$. For $\operatorname{Re}z>0$, $\tau=h/k+iz/k^2$ has $\operatorname{Im}\tau=\operatorname{Re}z/k^2>0$, so $G_5(e(\tau))$ and the functions derived from it in §3.4 are holomorphic there.

Cauchy's theorem therefore allows every integral over a Rademacher arc to be replaced by the integral over its chord. This is used in §3.6 and §3.7 alike.

### 3.4 The modular transformation at a cusp

**Fact 3.5 [standard, Fact C.1].** Let $k\ge1$, $(h,k)=1$, $\operatorname{Re}z>0$, and $hH\equiv-1\pmod k$. Put $x=\exp(2\pi ih/k-2\pi z/k^2)$ and $x'=\exp(2\pi iH/k-2\pi/z)$. Then
$$F(x)=\omega(h,k)(z/k)^{1/2}\exp\Big(\frac\pi{12z}-\frac{\pi z}{12k^2}\Big)F(x'),$$
with the principal square root. Equivalently, $(x;x)_\infty=\omega(h,k)^{-1}(z/k)^{-1/2}\exp(-\frac\pi{12z}+\frac{\pi z}{12k^2})(x';x')_\infty$.

Let $x=e(\tau)$ with $\tau=h/k+iz/k^2$. Then $G_5(x)=(x;x)_\infty F(x^5)$.

**Case $5\mid k$ (growing cusps), $k=5k_1$.** Write $x^5=\exp(2\pi ih/k_1-2\pi z_1/k_1^2)$ with $z_1=z/5$, and apply Fact 3.5 to $(h,k_1,z_1)$. Then
$$G_5(x)=\Phi_{h,k}\exp\Big(\frac\pi{3z}-\frac{\pi z}{3k^2}\Big)\hat G_{h,k}(z), \tag{3.2}$$
where

- $\Phi_{h,k}=\omega(h,k_1)/\omega(h,k)$, which is unimodular;
- $\hat G_{h,k}=(x';x')_\infty F(x_5')$, with $x_5'=\exp(2\pi iH_1/k_1-10\pi/z)$.

The square roots cancel exactly, because $(z_1/k_1)^{1/2}=(z/k)^{1/2}$ is the same principal root. So no branch choice enters. For the exponents: $-\frac\pi{12z}+\frac{5\pi}{12z}=\frac\pi{3z}$ and $\frac{\pi z}{12k^2}-\frac{5\pi z}{12k^2}=-\frac{\pi z}{3k^2}$. For $k=5$: $\Phi_{h,5}=w_h$.

**Case $5\nmid k$ (non-growing cusps).** Apply Fact 3.5 to $(5h,k,5z)$ with $5hH_5\equiv-1\pmod k$. Then
$$G_5(x)=\Phi'_{h,k}\sqrt5\exp\Big(-\frac\pi{15z}-\frac{\pi z}{3k^2}\Big)(x';x')_\infty F(x_5''),\qquad x_5''=\exp(2\pi iH_5/k-2\pi/(5z)), \tag{3.3}$$
with $|\Phi'_{h,k}|=1$. Here $(z/k)^{-1/2}(5z/k)^{1/2}=\sqrt5$, because $z$ and $5z$ have the same argument.

**The shift in $C$.** In (3.1), the factor $e^{-\pi z/(3k^2)}$ combines with $e^{2\pi\nu z/k^2}$ to give $e^{\pi Cz/k^2}$, where $C=2\nu-1/3$.

(3.2) and (3.3) were checked against direct $q$-series, phases included:

- `SEC3_PROOF.md` §3.4: 8 cusps at 60 digits, relative error $\le4\cdot10^{-60}$;
- the referee: 10 cusps at two values of $z$, relative error $\le4.6\cdot10^{-49}$.

**Majorant of $\hat G$.** Put $u=e^{-2\pi/z}$ and $w=\operatorname{Re}(1/z)$. Expanding $(x';x')_\infty=\prod(1-e(mH/k)u^m)$ and $F(x_5')=\prod(1-e(mH_1/k_1)u^{5m})^{-1}$ gives $\hat G=1+\sum_{j\ge1}g_ju^j$. Coefficientwise,
$$|g_j|\le[y^j]\prod(1+y^m)f(y^5)\le[y^j]f(y)f(y^5)=:\phi_j,$$
because distinct partitions are partitions. Hence
$$|\hat G-1|\le\sum_{j\ge1}\phi_je^{-2\pi jw}. \tag{3.4}$$

### 3.5 Growing cusps: main term

By (3.1) and (3.2), the cusp $h/k$ with $5\mid k$ contributes $\frac i{k^2}e(-\nu h/k)\Phi_{h,k}[M_{h,k}+R_{h,k}]$, where
$$M_{h,k}=\int_{z'}^{z''}e^{\pi/(3z)+\pi Cz/k^2}dz,\qquad R_{h,k}=\int_{z'}^{z''}e^{\pi/(3z)+\pi Cz/k^2}(\hat G_{h,k}-1)\,dz.$$

**Fact 3.6 (Bessel integral) [standard, Fact C.4].** For $x>0$ and $c>0$,
$$I_1(x)=\frac{x/2}{2\pi i}\int_{c-i\infty}^{c+i\infty}e^{s+x^2/(4s)}s^{-2}\,ds .$$

**[gap (a)] From the Hankel loop to the vertical line.** DLMF 10.9.19 gives $I_1(x)=\frac{x/2}{2\pi i}\int_{-\infty}^{(0+)}e^{s+x^2/(4s)}s^{-2}\,ds$, over a Hankel loop that starts and ends at $-\infty$ and encircles the negative real axis. For integer order the only singularity of the integrand is the essential singularity at $s=0$; the loop has no branch cut to respect.

Fix $R>c$ and close the region between the loop and the line $\operatorname{Re}s=c$ with the two arcs $\Gamma_R^\pm$ of $|s|=R$ in $\{\operatorname{Re}s\le c\}$, running from $c\pm i\sqrt{R^2-c^2}$ to the loop. The region enclosed contains no singularity, so by Cauchy's theorem the loop integral equals the line integral truncated at height $\sqrt{R^2-c^2}$, plus the integrals over $\Gamma_R^\pm$ (and minus the loop's own parts beyond $|s|=R$).

On $\Gamma_R^\pm$:

- $|e^s|\le e^c$;
- $|e^{x^2/(4s)}|\le e^{x^2/(4R)}\le e^{x^2/(4c)}$;
- $|s^{-2}|=R^{-2}$;
- the total length is at most $2\pi R$.

So the arcs contribute $O(R^{-1})\to0$. The loop's parts beyond $|s|=R$ vanish too, since $e^s$ decays along the negative axis. The line integral converges absolutely, because the integrand is $O(|s|^{-2})$ on it. Letting $R\to\infty$ proves Fact 3.6. The referee also confirmed the vertical-line form by direct quadrature.

**Lemma 3.7. [proved]** Let $K$ be oriented clockwise. Then $\frac i{k^2}\oint_Ke^{\pi B/z+\pi Cz/k^2}dz=P_k$.

*Proof.* Put $w=1/z=1+i\rho$, $\rho\in\mathbb R$. The clockwise traversal of $K\setminus\{0\}$ corresponds to $\rho$ running from $-\infty$ to $\infty$. With $dz=-dw/w^2$, the integral becomes $-\int_{1-i\infty}^{1+i\infty}e^{\pi Bw+\pi C/(k^2w)}w^{-2}dw$. Substituting $s=\pi Bw$ gives $-\pi B\int e^{s+x^2/(4s)}s^{-2}ds$, with $x=2\pi\sqrt{BC}/k$. By Fact 3.6 this equals $-\pi B\cdot 2\pi i\cdot(2/x)I_1(x)$. Multiplying by $i/k^2$ gives $4\pi^2BI_1(x)/(k^2x)=P_k$. ∎

Numerical check at $k=10$, $\nu=7$: 0.08174002496770177 on both sides. The referee's check agrees to 30 digits.

**Completing the arc.** On $K\setminus\{0\}$ the integrand of $M_{h,k}$ has modulus $e^{\pi/3}|e^{\pi Cz/k^2}|$, which is bounded, so $\oint_K$ converges absolutely. The complement of the Rademacher arc is the arc $z''\to0\to z'$. By Lemma 3.3(d) and Corollary 3.4, the integrand there is at most $e^{1+\pi/3}$, and the total length is at most $\pi\sqrt2k/(N+1)$. Hence
$$\frac i{k^2}M_{h,k}=P_k(\nu)+\rho_{h,k},\qquad|\rho_{h,k}|\le\frac{\pi\sqrt2e^{1+\pi/3}}{k(N+1)}. \tag{3.5}$$

**Amplitudes.** Put $a_k(\nu)=\sum_{(h,k)=1}\Phi_{h,k}e(-\nu h/k)$ for $5\mid k$. This is a sum of $\varphi(k)$ unimodular terms, so $|a_k|\le\varphi(k)$, and $a_5=a$.

### 3.6 Growing cusps: remainder

By the Remark in §3.3, $R_{h,k}$ may be taken along the chord, where $w\ge1$. By (3.4),
$$|e^{\pi/(3z)}(\hat G-1)|\le\sum_{j\ge1}\phi_je^{(\pi/3-2\pi j)w}.$$
Each exponent is negative, so each term decreases in $w$. This is where $B=1/3<2$ is used. The sum is therefore at most its value at $w=1$, namely $e^{\pi/3}(f(\tau)f(\tau^5)-1)$ with $\tau=e^{-2\pi}$. With the chord length and Corollary 3.4,
$$\Big|\frac i{k^2}R_{h,k}\Big|\le\frac1{k^2}\cdot\frac{2\sqrt2k}{N+1}e^{1+\pi/3}\big(f(\tau)f(\tau^5)-1\big). \tag{3.6}$$
The $j=1$ term decays, so, unlike for $Q_8$, no correction term needs to be extracted.

### 3.7 Non-growing cusps

Here $5\nmid k$. By the Remark in §3.3 (the region between arc and chord lies in $\operatorname{Re}z>0$, away from $z=0$), the integral may be moved to the chord. By (3.3), for $w\ge1$ on the chord:

- $|(x';x')_\infty|\le f(e^{-2\pi w})\le f(e^{-2\pi})$;
- $|F(x_5'')|\le f(e^{-2\pi w/5})\le f(e^{-2\pi/5})$;
- $e^{-\pi w/15}\le e^{-\pi/15}$.

Each bound is non-increasing in $w$. With Corollary 3.4, the integrand is at most $K_5e$, so each such cusp contributes at most
$$\frac1{k^2}\cdot\frac{2\sqrt2k}{N+1}\cdot eK_5 . \tag{3.7}$$

### 3.8 Assembly

By (3.1) and (3.5)–(3.7),
$$s(\nu)-a(\nu)P_5(\nu)-\sum_{5\mid k,\,10\le k\le N}a_kP_k=E(\nu),$$
where
$$|E(\nu)|\le\sum_{5\mid k\le N}\frac{\varphi(k)\big[\pi\sqrt2+2\sqrt2(f(\tau)f(\tau^5)-1)\big]e^{1+\pi/3}}{k(N+1)}+\sum_{5\nmid k\le N}\frac{\varphi(k)\,2\sqrt2eK_5}{k(N+1)}.$$

Two counts finish the bound.

- For $5\mid k$: $\varphi(k)/k\le4/5$, and there are at most $N/5$ such $k$. So the first weight sum is less than $4/25$.
- For $5\nmid k$: $\varphi(k)/k\le1$, and there are at most $N$ values of $k$. So the second weight sum is less than 1.

This gives $|E|\le E_{\rm arc}+E_{\rm rem}+E_{\rm ng}$. Using $|a_k|\le\varphi(k)$ and $N=K(\nu)$ proves Theorem 3.1. ∎

### 3.9 Constants

Computed at 60 digits (`SEC3_PROOF.md` §3.9) and independently at 50 digits (`REFEREE_SEC3.md` item 6):

| constant | value | rounded up |
|---|---|---|
| $K_5$ | 2.85499540499087… | 2.8549955 |
| $E_{\rm arc}$ | 5.50644686778894… | 5.506447 |
| $E_{\rm rem}$ | 0.00657086323166… | 0.006571 |
| $E_{\rm ng}$ | 21.9505238422352… | 21.950524 |
| $E_{\rm con}$ | 27.4635415732558… | **27.4636** |

The Dedekind sums $s(1,5)=\tfrac15$, $s(2,5)=s(3,5)=0$ and $s(4,5)=-\tfrac15$ were checked at 60 digits.

### 3.10 Numerical confirmation and tightness (not used in the proof)

- **Referee** (`REFEREE_SEC3.md` item 7). Exhaustive over $\nu\in[1,5000]$ and at 1509 further $\nu\le200000$: no violation, maximum LHS/RHS 0.9045085.
- **Earlier survey** (PROOF.md §3). Ball arithmetic, every $\nu\le60000$ and every 7th $\nu\le200000$: no failures.
- **Tightness.** The ratio tends to $(5+\sqrt5)/8=0.9045085$ in the class $\nu\equiv2\pmod 5$. There, $|a_{10}|=(5+\sqrt5)/2$ is set against $\varphi(10)=4$. So the $\sum\varphi(k)P_k$ term must not be weakened in any later use.

**Uses of Theorem 3.1.** Write $R(i)$ for its right-hand side. The certificates use:

- $|s(i)-a(i)P_5(i)|\le R(i)$ for every $i\ge1$;
- exact values of $s(i)$ for $i\le200000$ (`g5_coef.py` → `g5_coeffs.pkl`, cross-checked by an independent recomputation, `bandfix_coeffs.py`).

---
## 4. Exact computation

**Theorem 4.1. [exact]** The conjecture holds, for every coefficient and in every class, for all $n\le 979$.

The whole range $n\le979$ is computed by one script, `centre3_exact.py NA NB`, run in segments.

- **Method.** The script works modulo $x^{5NB^2+1}$ and builds $S_{NA}$ from $S_0=1$. It obtains each next $S_n$ by four exact shift-subtracts $S\leftarrow S-x^kS$, for $k=5n+1,\dots,5n+4$ with $5\nmid k$. Truncation is exact, because each factor feeds only higher degrees. It then checks every coefficient with $0\le m\le5n^2$, in every class. A strict wrong sign counts as a violation; zeros are allowed.
- **Segments.** Every segment is logged:

  | segment | log | result |
  |---|---|---|
  | 1–150 | `exact_1_150.log` | `DONE n=1..150 total violations 0, 4s` |
  | 150–291 | `exact_150_291.log` | `DONE n=150..291 total violations 0, 34s` |
  | 291–410 | `centre3_exact_291_410.log` | DONE, 0 violations (114 s) |
  | 410–580 | `centre3_exact_410_580.log` | DONE, 0 violations (419 s) |
  | 580–820 | `centre3_exact_580_820.log` | DONE, 0 violations (1357 s) |
  | 820–920 | `centre3_exact_820_920.log` | `DONE n=820..920 total violations 0, 2825s` |
  | 920–1000 | `centre3_exact_920_1000.log` | lines $n=920,\dots,979$ all `wrong_sign=0`, then `Killed` (out of memory) |

  The segments abut or overlap, so the exact range is $1\le n\le979$. (AUDIT_CENTRE.md item 6 checked the line counts of the segments from 291 on; the 920–1000 log has 60 lines `wrong_sign=0`, for $n=920,\dots,979$.) An earlier run of `verify_chunk.py` for $n\le290$, recorded in NOTES.md, is superseded by the two logged segments. Zeros occur only at $m\lesssim5n$, i.e. $\mu<2\cdot10^{-3}$, which is inside the band. There the conjecture's signs are non-strict.
- **Palindromy** (Lemma 1.1) extends each check from $m\le5n^2$ to all $m$.

**Further exact checks (used in §5).**

- `bandfix_exact.py`: every class and every $m<4N$, for every $n\in[291,2000]$. The method is $c_n=G_5E_n\bmod q^{4N}$ with $6E_n=6+6p_1+3(p_1^2+p_1(q^2))+(p_1^3+3p_1(q^2)p_1+2p_1(q^3))$ (the cycle index of $S_2$ and $S_3$), where $p_1=\sum_{N\le k<4N,\,5\nmid k}q^k$, computed with `fmpz_poly.mul_low`. The script was validated against the full product for $n\le14$, and it reports zero sign failures (`bandfix_exact.log`, GAP_BAND.md §1). AUDIT_BAND.md recomputed 16 values of $n$ in $[291,2000]$ by a direct product that does not use $E_n$ or the coefficient pickle, with identical results for all $m<4N$.
- The coefficients $s(i)$ for $i\le200000$ come from `g5_coef.py` (pentagonal series times $p(k)$ at $q^5$). They were recomputed independently by `bandfix_coeffs.py`, which prints `g5_coeffs.pkl agrees with independent recomputation for i <= 200000` (`bandfix_coeffs.log`).

**Remark.** The analytic certificates below hold from $n\ge821$ (Laplace classes 3–4 and the centre) or from $n\ge291$ (Laplace classes 0–2). The overlap with the exact range is therefore $821\le n\le979$. From here on, **$n\ge980$**.

---

## 5. The band, $t<4$ (i.e. $m<4N$)

Every part of $E_n$ is $\ge N$. So for $m<4N$ at most three parts contribute, and
$$c(m)=s(m)+b_1+b_2+b_3,$$
where $b_j$ collects the $j$-part terms. For $j\in[N,2N)$, $e(j)=1$ if $5\nmid j$ and $e(j)=0$ otherwise.

### 5.1 Classes 3 and 4

*(Source: GAP_BAND.md.)*

**One-part block.** Let $m\equiv3$ or $4\pmod 5$ and write $m=5M+r$. Then $s(m)=0$, and
$$b_1=\sum_{\substack{j\in[N,m]\\5\nmid j}}s(m-j)=\sum_{\substack{i\le m-N\\ i\not\equiv m}}s(i).$$
Since $s(i)=0$ for $i\equiv3,4$, only $i\equiv0,1,2$ survive. Also $m-N=5(M-n)+r-1$ with $r-1\in\{2,3\}$, so every block $\{5t,5t+1,5t+2\}$ with $t\le T:=M-n$ is included in full. Hence
$$b_1=G(T)=\sum_{t\le T}g(t),\qquad g(t)=s(5t)+s(5t+1)+s(5t+2).$$

**Lemma 5.1. [proved, given §3]** The main term satisfies $\mathrm{main}(g(t))\le-5P_5'(5t)$. Moreover $g(t)<-4P_5'(5t)<0$ for $t\ge200$.

*Proof.* The main term of $g(t)$ is $a(0)P_5(5t)+a(1)P_5(5t+1)+a(2)P_5(5t+2)$. Since $a(0)=|a(1)|+|a(2)|$, this equals $-|a(1)|\,[P_5(5t+1)-P_5(5t)]-|a(2)|\,[P_5(5t+2)-P_5(5t)]$. By the series of $I_1$, $P_5(i)=\sum_{r\ge1}a_5^r x^{r-1}/(r!\,(r-1)!)$ is a power series with positive coefficients in $x=i-\tfrac16$ (not in $i$ itself). So $P_5'$ is increasing for $x>0$, and the main term is at most $-(|a(1)|+2|a(2)|)P_5'(5t)=-5P_5'(5t)$.

By Theorem 3.1, $|g(t)-\mathrm{main}|\le3R_{\rm grp}(5t+2)$, where $R_{\rm grp}$ is the (increasing) right-hand side of Theorem 3.1. Lemma L1 of `bandfix_analytic.py` states $3R_{\rm grp}(X+2)\le0.018\,P_5'(X)$ for all $X\ge1000$. Its proof has four steps:

1. write $P_5'(X)=4a_5^2I_2(Y)/Y^2$ with $Y=2\sqrt{a_5(X-1/6)}$;
2. use $I_2\ge I_{5/2}$ (monotonicity in the order, Fact C.5) and the closed form of $I_{5/2}$, which give $I_2(Y)\ge e^Y(1-3/Y-2e^{-2Y})/\sqrt{2\pi Y}$;
3. use $I_1\le I_{1/2}$, $P_k\le P_{10}$ for $k\ge10$, and $\sum_{5\mid k\le K}\varphi(k)\le K^2/10+K/2$;
4. bound the ratio by $c_1X^{3/2}e^{-Y/2}+c_2Y^{5/2}e^{-Y}$. Both terms decrease for $X\ge1000$ (sign of the log-derivative), so the bound evaluated at $X=1000$ in Arb holds for all larger $X$. The log prints `L1: sup_{X>=1000} 3 Rgrp(X+2)/P5'(X) <= [0.017972 +/- 1.12e-7]` (`bandfix_analytic.log`).

This gives $g(t)\le-5P_5'+0.018P_5'<-4P_5'$ for $5t\ge1000$. ∎

**Lemma 5.2. [exact]** For $t\le39999$, $g(t)<0$ except at $g(1)=1$ and $g(5)=0$. Every partial sum satisfies $G(T)\le0$, with equality only at $T=1$ (`bandfix_coeffs.py`, `bandfix_coeffs.log`). The computation uses the pickle of exact $s(i)$, which the same script recomputes independently first.

**Corollary 5.3.** $G(T)\le0$ for every $T$, and $G(T)\le g(T)<0$ for $T\ge200$. Hence, for classes 3 and 4:

- for $m<N$: $c(m)=s(m)=0$;
- for $N\le m<2N$: $c(m)=G(T)\le0$.

Both hold for every $n$.

**$2N\le m<4N$.**

- *For $980\le n\le2000$: [exact]* by `bandfix_exact.py`, which covers $291\le n\le2000$. Over that whole run the log prints the largest ratio $|b_2+b_3|/|b_1|$ as `0.0001`, at $n=291$, $m=5819$; AUDIT_BAND.md's rerun gives the unrounded value $1.1302\cdot10^{-4}$.
- *For $n>2000$ (so $N\ge10001$): [proved, given §3].* Write $I=m-2N$. Then $5T\ge I+N-3$.
  - **Denominator.** $|b_1|\ge|g(T)|\ge4P_5'(I+N-3)$, by Corollary 5.3 and the monotonicity of $P_5'$.
  - **Numerator.** $|b_2|+|b_3|\le\sum_{d\le I}w(2N+d)|s(I-d)|$. Here the number of pairs is at most $d/2+1$, and the number of triples at $3N+e$ is at most $(e+3)^2/12+1/2$. Also $|s(i)|\le e^{2\sqrt{a_5i}}$: this is exact for $i\le200000$, and beyond that Theorem 3.1 gives the stronger constant $7.7\cdot10^{-5}$.
  - **Split.** Split the sum at $D=\lfloor3\sqrt{2N/a_5}\log N\rfloor$. The exponent $2\sqrt{a_5}(\sqrt I-\sqrt{I+N-19/6})$ increases in $I$, so it is largest at $I=2N$. This gives
    $$\frac{|b_2|+|b_3|}{|b_1|}\le\Psi(N)=C\,N^{9/4}(\log N)^2e^{-2\sqrt{a_5}(\sqrt3-\sqrt2)\sqrt N}.$$
  - **Conclusion.** $d\log\Psi/dN<0$ exactly when $(c_c/2)\sqrt N>9/4+2/\log N$, where $c_c=2\sqrt{a_5}(\sqrt3-\sqrt2)$ and $c_c/2=0.16304$ (AUDIT_BAND.md (b)). This first holds at $N=257$ and persists for larger $N$, so $\Psi$ is decreasing for $N\ge257$. In Arb, $\Psi(10001)\le0.19769$ (`bandfix_analytic.log`: `BAND (analytic): (|b2|+|b3|)/|g(T)| < 0.1977 for all N >= 10001`). Hence $c(m)\le b_1(1-0.1977)<0$.

**Conclusion 5.4.** For classes 3 and 4, every $n\ge980$ and every $m<4N$: $c(m)\le0$.

### 5.2 Classes 0, 1, 2

- **$980\le n\le2000$: [exact]**, by `bandfix_exact.py` (§4).
- **$n>2000$: [proved, given §3 and exact coefficients]** (GAP_CLASS012.md §A, audited in GAP_AUDIT012.md §6). Write $c(m)=s(m)+\sum_{j\ge N}e(j)s(m-j)$. The weights satisfy $e(j)\le1$ for $j<2N$ and $e(j)\le1+(N+1)+((N+3)^2/12+1/2)\le N^2\le m^2$ for $j<4N$. With $A(I)=\sum_{i\le I}|s(i)|$, this gives
  $$|c(m)-s(m)|\le A(m-N)+m^2A(m-2N)\le A(m-N_{\min})+m^2A(m-2N_{\min}),\qquad N_{\min}(m)=\max(10001,\lfloor m/4\rfloor+1).$$
  The second inequality uses that $A$ is non-decreasing and that $m<4N$ is equivalent to $N\ge\lfloor m/4\rfloor+1$. So one inequality per $m$ covers every $n>2000$ at once.
  - *Part A, $m\le200000$ [exact].* $\operatorname{sgn}s(m)=\operatorname{sgn}a(m)$ and $|s(m)|>A(m-N_{\min})+m^2A(m-2N_{\min})$ for every $m$ in classes 0, 1, 2, except $s(7)=0$. That coefficient has $m<N$, so $c(7)=0$, which is allowed. The worst ratio is $5.06\cdot10^{-10}$.
  - *Part B, $m>200000$ [certified, given §3].* The inputs are:
    - $|s(i)|\le e^{c\sqrt i}$ with $c=2\sqrt{a_5}$;
    - $\sum_{i\le I}e^{c\sqrt i}\le e^{c\sqrt I}(1+2\sqrt I/c)$;
    - $P_5(m)\ge h_0m^{-3/4}e^{c\sqrt m}$, via $I_1\ge I_{3/2}$;
    - $R(m)\le1.3\,m\,e^{c\sqrt m/2}+E_{\rm con}$.

    Together these give $[\text{corrections}+R(m)]/((5-\sqrt5)/2\cdot P_5(m))\le Q\le5.14\cdot10^{-20}$. Each term has the form $Cm^pe^{-\kappa\sqrt m}$ with $p\in\{1.25,3.25,1.75,0.75\}$, and is decreasing once $\sqrt m>2p/\kappa$ (`cls012_band.py`, `cls012_band.log`).

**Conclusion 5.5.** For classes 0, 1, 2, every $n\ge980$ and every $m<4N$, $c(m)$ has the conjectured sign.

---
## 6. The Laplace region, $t\ge4$, $\mu<0.105$

*(Sources: GAP_CLASS012.md §B (classes 0–2), audited in GAP_AUDIT012.md and AUDIT2_CLASS012.md; GAP_LAP34.md (classes 3–4), audited in AUDIT_LAP34.md and AUDIT2_LAP34.md.)*

Both certificates share one exact representation (§6.1) and one Euler–Maclaurin lemma (§6.2). They differ in the saddle, because classes 3 and 4 are thin. Neither certificate uses any withdrawn component: `B_uniform5.py`, the thin-factor constant $K_0$, the Cauchy constant $M$, the $t_{\rm eff}$ device, monotonicity in $t$, `mono_cert.py`, `zone2_cert.py` or `tail3_cell.py`.

### 6.1 Exact Bessel–Laplace representation [proved]

Write $a(i)=\sum_{h=1}^4A_h\zeta^{hi}$ with $A_1=e^{i\pi/5}$, $A_4=e^{-i\pi/5}$ and $A_2=A_3=1$. By Lemma 2.1 and Theorem 3.1,
$$c(m)=e(m)+\sum_hA_h\zeta^{hm}J_h+\rho,\qquad J_h=\sum_{j<m}e(j)\zeta^{-hj}P_5(m-j),\qquad|\rho|\le\sum_{j<m}e(j)R(m-j).$$
Here $e(m)=e(m)s(0)$ is the $j=m$ term, and $\zeta^{-hj}\zeta^{hm}=\zeta^{h(m-j)}$ reproduces $a(m-j)$.

For $x=i-1/6>0$ we have $P_5(i)=\sqrt{a_5/x}\,I_1(2\sqrt{a_5x})$, which is the inverse Laplace transform of $e^{a_5/u}-1$. Put $r_1(u)=e^{a_5/u}-1-a_5/u=O(|u|^{-2})$. Then
$$\frac1{2\pi i}\int_{(W)}r_1(u)e^{xu}du=\begin{cases}P_5-a_5,&x>0,\\0,&x<0.\end{cases}$$
Since $r_1=O(|u|^{-2})$, Fubini applies, and for every $W>0$, exactly,
$$J_h=a_5H_h+K_h,\qquad H_h=\sum_{j\le m-1}e(j)\zeta^{-hj},\qquad K_h=\frac1{2\pi}\int_{\mathbb R}r_1(W+iy)\,e^{(m-1/6)(W+iy)}E_n(\zeta^{-h}e^{-W-iy})\,dy.$$

### 6.2 Euler–Maclaurin lemma for the twisted $E_n$ [proved]

**Lemma 6.1.** Let $\theta\in(0,1)$, and let $f$ be such that $f,f'\to0$ at $\infty$ and $f''$ is integrable. Then
$$\sum_{t\ge n}f(t+\theta)=\int_n^\infty f-B_1(\theta)f(n)+R,\qquad|R|\le\tfrac1{12}\Big(|f'(n)|+\int_n^\infty|f''|\Big).$$

*Proof.* On each unit interval $[t,t+1]$, write $f(t+\theta)-\int_t^{t+1}f=\int_t^{t+1}k(x-t)f'(x)\,dx$, where $k(x)=x-H(x-\theta)$ and $H$ is the Heaviside function. The mean of $k$ is $B_1(\theta)$, and $k-B_1(\theta)=\tilde B_1(x-\theta)$. The mean part telescopes to $-B_1(\theta)f(n)$. Integrate the remaining part by parts once more, with kernel $P(x)=(\tilde B_2(x-\theta)-B_2(\theta))/2$. Its constant part contributes $(B_2(\theta)/2)f'(n)$, with $|B_2(c/5)|/2\le 11/300<1/12$. Its periodic part is bounded by $\sup|\tilde B_2|/2=1/12$ times $\int|f''|$. ∎

The audit (GAP_AUDIT012.md §3) notes that the constant part must be split off: without the split, $\sup|P|$ reaches $0.12>1/12$ at $\theta=2/5$.

**Application.** Write $k=5t+c$ with $c\in\{1,2,3,4\}$. Put $\Phi_c(s)=\log(1-\zeta^{-hc}e^{-s})$, $u=W\zeta_\eta$ with $\zeta_\eta=1+i\eta$, $\sigma=5nu$, $V=5nW$ and $\xi=e^{-V}$. Apply Lemma 6.1 to $f(t)=\Phi_c(5ut)$ with $\theta=c/5$, and sum over $c$. This gives
$$\log E_n(\zeta^{-h}e^{-u})=-\gamma(\sigma)/u+\Omega_h(\sigma)+R_h,$$
where

- $\gamma(z)=\big[\mathrm{Li}_2(e^{-z})-\mathrm{Li}_2(e^{-5z})/5\big]/5$;
- $Q(s)=\sum_c\Phi_c(s)=\log(1-e^{-5s})-\log(1-e^{-s})$, which is independent of $h$ because $\prod_c(1-\zeta^{-hc}w)=(1-w^5)/(1-w)$;
- $\Omega_h(\sigma)=\sum_cB_1(c/5)\Phi_c(\sigma)$;
- $|R_h|\le\bar R:=(5W/3)|\zeta_\eta|\xi\big[1/(1-\xi)+|\zeta_\eta|/(1-\xi)^2\big]$.

*Derivation of $\bar R$.* Here $f(t)=\Phi_c(5ut)$, so $f'=5u\,\Phi_c'$ and $f''=25u^2\Phi_c''$, with $|u|=W|\zeta_\eta|$ and $\operatorname{Re}(5ut)=5Wt$. From $\Phi_c'(s)=\zeta^{-hc}e^{-s}/(1-\zeta^{-hc}e^{-s})$, $\Phi_c''(s)=-\zeta^{-hc}e^{-s}/(1-\zeta^{-hc}e^{-s})^2$ and $|1-\zeta^{-hc}e^{-s}|\ge1-e^{-\operatorname{Re}s}$ we get:

- $|\Phi_c'(s)|\le e^{-\operatorname{Re}s}/(1-e^{-\operatorname{Re}s})$, so $|f'(n)|\le5W|\zeta_\eta|\,\xi/(1-\xi)$, since $\operatorname{Re}(5un)=V$ and $e^{-V}=\xi$.
- $|\Phi_c''(s)|\le e^{-\operatorname{Re}s}/(1-e^{-\operatorname{Re}s})^2$. For $t\ge n$ the denominator is at least $(1-\xi)^2$, so
  $$\int_n^\infty|f''|\,dt\le\frac{25W^2|\zeta_\eta|^2}{(1-\xi)^2}\int_n^\infty e^{-5Wt}dt=\frac{5W|\zeta_\eta|^2\xi}{(1-\xi)^2}.$$

Lemma 6.1 gives $\tfrac1{12}(|f'(n)|+\int|f''|)$ for each $c$. Summing over the four values of $c$ gives $\bar R$.

The pairs $c\leftrightarrow5-c$ have opposite $B_1$ and conjugate roots of unity, so $\operatorname{Re}\Omega_h(W)=0$, i.e. $|e^{\Omega_h(W)}|=1$. Also $|\Omega_h|\le\bar\omega:=0.8\,\mathrm{Li}_1(\xi)$, since $\sum_c|B_1(c/5)|=0.8$.

*Checks.* In `cls012_em_check.py` the bound exceeds the true remainder by a factor 5–20, over 12 cases. In `lap34_emchk.py`, at $n=821$ with complex $u$ and $\eta$ up to 2, the factor is also 5–20.

### 6.3 Classes 0, 1, 2

**Saddle.** Define $\alpha(V)=a_5-\gamma(V)-VQ(V)/5$ and choose $W$ by $(m-1/6)W^2=\alpha(V)$, where $V=5nW$. With $X=1/W$ the exponent becomes
$$\frac{a_5}u+(m-\tfrac16)u-\frac{\gamma(\sigma)}u=X\,\Xi_V(\eta),\qquad\Xi_V(\eta)=\frac{a_5-\gamma(V\zeta_\eta)}{\zeta_\eta}+\alpha(V)\zeta_\eta .$$
$\Xi_V$ depends on $V$ only. It satisfies $\Xi_V'(0)=0$ identically (this is the definition of $\alpha$), and $\Xi(-\eta)=\overline{\Xi(\eta)}$. Put $\Xi_2=-\Xi''(0)>0$ and
$$\mathrm{Amp}=\frac{W^{3/2}e^{X\Xi(0)}}{\sqrt{2\pi\Xi_2}}.$$
Then
$$c(m)=\mathrm{Amp}\Big[M+\sum_hA_h\zeta^{hm}e^{\Omega_h}\delta_h\Big]+E_2',\qquad M=\sum_hA_h\zeta^{h\cdot\mathrm{cls}}e^{\Omega_h(W)}\in\mathbb R .$$
As $V\to\infty$, $M\to a(\mathrm{cls})$. $M$ is the main term with the correct normalisation $|E_n(\zeta_5e^{-W})|$.

**Error terms.** All are relative to Amp, and each is either independent of $X$ or non-increasing in $X$. Set $\eta_0=cc/\sqrt{X\Xi_2}$ with $cc=4$, and $s=\eta\sqrt{X\Xi_2}$.

- **Window, $|\eta|\le\eta_0$.** Write $\psi=X(ic_3\eta^3+R_4)+\Delta\Omega+R_h$. The cubic term is a pure phase, so
  $$|e^\psi-1|\le\min(k_3s^3,2)e^r+e^r-1,\qquad r=k_4s^4+k_1s+\bar R,$$
  with $k_3=|c_3|/(\Xi_2^{3/2}\sqrt X)$, $k_4=\sup|\Xi''''|/(24\Xi_2^2X)$ and $k_1=\beta_1/\sqrt{X\Xi_2}$. Here $\beta_1=0.8V\xi/(1-\xi)$ bounds $|d\Omega/d\eta|$. The supremum of $|\Xi''''|$ is enclosed on sub-balls of width 0.01 covering $[0,\eta_0(X_{\min})]$. The Gaussian truncation contributes $\operatorname{erfc}(cc/\sqrt2)$.
- **Mid region, $\eta_0<|\eta|\le2$.** $D(\eta)=\operatorname{Re}(\Xi(0)-\Xi(\eta))\ge d_k\eta^2$ on 170 fixed $\eta$-tiles, using the better of the Taylor bound $\Xi_2/2-\sup|\Xi''''|\eta^2/24$ and a direct ball evaluation. The contribution is
  $$e^{2\bar\omega+\bar R}\sum_k\sqrt{\Xi_2/(2d_k)}\,\operatorname{erfc}\big(\max(cc\sqrt{d_k/\Xi_2},\ \eta_k\sqrt{X_{\min}d_k})\big).$$
- **Crude terms $E_2$.** These are:
  - the $(1+a_5/u)$ part of $r_1$;
  - $|y|\in[2W,1]$, with $|r_1|\le e^{a_5X/(1+\eta_1^2)}+1+a_5X/\eta_1$;
  - $|y|>1$, with $|r_1|\le(a_5^2/2y^2)e^{a_5}$;
  - $e(m)\le e^{mW}E_n(e^{-W})$;
  - $a_5|H_h|\le a_5e^{(m-1)W}E_n(e^{-W})$;
  - the $\rho$-term, by Rankin's trick at the same $W$, with $P_k(i)\le W(e^{a_kX}-1)e^{(i-1/6)W}$ and $k\le K(m)$.

  *The $\rho$-term in detail.* By §6.1, $|\rho|\le\sum_{j<m}e(j)R(m-j)$, where $R(i)=\sum_{5\mid k,\,10\le k\le K(i)}\varphi(k)P_k(i)+E_{\rm con}$. Every $i=m-j$ here is $\ge1$, so Theorem 3.1 applies. Write $x=i-\tfrac16$ and $a_k=2\pi^2/(3k^2)$, so that $a_{10}=a_5/4$ and $a_k\le a_5/9$ for $k\ge15$.
  - *Bessel part.* The series of $I_1$ gives $P_k(i)=\sum_{r\ge1}a_k^rx^{r-1}/(r!(r-1)!)$. With $x^{r-1}/(r-1)!\le e^{xW}/W^{r-1}$ this gives $P_k(i)\le W(e^{a_kX}-1)e^{(i-1/6)W}$.
  - *Counting.* $\varphi(10)=4$, and $\sum_{j=3}^{K/5}\varphi(5j)\le\sum4j\le2K_5(K_5+1)$, where $K_5=K/5$ and $K$ bounds $K(i)$ for every $i\le m$. The code uses $K=\sqrt{4\pi(a_5X^2+1)}+2$, since $m-\tfrac16=\alpha X^2\le a_5X^2$.
  - *Rankin.* Since $e(j)\ge0$, $\sum_{j<m}e(j)e^{(m-j-1/6)W}\le e^{(m-1/6)W}E_n(e^{-W})$. So the Bessel part contributes at most $W\big[4e^{a_5X/4}+2K_5(K_5+1)e^{a_5X/9}\big]e^{(m-1/6)W}E_n(e^{-W})$.
  - *$E_{\rm con}$ part.* $e^{(m-j)W}\ge1$ for $j\le m$, so $\sum_{j<m}e(j)E_{\rm con}\le E_{\rm con}e^{mW}E_n(e^{-W})=E_{\rm con}e^{W/6}\cdot e^{(m-1/6)W}E_n(e^{-W})$.

  Together these give the script's term `Trho` $=F_m e^{W/6}\big(W[4e^{a_5X/4}+e^{a_5X/9}2K_5(K_5+1)]+E_{\rm con}\big)$, where $F_m$ is an upper bound for $e^{(m-1/6)W}E_n(e^{-W})/\mathrm{Amp}$. The script also applies the factor $e^{W/6}$ to the Bessel part, where it is not needed; this only enlarges the bound. (AUDIT2_LAP34.md §2 checked these constants against Theorem 3.1.)

  All the crude terms use $E_n(e^{-W})\le\exp(0.8X\mathrm{Li}_2(\xi)+4\mathrm{Li}_1(\xi))$, which holds because $E_n$ has non-negative coefficients, so $|E_n(\zeta^{-h}e^{-u})|\le E_n(e^{-W})$. The total is at most $\sqrt{2\pi\Xi_2}X^{3/2}e^{-cX}\mathrm{poly}(X)$, with $c\ge\kappa-a_5/4$ and $\kappa=a_5-\gamma-0.8\mathrm{Li}_2(\xi)$. The rate $a_5/4$ comes from the $k=10$ Rankin term; the original script asserted $a_5/5$, and the audit corrected it. The script asserts $c>0$ and $X_{\min}>4.5/c$, and $E_2\le2.3\cdot10^{-10}$.

**Certified condition.** $\operatorname{sgn}(a)M-4\delta-E_2>0$, where $\delta=\varepsilon_{\rm win}+\operatorname{erfc}(cc/\sqrt2)+\varepsilon_{\rm mid}$. Each factor $A_h\zeta^{hm}e^{\Omega_h(W)}$ is unimodular, because $\operatorname{Re}\Omega_h(W)=0$, and $|\delta_h|\le\delta$. So the bracket is $M$ plus at most $4\delta$, and the condition implies that $c(m)$ has the sign of $a(\mathrm{cls})$.

**Range of $X$.** $X=5n/V\ge1455/V$ for $n\ge291$. $t\ge4$ gives $\alpha X^2\ge4VX+23/6$, so $X\ge4V/a_5$. Per tile, $X_{\min}=\max(1455/V_{hi},4V_{lo}/a_5)$.

**Monotonicity in $X$** (verified term by term in GAP_AUDIT012.md §4):

- $k_3,k_1\propto X^{-1/2}$ and $k_4\propto X^{-1}$;
- $\bar R\propto W$;
- the erfc arguments are non-decreasing in $X$;
- the window sets are taken at $\eta_0(X_{\min})$;
- the crude terms decay.

So one evaluation at $X_{\min}$ per tile covers every $X\ge X_{\min}$, that is, every $(n,m)$ whose saddle lies in the tile. **No monotonicity in $t$ is used.**

**Certificate. [certified]** `cls012_lap.py 2.0 12.0 0.02` checks 500 tiles on $V\in[2,12]$ with 0 failures. The tiles are contiguous. The script steps with `v2=min(v+w,Vb)` and snaps the last tile to $V_b=12$ when float drift leaves it short, so the last tile is $[11.98,12]$. (The outside referee found that, before the snap, the last tile ended at 11.999999999999831 while the tail started at 12; see Appendix B.) The minimum margins are:

| class | min margin | $\|a\|$ |
|---|---|---|
| 0 | 3.0700 | 3.618 |
| 1 | 1.6879 | 2.236 |
| 2 | 0.8339 | 1.382 |

The worst $\delta$ is 0.137, at $V=9.78$, $X_{\min}=148.6$. These numbers come from `cls012_lap_rerun.log`, rerun after all fixes (the audit fixes and the snap). The tile margins are identical to the original run.

**Tail cell $V\ge11.99$. [proved + certified]** Write $\Xi=a_5/\zeta_\eta+a_5\zeta_\eta+\Delta$. On the strip $|\operatorname{Im}\eta|\le1/2$ we have $\operatorname{Re}\zeta_\eta\ge1/2$, so $|\gamma(V\zeta_\eta)/\zeta_\eta|\le\frac{12}{25}\mathrm{Li}_2(e^{-V/2})$. Cauchy estimates on discs of radius 1/2 give $\Xi_2\in2a_5\pm8\varepsilon$, $|c_3|\le a_5+8\varepsilon$ and $|\Xi''''|\le24a_5+384\varepsilon$. The mid-region bound is $a_5\eta^2/(1+\eta^2)-\frac{12}{25}\mathrm{Li}_2$. Also $M\ge|a|-4(e^{\bar\omega}-1)$ and $X\ge4V/a_5\ge4V_T/a_5$. Every $V$-dependent quantity moves in the safe direction as $V$ grows, so a single evaluation at $V_T$ covers all $V\ge V_T$. The cell is run from $V_T=11.99$, so it overlaps the last tile on $[11.99,12]$. Coverage at $V=12$ therefore does not depend on the snap. The margins are **3.0950 / 1.7130 / 0.8589**, with $X_{\min}=182.2$ (`cls012_lap.py tail 11.99`, `cls012_lap_rerun.log`: `tail cell V >= 11.99`).

**Validation.**

- `cls012_asym_check.py`: $c/(\mathrm{Amp}\cdot M)=0.990$–$0.998$ at $n=60,100$.
- The audit's independent check (GAP_AUDIT012.md §1): 45 exact points, $n\in\{291,437,500,700,1500,2500\}$, $V=2.08$–$40$, including $t\approx4$ and $\mu\approx0.105$. Every sign was correct. $|c/\mathrm{Amp}-M|$ ranged over 0.0005–0.0185, against a certified bound of 0.27–0.54, a margin of 15–30×.

### 6.4 Classes 3 and 4: the re-centred (thin) saddle

For classes 3 and 4, $\sum_hA_h\zeta^{h\cdot\mathrm{cls}}=a(\mathrm{cls})=0$, so the constant term of the twisted sum cancels. The survivor carries $e^{-5nu}=e^{-\sigma}$, which is an $O(N)$ term *in the exponent*. The withdrawn certificate treated this term as an amplitude (Appendix B.3). Here it is moved into the exponent.

**The thin function. [proved]** Put
$$g(s)=e^{s}\sum_{h=1}^4A_h\zeta^{h\cdot\mathrm{cls}}\exp(\Omega_h(s)),\qquad s=\sigma=V\zeta_\eta .$$
Since $\zeta^{hm}=\zeta^{h\cdot\mathrm{cls}}$, we have $\sum_hA_h\zeta^{hm}e^{\Omega_h(\sigma)}=e^{-\sigma}g(\sigma)$ exactly. The function $g$ has three properties.

1. **Leading coefficient.** Expand $e^{\Omega_h}=1+\Omega_h+\dots$. The constant term cancels. The coefficient of $e^{-s}$ is
   $$f_\infty=-\sum_hA_h\zeta^{h\cdot\mathrm{cls}}\sum_cB_1(c/5)\zeta^{-hc}=-1$$
   exactly, for both classes. This is checked in Arb to $10^{-37}$ (`lap34_lib.finf`).
2. **Uniform bound.** For every $s$ with $\operatorname{Re}s=V$,
   $$|g(s)-f_\infty|\le\varepsilon_g(V):=4e^V\big[0.8(\mathrm{Li}_1(\xi)-\xi)+\bar\omega^2e^{\bar\omega}/2\big].$$
   This follows from $|\Omega_h-\Omega_h^{(1)}|\le0.8\sum_{r\ge2}\xi^r/r$ and $|e^\Omega-1-\Omega|\le|\Omega|^2e^{|\Omega|}/2$, where $\Omega_h^{(1)}$ is the first-order term. $\varepsilon_g$ is $\xi^{-1}$ times a series in $\xi$ with non-negative coefficients starting at $\xi^2$, so it decreases in $V$. For example $\varepsilon_g(11.99)=1.787\cdot10^{-5}$ (`lap34_cert.log`, `epsg`). The audit tested this at 6000 points and found $\max|g+1|/\varepsilon_g=0.233$.
3. **Reality.** $g(V)$ is real for real $V$, because $h\leftrightarrow5-h$ are conjugate. The Arb enclosures of $g$ on the $V$-tiles lie in $[-1.13,-0.96]$.

**Saddle.** Choose $W$ by
$$\alpha(V)X^2=m-\tfrac16-5n .$$
The exponent $a_5/u-\gamma/u+(m-1/6)u-5nu$ then equals $X\Xi_V(\eta)$, with the same $\Xi_V$ as in §6.3. So
$$c(m)=\mathrm{Amp}\,[g(V)+\Delta],\qquad \mathrm{Amp}=W^{3/2}e^{X\Xi_V(0)}/\sqrt{2\pi\Xi_2}>0,$$
and the certificate proves $|\Delta|<-g(V)$, which gives $c(m)<0$. All the $O(N)$ parts of $\log F_m$ are inside $X\Xi_V$: the kernel, the twisted-product term $-\gamma/u$ (the "$e^{Dz}$" factor), and the thin shift $-5nu$. The remaining $\log g(V\zeta_\eta)$ is $O(1)$. It is kept out of the Gaussian and paid for in the window error through $\sup|g'|$.

**Error terms** (relative to Amp, all non-increasing in $X$):

- **Window.** With $G(\eta)=e^{V\zeta_\eta}\sum_hA_h\zeta^{hm}e^{\Omega_h+R_h}$, $\psi=X(c_3\eta^3+R_4)$ and $c_3$ purely imaginary,
  $$|e^\psi G-g(V)|\le e^{k_4s^4}(V|\eta|\sup_{\rm win}|g'|+\varepsilon_R)+|g(V)|\big(\min(k_3s^3,2)e^{k_4s^4}+e^{k_4s^4}-1\big).$$
  In this bound:
  - $\varepsilon_R=4e^{\bar\omega}e^V(e^{\bar R}-1)\ge|G-g(V\zeta_\eta)|$;
  - $\sup|g'|$ is the smaller of a ball evaluation (sub-balls of width at most 0.05 in $\operatorname{Im}s$) and the Cauchy bound $\varepsilon_g(V-1)$;
  - the integral is an upper Riemann sum over 800 $s$-cells;
  - the Gaussian truncation contributes $|g|\operatorname{erfc}(cc/\sqrt2)$.
- **Mid region, $\eta_0<|\eta|\le2$.** The bound $\operatorname{Re}(\Xi(0)-\Xi(\eta))\ge d_k\eta^2$ on 170 $\eta$-tiles. On each tile, $\sup|G|\le\sup|g(V\zeta_\eta)|+\varepsilon_R$. Tiles whose erfc argument exceeds 40 are bounded jointly, via $4e^{V+\bar\omega+\bar R}\cdot170\sqrt{\Xi_2/2d_{\min}}\operatorname{erfc}(40)$.
- **Crude terms.** As in §6.3, but with an extra factor $e^V$, because Amp carries $e^{-V}$. The rate is $\kappa=a_5-\gamma-0.8\mathrm{Li}_2(\xi)$. Tiles use $e^{V_{hi}}$; the tail cell uses $V\le a_5X/3$. The script asserts decay $c\ge\kappa-\max(a_5/5,a_5/4)$ (tiles), or $c\ge\kappa-a_5/3-a_5/4$ (tail), and $X_{\min}>4.5/c$. $E_2\le10^{-8}$ everywhere.

**Range of $X$.** $X=5n/V\ge4105/V$ for $n\ge821$. $t\ge4$ gives $m-1/6-5n\ge15n+23/6>3VX$, so $X>3V/a_5$. Per tile, $X_{\min}=\max(4105/V_{hi},3V_{lo}/a_5)$. For the tail cell, $X\ge X^*=\sqrt{15\cdot821/a_5}=216.3$.

**Certificate. [certified]** `lap34_cert.py`:

| region | cells | failures | min margin $-g-\|\Delta\|$ (class 3 / class 4) |
|---|---|---|---|
| $V\in[2,12]$, width 0.02 (last tile snapped to $[11.98,12]$) | 500 | 0 | **0.8363 / 0.8361** (tile $[11.98,12]$: win 0.119, mid 0.005, $\varepsilon_R$ 0.042) |
| $V\ge11.99$ (analytic tail cell, $V_T=11.99$) | 1 | 0 | **0.8146 / 0.8146** (0.81456) |

The tail cell uses the strip bounds of §6.3 together with $|g(V\zeta_\eta)-g(V)|\le2\varepsilon_g(V_T)$ and $|g(V)|\ge1-\varepsilon_g(V_T)$, where $\varepsilon_g(V_T)=\varepsilon_g(11.99)=1.787\cdot10^{-5}$. Since $\varepsilon_g$ decreases, these bounds hold for every $V\ge11.99$. The tail cell overlaps the last tile on $[11.99,12]$.

The certified relative error is therefore at most 0.19.

All these numbers are in `lap34_cert.log`, which was produced by the current script. That script has the audit fixes: the snapped last tile, the tail from $V_T=11.99$, the `_fin` finiteness guards and the NaN-safe pass test. The second audit's own rerun reproduced them (AUDIT2_LAP34.md §6).

**Validation.**

- `lap34_validate.py`: 70 exact points at $n=150,300,821,2000$, with $\mu$ up to 0.1050 and $t$ down to 4.0. The true $|c/\mathrm{Amp}-g|$ lies in $1.8\cdot10^{-4}$–$1.6\cdot10^{-2}$, which is 24–435× below the certified bound. Every sign is negative.
- The first audit's independent check (AUDIT_LAP34.md): 36 points at $n=821,1000,1200$, with $\mu$ up to 0.10500. The error is 28–579× below the bound of the actual certified tile, and every sign is negative.
- The second audit (AUDIT2_LAP34.md §4) used a third, independent exact method: $S_n=G_5E_n$ with $E_n$ from two Euler identities, computed modulo 31-bit primes and combined by CRT. It shares no code with `lap34_exact.py` or with the first audit. It checked 31 points at $n=821$, 2000 and 5000, including $\mu=0.10500$ and $t=4.000$ at $n=821$, pairs straddling the tile/tail junction at $V=11.99$ and 12.01, and $V$ up to 46.8. Every sign is negative. The true error is 28–579× below the bound of the actual certified tile or tail cell.

### 6.5 Coverage lemma for the saddle [certified]

Write $A(V):=\alpha(V)/V^2$. Then:

- $A$ is continuous, $A\to1/5$ as $V\to0$, and $A\to0$ as $V\to\infty$;
- $\alpha$ is increasing, since $\alpha'=Vq(V)/5$ with $q=1/(e^V-1)-5/(e^{5V}-1)>0$.

`cls012_cov.py` certifies $A\ge0.04424$ on $[0.001,2]$ (5858 tiles, using $A\ge\alpha(v_1)/v_2^2$) and $A\ge0.1998$ on $(0,0.001]$.

- **Classes 0–2.** The saddle equation is $A(V)=(m-1/6)/(25n^2)<0.4\mu$.
- **Classes 3–4.** The saddle equation is $A(V)=(m-1/6-5n)/(25n^2)$, which is smaller still.

When $\mu<0.105$, the right-hand side is $<0.042$. By the intermediate value theorem a root exists, and every root has $V>2$. The representation is exact for every $W>0$, so any root may be used.

---
## 7. The centre, $\mu\ge0.105$

*(Sources: GAP_CENTRE2.md §2 for the near-peak part; GAP_OFFPEAK.md for the off-peak lemma; AUDIT_CENTRE.md, which audited the near-peak part and replaced the combination step; AUDIT_OFFPEAK.md, which audited the off-peak lemma.)*

Throughout this section $n\ge821$ and $0.105\le\mu\le\tfrac12$. We use:

- $r=e^{-L/(5n)}$ and $q=\zeta^he^{-\sigma/(5n)}$, with $\sigma=L-iy$, so that $y=5n(\theta-2\pi h/5)$;
- $\Phi_c(s)=\log(1-\zeta^{hc}e^{-s})$ (principal branch) and $Q=\sum_{c=1}^4\Phi_c=\log\frac{1-e^{-5s}}{1-e^{-s}}$;
- $p(\sigma)=\int_0^1Q(t\sigma)\,dt$;
- $C_h(\sigma)=\sum_c(c/5-\tfrac12)(\Phi_c(\sigma)-\Phi_c(0))$.

The starting point is Cauchy's formula for the polynomial $S_n$:
$$c(m)=\frac1{2\pi}\int_{-\pi}^{\pi}S_n(re^{i\theta})\,r^{-m}e^{-im\theta}\,d\theta,\qquad r^{-m}=e^{2n\mu L}.$$

### 7.1 Near-peak piece (A): $|y|\le1.1$ around $\theta=2\pi h/5$, $h=1,\dots,4$ [certified]

**Lemma 7.1 (finite Euler–Maclaurin). [proved]** For $0\le a<1$,
$$\sum_{t=0}^{n-1}f(t+a)=\int_0^nf+B_1(a)(f(n)-f(0))+\frac{B_2(a)}2(f'(n)-f'(0))-\int_0^n\frac{\tilde B_2(x-a)}2f''(x)\,dx .$$

*Proof.* Apply the Euler–Maclaurin formula of order 2 with offset $a$ on each $[t,t+1]$ and sum over $t$; the argument is the same as for Lemma 6.1, but on a finite range. ∎

The constants are $B_2(1/5)=B_2(4/5)=1/150$ and $B_2(2/5)=B_2(3/5)=-11/150$. Since $\tilde B_2\in[-1/12,1/6]$, the last term is at most $\frac1{12}\int|f''|$. The identity was checked numerically to $10^{-15}$.

**Lemma 7.2. [proved]** Apply Lemma 7.1 to $f(x)=\Phi_c(x\sigma/n)$ with $a=c/5$, and sum over $c$. The terms $k=5t+c$, $t<n$, are exactly the factors of $S_n$. This gives
$$\big|\log S_n(q)-np(\sigma)-C_h(\sigma)\big|\le\frac{K_h(\sigma)}n,$$
$$K_h=\sum_c\Big[\frac{|B_2(c/5)|}2|\sigma|\big(|\Phi_c'(\sigma)|+|\Phi_c'(0)|\big)+\frac{|\sigma|^2}{12}\sup_{t\in[0,1]}|\Phi_c''(t\sigma)|\Big],$$
valid for $|y|<2\pi/5$ and $L\ge0$, uniformly down to $L=0$. The reason is that $\Phi_c$ has no singularity on the segment $[0,\sigma]$, because $\zeta^{hc}\ne1$.

Two exact facts are used in place of the old arctan law and its equal-modulus assumption:

- $|e^{C_h(L)}|=1$ for real $L$. This is enclosed in the certificate, not assumed.
- The four peaks have equal modulus exactly.

The class margin is then
$$M_s(L)=\operatorname{sgn}_s\operatorname{Re}\sum_h\zeta^{-hs}e^{C_h(L)},\qquad \operatorname{sgn}_0=+1,\ \operatorname{sgn}_{1..4}=-1 .$$

**Saddle.** $\mu=-p'(L)/2$. Since $Q'(0)=\sum_c\zeta^c/(1-\zeta^c)=-2$, we get $p'(0)=-1$ and $\mu(0)=1/2$. `centre3_saddle_end.py` certifies $\mu(2.1)<0.105$, by the real mean-value theorem on 4096 $t$-pieces. Its log `centre3_saddle_end.log` prints `mu(2.1) = -p'/2 in [0.104 +/- 3.02e-4]  < 0.105: True`. (The second audit computes $\mu(2.1)=0.103772$ from the closed form, so $L(0.105)=2.0814$.) The certificate asserts $a=p''(L)>0$ on every $L$-box, so $\mu(L)$ is strictly decreasing. Therefore every $\mu\in[0.105,1/2]$ has a unique saddle $L\in[0,2.1]$.

**Laplace step.** Put
$$\mathcal N(n,L)=\frac{e^{np(L)+2n\mu L}\sqrt{2\pi/(na)}}{10\pi n}.$$
Then
$$c(m)=\mathcal N(n,L)\Big[\sum_{h=1}^4\zeta^{-hm}e^{C_h(L)}(1+\rho_h)+\mathrm{strip}+\mathrm{OFF}\Big].$$
The ingredients are:

- $d\theta=dy/(5n)$;
- the linear phase cancels at the saddle, since $-inp'(L)y=2in\mu y$;
- each peak contributes $\mathcal N\zeta^{-hm}e^{C_h}$.

The audit re-derived this representation (AUDIT_CENTRE.md item 5). The three error pieces are bounded as follows.

- **Inner window $|y|\le0.35$.** Write $\psi(y)=p(L-iy)-p(L)+ip'(L)y$ and $Z=n\psi+nay^2/2+\Delta C$. The odd part $T=n\psi_3y^3+C_1y$ integrates to exactly 0 against the Gaussian. The rest is bounded using:
  - $|e^Z-1-Z|\le\frac{|Z|^2}2\max(1,e^{\operatorname{Re}Z})$;
  - $\operatorname{Re}\psi\le-by^2/2$, where $b$ is a certified lower bound for $\operatorname{Re}p''$ on the window;
  - Taylor constants from $p^{(k)}(\sigma)=\int_0^1t^kQ^{(k)}(t\sigma)\,dt$, enclosed on $t$-pieces;
  - the factor $e^{K/n}-1$ and the Gaussian tail.
- **Strip $0.35\le|y|\le1.1$.** Lemma 7.2 is combined with $\operatorname{Re}\psi(y)=-\int_0^{|y|}(|y|-u)\operatorname{Re}p''(L-iu)\,du$, with $\operatorname{Re}p''$ bounded below on boxes. This gives terms of the form $\sqrt n\,e^{-nD}$, which are non-increasing for $n\ge1/(2D)$; that condition is asserted. The strip and the Gaussian tail are normalised with the upper curvature $a_{up}$, which is the safe direction. (The withdrawn certificate had used a lower bound here.)
- **Uniformity in $n$.** $n$ enters only through $K/n$, the $1/n$ moment terms, the Gaussian tail $e^{-nAY^2/2}/\sqrt n$, and the strip. Each of these is non-increasing in $n$. So a bound at $n=N_0=821$ holds for all $n\ge821$.

**Grid.** 210 $L$-boxes of width 0.01 on $[0,2.1]$. Each has 28 inner $y$-boxes, 60 strip $y$-boxes and 64 $t$-pieces, and every quantity is enclosed over its full box (`centre3_peak.py 821 210 64`). Since $L(0.105)=2.0814$, only the boxes up to $[2.08,2.09]$ are needed; the last box $[2.09,2.10]$ is certified as well.

**Result at $N_0=821$** (`centre3_peak_821.log`; reproduced in `centre3_comb_821.log`). For every box and every class, $M_s>\sum_h|e^{C_h}|(\rho_h+\mathrm{strip}_h)=:\mathrm{err}_s$. The worst ratios err/margin are:

| class | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| worst err/margin | 0.0566 | 0.2225 | 0.2241 | 0.7977 | 0.8469 |

The minimum slack $\min_s(M_{s,lo}-\mathrm{err}_{s,up})$ is 0.018703, attained for class 4 on the box $L\in[2.09,2.10]$.

### 7.2 The off-peak lemma (B) [certified]

**Lemma 7.3 (B).** Let $n\ge821$ and $L\in[0,2.1]$, and let $\theta$ satisfy $|5n(\theta-2\pi h/5)|>1.1$ for $h=1,\dots,4$, taken mod $2\pi$. Then
$$\log|S_n(e^{-L/(5n)}e^{i\theta})|-np(L)\le-E(n),\qquad E(n)=1.5\log n+6.3\quad(E(821)=16.36).$$

Every certificate uses the constant $6.3$ as the exact ball `arb('6.3')`. (An earlier version of (E), (E0) and (N-a) used the float $\mathrm{fl}(6.3)$, which is $1.8\cdot10^{-16}$ low. Those pieces were re-run with the exact constant, and the logs quoted below come from the re-run.)

**Closure.** For fixed $n$, $\log|S_n(e^{-L/(5n)}e^{i\theta})|$ is continuous in $(\theta,L)$, with values in $[-\infty,\infty)$. So each bound "$\le np(L)-E(n)$" that is certified on a box also holds on the closure of that box. This closes the float-endpoint slivers the audits found:

- $|y|\in(1.1,\mathrm{fl}(1.1))$, at the start of (E) (AUDIT2_OFFPEAK.md N1);
- $(\mathrm{fl}(\pi),\pi]$, at the end of (F) (AUDIT_OFFPEAK.md D1, which the Arb run also closes by ending at `nextafter(pi, 4)`);
- $L\in[0.1,\mathrm{fl}(0.1))$, at the lower end of the $L$-boxes of (E), (E0) and (N), where $\mathrm{fl}(0.1)=0.1+5.5\cdot10^{-18}$ (AUDIT2_OFFPEAK.md N2).

The first two are closed by this argument. The third is closed directly: the lowest $L$-box of (E), (E0) and (N) now starts at $0.0999$, below the exact $0.1$, so $L=0.1$ is inside a certified box and the Lemma M reduction from $L<0.1$ to $L'=0.1$ costs exactly $n(0.1-L)\le 0.1\,n$, which is what each piece charges (`arb('0.1')`). All three pieces were re-run after this change and remain certified.

Since $S_n$ has real coefficients, $|S_n(\bar q)|=|S_n(q)|$, so it suffices to take $\theta\in[0,\pi]$. The pieces (E) and (N) use both signs of $y$ near $h=1,2$, which also covers $h=3,4$.

**Closed forms. [standard]** For $h\ne0$,
$$p(\sigma)=\frac{\mathrm{Li}_2(e^{-5\sigma})/5-\mathrm{Li}_2(e^{-\sigma})+4\zeta(2)/5}{\sigma},$$
using $\sum_{c=0}^4\mathrm{Li}_2(\zeta^cx)=\mathrm{Li}_2(x^5)/5$. For $h=0$, $p_0(\sigma)=4(\mathrm{Li}_2(e^{-\sigma})-\zeta(2))/\sigma$. $Q$ is decreasing and convex with $|Q'|\le2$: indeed $Q''=[g(s/2)-g(5s/2)]/s^2$ with $g(x)=x^2/\sinh^2x$ decreasing. Since $p'=(Q-p)/L$, $p$ is decreasing on $L>0$.

**Lemma M (monotonicity in $L$). [proved]** For $0\le L\le L'$,
$$\log|S_n(re^{i\theta})|\le\log|S_n(r'e^{i\theta})|+n(L'-L),\qquad r=e^{-L/(5n)},\ r'=e^{-L'/(5n)}.$$

*Proof.* We have $|1-\rho e^{ia}|^2=(1-\rho)^2+4\rho\sin^2(a/2)$. For $\rho'\le\rho\le1$, this gives $|1-\rho e^{ia}|^2\le(\rho/\rho')|1-\rho'e^{ia}|^2$. Take $\rho=r^k$, $\rho'=r'^k$, and use $\sum_{k\le5n,\,5\nmid k}k=10n^2$: the total factor on $\log|S_n|$ is $\frac12\cdot\frac{L'-L}{5n}\cdot10n^2=n(L'-L)$. ∎

Pieces (E), (E0) and (N) are certified on $L\in[0.1,2.1]$. For $L<0.1$ they use Lemma M at $L'=0.1$ together with $np(L)\ge np(0.1)$ ($p$ decreasing, and $y$ is unchanged because $\theta$ is fixed). The cost is $0.1n$, which the boxes touching $L=0.1$ are required to absorb.

**Piece (E): $h\in\{1,2\}$, $1.1\le|y|\le6$. [certified]** $\Phi_c$ is analytic on $\operatorname{Re}s>0$ and at $s=0$. So Lemma 7.1 applies for every $y$, even past $|\sigma|=2\pi/5$; the disc was needed only for the Taylor step in §7.1. This gives
$$\log|S_n|\le n\operatorname{Re}p(\sigma)+\operatorname{Re}C_h(\sigma)+K_h(\sigma)/n,$$
with the remainder in integral form. The certificate uses 20000 adaptive boxes. On each box it checks two things:

- $V(N_0)=-N_0g+\operatorname{Re}C+K/N_0+E(N_0)<0$, where $g=p(L_{hi})-\operatorname{Re}p(\sigma)-0.1\chi$;
- $g>1.5/N_0$, so that $V(n)$ is non-increasing in $n$.

The worst values are $V=-2.49,-2.27$ ($h=1$) and $-2.55,-2.18$ ($h=2$); an earlier run quoted $-12.03$, before the lowest $L$-box was widened to start at 0.0999. The largest $K$ is $5.6\cdot10^3$, near a singular crossing at $L=0.1$ (`offpeak_em.py 821 6.0 32`, `offpeak_em_821.log`).

**Piece (E0): $h=0$, $|y|\le20$. [certified]** Here $\Phi=\log(1-e^{-s})$ is singular at 0, so the $t=0$ terms are split off and Lemma 7.1 is applied on $t=1,\dots,n-1$; note $\sum_cB_1(c/5)=0$. This gives
$$\log|S_n|\le n\operatorname{Re}p_0(\sigma)+\Psi(\omega)+0.08\big(|\omega\Phi'(\omega)|+|\omega||\Phi'(\sigma)|\big)+\frac{|\omega||\sigma|}3\int_{1/n}^1|\Phi''(t\sigma)|\,dt,$$
where $\omega=\sigma/n$, $\Psi(\omega)=4+\log(24/625)+\sum_c\operatorname{Re}\ell(c\omega/5)-4\operatorname{Re}[\omega^{-1}\int_0^\omega\ell]$ and $\ell(s)=\log((1-e^{-s})/s)$. The $\log|\omega|$ terms cancel exactly. For $\operatorname{Re}s\ge0$ and $|s|\le0.1$ the following are proved:

- $|(1-e^{-s})/s-1|\le|s|e^{|s|}/2$, hence $|\ell(s)|\le|s|$;
- $|s\Phi'(s)|\le1/(1-|s|e^{|s|}/2)$;
- $|\Phi''(s)|\le1/(|s|^2m^2)$, with $m=1-re^r/2$.

Consequently $\Psi\le0.7403+6|\sigma|/N_0$. The certificate uses 4945 boxes, and the worst value is $V=-0.23$; this is a bisection artefact, as the true slack is in the hundreds (`offpeak_em0.py 821 20 32`, `offpeak_em0_821.log`).

**Piece (N): $6\le|y|\le0.25n$ ($20\le y$ for $h=0$). [certified]** Since $|q|<1$,
$$\log|S_n|=-\operatorname{Re}\sum_jG_j/j\le\sum_j|G_j|/j,\qquad G_j=\frac{(1-e^{-j\sigma})A_j}{1-e^{-j\omega}},\qquad A_j=\sum_c\zeta^{hcj}e^{-cj\omega/5}.$$
The bounds on $A_j$ are $|A_j|\le\min(4,\bar A_j+2j|\omega|)$ and $|A_j|\le\sum_ce^{-cjL/(5n)}$, where $\bar A_j=|\sum_c\zeta^{hcj}|$ is 1 or 4.

Here $\omega=\sigma/n$, and $A_j$ is the sum over the four $c$ for fixed $j$ (the factor $1-e^{-j\sigma}$ comes from summing the geometric series over $t<n$). Write $\alpha(a)=\sum_ce^{-ca/5}$, so the second bound reads $|A_j|\le\alpha(jL/n)$.

- **Part I, $j\le\pi n/y$.** Use $1/|1-e^{-z}|\le1/|z|+c_1$ with $c_1=1+1/\pi$. This is the maximum modulus of $1/(1-e^{-z})-1/z$ on the rectangle $[0,20]\times[-\pi,\pi]$: it is at most 0.593 on the left edge, $1+1/\pi$ on the top and bottom edges, and 1.06 on the right edge. Every $z=j\omega$ in Part I lies in the rectangle, since $|\operatorname{Im}z|=jy/n\le\pi$ and $\operatorname{Re}z=jL/n\le\pi L/y\le1.1$. This gives Part I $\le nM(\sigma)+R_1$, where
  $$M(\sigma)=\sum_j\frac{|1-e^{-j\sigma}|\bar A_j}{j^2|\sigma|},\qquad R_1\le B\big(1+\log(\pi n/y)\big)+O(1),$$
  with $B=2+1.6c_1$ for $h\ne0$ and $B=2+4c_1$ for $h=0$.
- **Part II, $j>\pi n/y$.** Split the $j$ into periods $k\ge1$, with $jy/n\in((2k-1)\pi,(2k+1)\pi]$. Write $a=jL/n$ and $b=jy/n$, so that $j\omega=a-ib$. On period $k$:
  - $1/j\le y/((2k-1)\pi n)$;
  - $|1-e^{-j\omega}|\ge\max\big(1-e^{-a},\ 2e^{-a/2}|\sin(b/2)|\big)$;
  - $|\sin(b/2)|\ge|b-2\pi k|/\pi$ on the period.

  The two $j$ nearest the pole $b=2\pi k$ are bounded by $G_k=1/(1-e^{-a_k^-})$, where $a_k^-$ is the smallest $a$ on the period. The other $j$ on the period are bounded through $|\sin(b/2)|$ and summed with
  $$\sum_{i\le M}\min(G,D/i)\le\min\big(MG,\ D(1+\log^+(2MG/D))\big).$$
  This has three consequences:
  - summed over $k$, the pole terms total at most $\beta u^2n/L$, with $u=y/n$ and $\beta=1+e^{-\pi L/u}$;
  - the remaining terms of the periods total $O(\log^2n)$;
  - the periods with $a_k^->2$ give a constant tail, bounded with $\int_2^\infty\phi(a)/a\,da$ in closed form.

  Part II is non-increasing in $n$ at fixed $(L,y)$: $M/n$, $D/n$ and $\log^+(2MG/D)$ do not increase, and $\beta_k$ decreases (AUDIT2_OFFPEAK.md §3). (N) stops at $u=0.25$ because the pole term $\beta u^2/L$ must stay below $p(L)-0.1$ at $L=0.1$. The triangle-bound majorant blows up at low-denominator rationals such as $\theta=\pi$ and $2\pi/3$, and that is where (F) takes over. *(Source: GAP_OFFPEAK.md §5.)*

The certificate has two parts:

- **(N-a)** $y\in[Y_E,200]$ on $(L,y)$-boxes, at $N_0$, with slope check $g\ge(B+1.5)/N_0$;
- **(N-b)** $y\in[200,0.25n]$ on $(L,u)$-boxes, in the form $-nG+a+b\log n+c\log^2n$ with $G\ge(b+2c\log N_0)/N_0$.

- **(N-a) in more detail.** It uses the exact $M(\sigma)$ (400 terms plus a tail), evaluated at $n=N_0$. Monotonicity in $n$ follows from $R_1=A+B\log n$, Part II being non-increasing, and the checked slope condition.

The worst values are $-147.9/-23.3$ for $h=1$ and $-71.3/-1.8$ for $h=0$ (`offpeak_tb.py 821 0.25`, `offpeak_tb_821.log`). The second audit compared the certified Part I and Part II with the exact $j$-series at 70 points and found 0 violations.

The script runs only $h\in\{0,1\}$ with $y>0$. This suffices: the majorant depends on $h$ only through $\bar A_j$, which equals 1 if $5\nmid j$ and 4 if $5\mid j$ for every $h=1,\dots,4$, and it depends on $y$ only through $|1-e^{-j\sigma}|$, $|\sigma|$ and $|y|$, all of which are invariant under $y\mapsto-y$ (AUDIT_OFFPEAK.md D2).

**Piece (F): $\operatorname{dist}(\theta,\{0,2\pi/5,4\pi/5\})\ge0.0499$, $\theta\le\pi$, $L\in[0,2.1]$. [certified in Arb; also float-interval]**

The chunk length is $T=40$ on the bands $0.05\le\operatorname{dist}<0.1$ and $T=20$ elsewhere. Put $M=\lfloor n/T\rfloor$, $\lambda=TL/n$, $\mu_r=L/(5n)$ and $P(V)=\prod_{k'\le5T,\,5\nmid k'}(1-Ve^{ik'\tilde\theta})$ with $\tilde\theta=\theta+i\mu_r$. Here the radial drift $r^{k'}$ appears as an imaginary shift of $\theta$.

*Chunk lemma [proved]:*
$$\log|S_n|\le\sum_{m<M}\mathcal B(m\lambda)+4(T-1)\log(1+e^{-(L-\lambda)}),\qquad\mathcal B(s)=\max_{|V|=e^{-s}}\log|P(V)|,$$
$$np(L)\ge\sum_{m<M}TQ(m\lambda)-T\big[\lambda+(\log5-Q(L))/2\big].$$
*Proof of the first line.* Write $n=MT+r$ with $0\le r<T$. For $m<M$, the factors with $5Tm<k\le5T(m+1)$ form one chunk: with $k=5Tm+k'$ and $V_m=q^{5Tm}$, the chunk equals $P(V_m)$ with $|V_m|=e^{-m\lambda}$, so its log-modulus is at most $\mathcal B(m\lambda)$. This holds for **every** $n$, whether or not $T\mid n$. The factors with $5TM<k\le5n$ and $5\nmid k$ are left over; there are $4r\le4(T-1)$ of them. Each satisfies $|1-q^k|\le1+e^{-kL/(5n)}$, and $kL/(5n)>TML/n>L-\lambda$ because $TM>n-T$. So the leftover factors cost at most $4(T-1)\log(1+e^{-(L-\lambda)})$.

*Proof of the second line.* It uses the tangent line of $Q$ at each $m\lambda$, convexity telescoping $-Q'(m\lambda)\lambda\le Q((m-1)\lambda)-Q(m\lambda)$, $|Q'|\le2$ for the $m=0$ term, $Q>0$, and $M\lambda\le L$ (AUDIT2_OFFPEAK.md §1).

*Reductions [proved]:*

- $\mathcal B$ is convex and non-increasing in $s$ (Hadamard three circles). So on an $s$-grid of step 0.1, $\mathcal B$ is below its chord and $TQ$ is above its tangent, which gives cell lower bounds $d_i$.
- $\sup_V\log|P|$ is subharmonic in $\tilde\theta$, since $P$ is entire in $\tilde\theta$. So on each rectangle $[\theta_a,\theta_b]\times[0,\mu_{\max}]$, with $\mu_{\max}=2.1/(5N_0)$, the maximum is attained on the boundary, and the drift needs no interior treatment.
- The assembly uses left Riemann sums of the running minimum of $d_i$; for small $L$ it uses $M\,d_{\min}$ instead. It is checked on the $L$-intervals $[i/100,(i+1)/100]$, $i=0,\dots,209$, which cover $[0,2.1]$ including $L=0$. The check is done at $N_0$, with the slope condition, so it holds for all $n\ge N_0$.

*Result:* 163 rectangles and 22 $s$-values, minimum assembled margin **47.9** (near $\theta=\pi$, $T=20$). The tightest $T=40$ band has margin 51.8 (`offpeak_chunk.py 821 all`, 470 s).

*Arithmetic:* float64 intervals, using only correctly rounded $+,-,\times,\div$ and exact frexp/ldexp. All transcendental inputs come from Arb, with explicit rounding margins. The margins are:

- $10^{-12}$ on phase half-widths and on cos;
- $32u$ on each $|1-\rho e^{i\phi}|^2$;
- a factor $(1+2u)^{2(K+3G+8)}$ on the product.

No libm function enters a bound.

*Arb certificate.* `offpeak_chunk_arb.py` redoes (F) with Arb balls over dyadic cells. It also fixes two defects of the float run: D1 (the last rectangle ends at `nextafter(pi, 4)`, not at float $\pi$, which is $1.2\cdot10^{-16}$ below $\pi$) and D3 (exact Arb comparisons in the assembly). All 163 rectangles of pieces 0–7 are OK. The summary `arb_F_summary_full.log` prints `(F) arb: CERTIFIED for all n >= N0 (all rects OK)`, with minimum margin **47.91** (piece 7, rectangle 26, ending at $\pi^+$), identical to the float run. It merges 13 logs:

- `arb_F_p0.log`, `arb_F_p0_rest.log`;
- `arb_F_p1.log` to `arb_F_p4.log`;
- `arb_F_p5a.log`, `arb_F_p5_1-13.log`, `arb_F_p5_14-19.log`, `arb_F_p5_20-26.log`, `arb_F_p5b.log`;
- `arb_F_p6.log`, `arb_F_p7.log`.

(The file `arb_F_summary.log` is an incomplete earlier summary and is not the certificate.)

*The early exit.* The first runs stalled. `ArbChunk.arb_ok` had an early `return False` once the running product exceeded the threshold, but later factors can be $<1$ and bring the product back below it (ARB_F_RUN.md). The exit could only reject cells, never accept them, so every rectangle logged OK before the fix is still valid. After the exit was removed, the cell counts match the float run (e.g. piece 5, rectangle 1: 61,952 cells, margin 185.04). 132 of the 163 rectangles come from runs made before the fix. The second audit re-ran 11 of them with the current code, including the worst one (47.91), and got identical margins with no Arb failures (AUDIT2_OFFPEAK.md §4).

**Coverage of the $\theta$-circle** (AUDIT_OFFPEAK.md (a)). The pieces meet or overlap at every handover:

| handover | how they meet |
|---|---|
| (A)–(E) | abut at $\|y\|=1.1$ (the statement has a strict $>$) |
| (E)–(N) | abut at $\|y\|=6$ |
| (E0)–(N) | abut at $y=20$ |
| (N-a)–(N-b) | meet at $y=200$, which needs $0.25n\ge200$, i.e. $n\ge800$ |
| (N)–(F) | (N) reaches $\|\theta-2\pi h/5\|=0.25n/(5n)=0.05$ for every $n$, and (F) starts at 0.0499; overlap $10^{-4}$ rad |

The root $6\pi/5$ is 0.63 away from $\pi$, so $[0,\pi]$ needs no exclusion around it. The (F) rectangles tile each piece exactly, and the pieces abut exactly. The float-endpoint slivers are closed by the closure remark above.

*Spot validation.*

- **First audit** (AUDIT_OFFPEAK.md (f)). 312 random points with exact Arb full products, at $n=821$, 1000, 3000 and three random $n\in[821,2500]$, with all $h$ and both half-circles. Every point lies below its piece's certified bound. At $n=821$ the worst true value is $-43.2$, where the certified value is $-12.0$.
- **Second audit** (AUDIT2_OFFPEAK.md §5). 5,570 adversarial Arb points for $n$ from 821 to $10^5$. These include window edges, box endpoints, $L$ values down to $10^{-9}$, and located local maxima. There are 0 violations. The smallest true slack is 43.18, at $n=821$, $L=2.1$, $|y|=1.1^+$, and it grows with $n$.

### 7.3 Combination [certified, per $L$-box]

On the off-peak set, Lemma 7.3 gives $|S_nr^{-m}|\le e^{np(L)+2n\mu L-E(n)}$. The set has measure at most $2\pi$, so
$$|\mathrm{OFF}|\le\frac{1}{\mathcal N}\cdot\frac{2\pi}{2\pi}e^{np+2n\mu L-E(n)}=10\pi n\sqrt{\frac{na}{2\pi}}\,e^{-E(n)}=10\pi\sqrt{\frac{a}{2\pi}}\,e^{-6.3}.$$
The powers of $n$ cancel exactly, so this bound is independent of $n$. It grows with $a=p''(L)$.

*Correction (AUDIT_CENTRE.md D1).* GAP_CENTRE2.md §3 and PROOF.md §16 claimed "$|\mathrm{OFF}|\le0.0187\le\min(M_s-\mathrm{err})$". That claim is false. The largest $a$ occurs at $L=0$, where $a=2/3$ (box enclosure 0.6843), and there $|\mathrm{OFF}|\le0.018791$ (0.019039 with the enclosure). The minimum slack is 0.018703, but it occurs at $L\approx2.1$.

The valid comparison is **per $L$-box**. `centre3_comb.py` is `centre3_peak.py` plus one step per box:

- it computes $\mathrm{OFF}_{up}=10\pi\sqrt{a_{up}/2\pi}\,e^{-6.3}$ in Arb, with the exact constant `arb('-6.3')`;
- for every class $s$ it **asserts** $M_{s,lo}-\mathrm{err}_{s,up}-\mathrm{OFF}_{up}>0$ in Arb.

The comparison is therefore part of the certificate. It is not read off printed decimals, which is what the second audit's item D-a had asked for. The grid is run in two halves, boxes 0–104 and 105–209, and `centre3_comb_821.log` concatenates both. Each half ends with the line `COMBINATION ASSERTED (slack > OFF, arb) on every box; N0=821  boxes=210  all classes certified on every box: True`. Results:

- Over all 210 boxes and all five classes, the minimum of slack minus OFF is **0.010043**, on the box $[2.09,2.10]$: slack 0.018703, OFF 0.008660.
- The boxes actually needed stop at $[2.08,2.09]$, since $\mu(2.1)<0.105$. On that box slack minus OFF is **0.0120** (0.020699 − 0.008696).

**Conclusion 7.4.** For $n\ge821$, $\mu\in[0.105,\tfrac12]$ and every class, $c(m)$ has the conjectured sign, since Lemma 7.3 is certified in Arb for all four pieces (E), (E0), (N) and (F).

*Validation (not part of the proof).* The second audit (AUDIT2_CENTRE.md §5) evaluated $c(m)/\mathcal N$ by quadrature on the windows. It first validated the quadrature against exact coefficients at $n=300$ (35 points, agreement to $4\cdot10^{-11}$). It then checked 120 points at $n=821$, 857 and 900, with $\mu\in\{0.105,0.15,0.2,0.3,0.4,0.45,0.49,0.5\}$ and all five classes. Every sign is correct. The realised error is at most $1.65\cdot10^{-3}$, and at most $4.7\cdot10^{-4}$ in classes 3 and 4, against a certified err of 0.06–0.21.

---

## 8. Coverage

**Theorem 8.1.** Let $n\ge980$ and $0\le m\le5n^2$. Then $(n,m)$ is covered by at least one of the following. (The three regions are disjoint, but the certificates themselves also hold for smaller $n$, down to 291 or 821, where they overlap the exact computation. "At least one" is all that is needed.)

- the band certificate of §5 (if $t<4$);
- the Laplace certificates of §6 (if $t\ge4$ and $\mu<0.105$);
- the centre certificate of §7 (if $\mu\ge0.105$).

Each certificate that applies proves the conjectured sign of $c_n(m)$.

*Proof.* The three conditions partition $\{0\le m\le5n^2\}$: every $m$ has either $t<4$ or $t\ge4$, and in the latter case either $\mu<0.105$ or $\mu\ge0.105$. Note that $t\ge4$ is automatic when $\mu\ge0.105$, since $m\ge1.05n^2>4(5n+1)$ for $n\ge22$, but it is not needed. It remains to check that each certificate's hypotheses hold on its whole region.

1. **Band, $t<4$, i.e. $m<4N$.**
   - Classes 3, 4: Corollary 5.3 covers $m<2N$ for every $n$. For $2N\le m<4N$, the exact computation covers $n\le2000$ and $\Psi$ covers $n>2000$.
   - Classes 0, 1, 2: the exact computation covers $n\le2000$, and `cls012_band` covers $n>2000$.

   The only threshold is $n\ge291$ for the exact band run, which holds.
2. **Laplace, $t\ge4$, $\mu<0.105$.** By §6.5 the saddle equation, for either class family, has a root, and every root satisfies $V>2$. Fix a root. Then $W=V/(5n)$, the representation of §6.1 is exact, and $X=5n/V$.
   - If $V\in[2,12]$, $V$ lies in one of the 500 tiles $[V_{lo},V_{hi}]$. For both class families the tiles are contiguous and the last tile is snapped to end exactly at 12.
     - Classes 0–2: $X\ge1455/V_{hi}$, since $n\ge291$, and $X\ge4V_{lo}/a_5$, since $t\ge4$. So $X\ge X_{\min}$ of the tile.
     - Classes 3–4: $X\ge4105/V_{hi}$, since $n\ge821$, and $X>3V_{lo}/a_5$, since $t\ge4$.
   - If $V\ge11.99$, the tail cell applies; for both class families it is run from $V_T=11.99$, so it overlaps the last tile. It needs $X\ge4V/a_5$ (classes 0–2) or $X\ge X^*=216.3$ together with $V\le a_5X/3$ (classes 3–4), and these follow from the same two inequalities.

   Every certificate term is non-increasing in $X$, so the tile's certified inequality holds at this $X$. The thresholds are $n\ge291$ (classes 0–2) and $n\ge821$ (classes 3–4), both below 980.
3. **Centre, $\mu\ge0.105$.** Since $m\le5n^2$, $\mu\le\tfrac12$. By §7.1, $\mu\mapsto L$ is a bijection from $[0.105,\tfrac12]$ onto $[0,L(0.105)]\subset[0,2.09]$ (since $L(0.105)=2.0814$), and every $L$-box in $[0,2.1]$ is certified at $N_0=821$ for all $n\ge821$. The integration circle is split into the four near-peak windows $|y_h|\le1.1$ and their complement. The complement is covered by (E), (E0), (N) and (F), with the handovers listed in §7.2 and the closure remark for the float endpoints, and Lemma 7.3 holds there. The per-box combination of §7.3 then gives the sign. The threshold is $n\ge821$.

So all thresholds (291, 821 and 2000 for the switch between exact and analytic in the band) are satisfied or bridged for $n\ge980$. The boundaries $t=4$ and $\mu=0.105$ are assigned to the Laplace region and the centre respectively. Each certificate includes its closed boundary: the Laplace certificates use $t\ge4$, and the centre uses $\mu\ge0.105$ through $\mu(2.1)<0.105$. So nothing is left out. ∎

*Numerical cross-check (not part of the proof).* PROOF.md §9 records an audit of 368056 points with $n$ up to $10^6$, which found 0 uncovered points for classes 3 and 4. That audit used the older $v$-based map. The coverage argument above does not depend on it.

---

## 9. Main theorem

**Theorem 9.1.** For every $n\ge1$ and every $m$, $c_n(m)\ge0$ if $5\mid m$ and $c_n(m)\le0$ otherwise. That is, the Third Borwein conjecture holds.

*Proof.* By Lemma 1.1 it suffices to take $0\le m\le5n^2$.

- For $n\le979$ the result is Theorem 4.1.
- For $n\ge980$ it follows from Theorem 8.1, together with Conclusions 5.4, 5.5, §6.3, §6.4 and Conclusion 7.4.

The analytic inputs are Theorem 3.1, the elementary lemmas proved above, and the standard facts of Appendix C. ∎

**Status.** Every region is certified, every certificate has a log in the directory (Appendix A), and nothing is pending. The audit record is in Appendix B:

- The band (classes 3–4), the Laplace region (classes 3–4), the centre near-peak and combination step, and the off-peak lemma have each passed two independent adversarial audits.
- The classes 0–2 band and Laplace certificates have two independent adversarial audits (GAP_AUDIT012.md, AUDIT2_CLASS012.md).
- §3 was refereed independently.
- An outside end-to-end read found one coverage gap, the $V=12$ sliver for Laplace classes 0–2. It is now fixed.

---
## Appendix A. Certificate index

All scripts are in `/home/claude/third-borwein` and use python-flint (FLINT/Arb). Runtimes are those printed in the logs or recorded by the audits, on the 2-core project box; "not recorded" means that neither the log nor an audit states one. "Must print" quotes the line of the log that constitutes the certificate, checked against the log in the directory.

| # | certifies | command(s) | log | must print | runtime |
|---|---|---|---|---|---|
| A.1 | exact, $1\le n\le150$ (§4) | `python3 centre3_exact.py 1 150` | `exact_1_150.log` | `DONE n=1..150 total violations 0, 4s` | 4 s |
| A.2 | exact, $150\le n\le291$ | `python3 centre3_exact.py 150 291` | `exact_150_291.log` | `DONE n=150..291 total violations 0, 34s` | 34 s |
| A.3 | exact, $291\le n\le410$ | `python3 centre3_exact.py 291 410` | `centre3_exact_291_410.log` | `DONE n=291..410 total violations 0, 114s` | 114 s |
| A.4 | exact, $410\le n\le580$ | `python3 centre3_exact.py 410 580` | `centre3_exact_410_580.log` | `DONE n=410..580 total violations 0, 419s` | 419 s |
| A.5 | exact, $580\le n\le820$ | `python3 centre3_exact.py 580 820` | `centre3_exact_580_820.log` | `DONE n=580..820 total violations 0, 1357s` | 1357 s |
| A.6 | exact, $820\le n\le920$ | `python3 centre3_exact.py 820 920` | `centre3_exact_820_920.log` | `DONE n=820..920 total violations 0, 2825s` | 2825 s |
| A.7 | exact, $920\le n\le979$ | `python3 centre3_exact.py 920 1000` | `centre3_exact_920_1000.log` | 60 lines `n=920` … `n=979 wrong_sign=0 [] …`, then `Killed` (out of memory; about 6 GB are needed for NB = 1000) | 3786 s to $n=979$ |
| A.8 | band, exact, $291\le n\le2000$, $m<4N$, all classes (§5) | `python3 bandfix_exact.py 291 2000` | `bandfix_exact.log` | `n in [291,2000]: sign failures (any class, m < 4N): 0;  max \|b2+b3\|/\|b1\| (classes 3,4, 2N<=m<4N) = 0.0001 at n=291, m=5819 (t=m/N=3.997)` | 192 s (log); 50 s in the audit rerun |
| A.9 | one-part block $g(t)$, $G(T)$, $t\le39999$; pickle check (Lemma 5.2) | `python3 bandfix_coeffs.py` | `bandfix_coeffs.log` | `g5_coeffs.pkl agrees with independent recomputation for i <= 200000` and `ONE-PART BLOCK (exact): g(t) < 0 for all t <= 39999 except t = 1 (g=1), t = 5 (g=0); G(T) <= 0 for all T <= 39999` | under 1 min |
| A.10 | Lemma L1 and $\Psi$, classes 3–4, $N\ge10001$ (§5.1) | `python3 bandfix_analytic.py 10001` | `bandfix_analytic.log` | `L1: sup_{X>=1000} 3 Rgrp(X+2)/P5'(X) <= [0.017972 +/- 1.12e-7]`; `N0=10001: D<N and Psi decreasing: True;  Psi(N0) <= [0.19769 +/- 4.78e-6]`; `BAND (analytic): (\|b2\|+\|b3\|)/\|g(T)\| < 0.1977 for all N >= 10001, 2N <= m < 4N, classes 3,4` | under 1 min |
| A.11 | band, classes 0–2, $n>2000$ (§5.2) | `python3 cls012_band.py` | `cls012_band.log` | `PART A: … sign failures 0 []  band failures 0 []  worst corr/\|s(m)\| = 5.058e-10`; `PART B: … Q <= [5.1390e-20 +/- 4.88e-25]` | not recorded |
| A.12 | Laplace classes 0–2, tiles $V\in[2,12]$ (§6.3) | `python3 cls012_lap.py 2.0 12.0 0.02` | `cls012_lap_rerun.log` (line 1) | `V in [2.0,12.0] width 0.02: 500 tiles, 0 failures; min certified margin sgn(a)M - 4 delta - E2 by class: 0: 3.0700, 1: 1.6879, 2: 0.8339` | about 8 s |
| A.13 | Laplace classes 0–2, tail $V\ge11.99$ | `python3 cls012_lap.py tail 11.99` | `cls012_lap_rerun.log` (line 2) | `tail cell V >= 11.99 : ({0: (3.0950…, …), 1: (1.7130…, …), 2: (0.8589…, …)}, …)` | seconds |
| A.14 | saddle coverage, $A(V)$ (§6.5) | `python3 cls012_cov.py` | `cls012_cov.log` | `… 5858 tiles, min A = [0.04424… ]` and `A(V) on (0,0.001] >=  [0.19979… ]` | not recorded |
| A.15 | Laplace classes 3–4, tiles and tail (§6.4) | `python3 lap34_cert.py 2.0 12.0 0.02` and `python3 lap34_cert.py tail 11.99` | `lap34_cert.log` | `V in [2.0,12.0] width 0.02: 500 tiles, 0 failures`; class 3 and class 4 minimum margins `0.8363` and `0.8361` at `(11.979999999999832, 12.0, …)`; `tail cell V >= 11.99 : ({3: (0.8145634857582285, {'epsg': 1.7873302435252118e-05, …` | 76 s for the tiles (audit rerun); tail seconds |
| A.16 | centre saddle endpoint $\mu(2.1)<0.105$ (§7.1) | `python3 centre3_saddle_end.py` | `centre3_saddle_end.log` | `p'(2.1) in [-0.208 +/- 6.04e-4]  mu(2.1) = -p'/2 in [0.104 +/- 3.02e-4]  < 0.105: True` | not recorded |
| A.17 | centre near-peak and combination (§§7.1, 7.3) | `python3 centre3_comb.py 821 210 64 0 105` and `python3 centre3_comb.py 821 210 64 105 210`, outputs concatenated | `centre3_comb_821.log` | on every box `ALL_OK`; after each half `COMBINATION ASSERTED (slack > OFF, arb) on every box; N0=821  boxes=210  all classes certified on every box: True` (twice in all); minimum `minslack − OFF63` = 0.018703 − 0.008660 = 0.010043 at `L=2.090` | not recorded for this run (the earlier single run took about 11 min) |
| A.18 | off-peak (E) (§7.2) | `python3 offpeak_em.py 821 6.0 32 > offpeak_em_821.log` | `offpeak_em_821.log` | `N0=821 YE=6.00: piece (E) CERTIFIED for all n >= N0` (worst V −2.18) | 57 s |
| A.19 | off-peak (E0) | `python3 offpeak_em0.py 821 20 32 > offpeak_em0_821.log` | `offpeak_em0_821.log` | `N0=821 YE0=20.0: piece (E0) CERTIFIED for all n >= N0  boxes=4945 worst V(N0)=-0.23 …` | 6 s |
| A.20 | off-peak (N) | `python3 offpeak_tb.py 821 0.25 > offpeak_tb_821.log` | `offpeak_tb_821.log` | `N0=821: piece (N) CERTIFIED for all n >= N0` (worst V −1.8) | 2 s |
| A.21 | off-peak (F), Arb | `for p in 0 … 7; do python3 offpeak_chunk_arb.py 821 $p -j 2 > arb_F_p$p.log; done` (piece 0 and piece 5 completed in the additional part-runs listed in §7.2), then `python3 offpeak_chunk_arb.py summarize <the 13 logs of §7.2> > arb_F_summary_full.log` | `arb_F_summary_full.log` | `(F) arb: CERTIFIED for all n >= N0 (all rects OK)`; first sorted margin line `47.91 piece=7 rect=26 T=20 …` | not recorded as a total; a full re-run is estimated at about 1.5 h on 2 cores (AUDIT2_OFFPEAK.md) |
| A.22 | off-peak (F), float-interval (superseded by A.21) | `python3 offpeak_chunk.py 821 all > offpeak_chunk_821.log` | `offpeak_chunk_821.log` | `N0=821 (F) all: CERTIFIED for all n >= N0` | 470 s |

Do not cite `arb_F_summary.log`: it is an incomplete earlier summary and ends with `(F) arb: INCOMPLETE`. `centre3_peak_821.log` (`python3 centre3_peak.py 821 210 64`) contains the near-peak part alone. It is reproduced inside A.17.

**Supporting checks, not certificates.**

- `cls012_em_check.py` and `cls012_asym_check.py` (classes 0–2).
- `lap34_exact.py n m…` (multi-modular CRT) then `lap34_validate.py file` → `lap34_validate*.log`. Every line must end in `OK` with `sign(c)=-`. Also `lap34_emchk.py`.
- `offpeak_crosscheck.py` and `audit_offpeak_spot*.py`. Shared off-peak code is in `offpeak_common.py`.
- §3: `g5_coef.py` computes $s(i)$ for $i\le200000$. The constants of §3.9 were computed at 60 digits (SEC3_PROOF.md) and at 50 digits (referee). `cusp10.py` and `cusp15b.py` concern the cusp-by-cusp vanishing, which is not used.

---

## Appendix B. Verification and audit record

### B.1 What was checked independently

**First round.**

| component | reviewer's method | verdict |
|---|---|---|
| §3 (SEC3_PROOF.md) | REFEREE_SEC3.md: independent re-derivation (without reading GAP_CITATION.md); own mpmath code at 50 digits; 10 cusps × 2 points for (3.2)/(3.3); quadrature for Lemma 3.7; exhaustive test $\nu\le5000$ plus 1509 values up to 200000 | correct; three presentational gaps, filled in §3 above |
| Band classes 0–2, Laplace classes 0–2, coverage | GAP_AUDIT012.md: re-derivation, NaN-guarded rerun of all tiles, 45 independent exact points | sound; two latent flaws fixed |
| Band classes 3–4 | The earlier band audit (items 19–21), whose findings GAP_BAND.md repairs | defects fixed by the repair (B.2) |
| Laplace classes 3–4 | AUDIT_LAP34.md: re-derivation, 6000-point test of $\varepsilon_g$, NaN-guarded rerun, 36 exact points computed with independent code | sound after two fixes |
| Near-peak and combination | AUDIT_CENTRE.md: re-derivation of the representation, EM constants, uniformity in $n$, saddle and exact-run logs; exact $S_{300}$ comparison | sound after replacing the combination step |
| Off-peak lemma | AUDIT_OFFPEAK.md: coverage in $\theta$, uniformity in $n$, Lemma M, chunk and maximum principle, 312 exact spot points; Arb version of (F) on 7 rectangles | sound after cosmetic fixes |

**Second round (independent audits; each auditor wrote their findings before reading the first-round report).**

| component | what was checked | findings | verdict |
|---|---|---|---|
| Band classes 3–4 (AUDIT_BAND.md) | exact run re-run; 16 values of $n$ recomputed by a direct product that avoids $E_n$ and the pickle; pair and triple counts brute-forced; Bessel directions; L1 and $\Psi$ recomputed in mpmath ($\Psi(10001)=0.19714$) | two cosmetic misstatements: $\Psi$ decreasing from $N=257$, not "about 230"; the §5.1 ratio is quoted at $n=291$, which is inside the run's range 291–2000 | sound |
| Laplace classes 3–4 (AUDIT2_LAP34.md) | re-derivation from scratch; complete term inventory, including the Trho constants against Theorem 3.1; monotonicity in $X$ at 6 levels; a third exact method (Euler identities mod $p$, CRT), 31 points up to $n=5000$; NaN handling and ball coverage | no hole; ulp-level cosmetic items (`e0` rounds to nearest); stale log | sound |
| Centre (AUDIT2_CENTRE.md) | Cauchy decomposition and $\mathcal N$ re-derived; saddle map, with $\mu(2.1)=0.103772$; thin box recomputed in independent Arb code; monotonicity in $n$ at $n$ up to $10^5$; quadrature validated to $4\cdot10^{-11}$ against exact $S_{300}$, then 120 points at $n=821$–900 for $\mu$ up to 0.5 | D-a: the comparison slack > OFF was printed but not asserted; D-b: the binding box $[2.09,2.10]$ is not needed, and the effective margin is 0.0120 | sound |
| Off-peak (AUDIT2_OFFPEAK.md) | EM, (E0), triangle bound, chunk lemma, subharmonicity and Lemma M re-derived; certified (N) pieces compared with the exact $j$-series at 70 points; certificates re-run at $N_0=2000$, $10^4$, $10^5$; every (F) rectangle of the summary matched bit for bit; 11 pre-fix (F) rectangles re-run; 5,570 adversarial Arb points | N1, N2: float slivers at $\|y\|=1.1$ and $L=0.1$; N3: the constant 6.3 used as a float, $1.8\cdot10^{-16}$ low, in (E), (E0), (N-a) | sound |

**Outside referee (REFEREE_FULL.md).** The referee read the whole document end to end, treating the certificates as black boxes, and spot-checked 13 logs and two scripts. The verdict was "correct with gaps to fill".

- **M1 (major).** For Laplace classes 0–2 the float-stepped tiling ended at 11.999999999999831, while the tail cell started at 12, so $V\in(11.999999999999831,12)$ was covered by nothing. This is the same float-drift defect as in lap34, but it had not been fixed for classes 0–2, and no audit had caught it.
- **m1–m9 (presentation).** Stale text left beside the addendum; the wrong (F) summary log named; $L\in(0,2.1]$ where $[0,2.1]$ is certified; $\varepsilon_g(12)$ where the tail uses $\varepsilon_g(11.99)$; three certificates without logs; the variable of the positive power series in Lemma 5.1; the ratio quoted with more digits than the log prints; steps only listed (the Rankin bound for $\rho$, $\int|\Phi_c''|$, (N) Part II, and the chunk lemma for $T\nmid n$); "exactly one" in Theorem 8.1.

### B.2 Defects found and fixed in the components used

- **§3.**
  - The Farey order must be $N=K(\nu)$, not $K(\nu)-1$; otherwise $k=5$ is missing at $\nu=1$.
  - $E_{\rm con}$ had been rounded down (27.4635); it is now rounded up to 27.4636.
  - $K_5=2.8549954$ was quoted truncated; the rounded-up value is 2.8549955.
  - The source manuscript contained no written $G_p$ proof, so the proof was transplanted and written out in full.
- **Band.**
  - The old majorant used a Bessel envelope in the wrong direction, and its crossing was near $n\approx480$ rather than 322.
  - `band_ext_cert.py` undercounted pairs by 4× and triples by 10×.
  - `band_mono.py` dropped NaN results silently.

  All three are superseded by the exact run to 2000 and by the majorant $\Psi$ (GAP_BAND.md). AUDIT_BAND.md corrected the threshold for "$\Psi$ decreasing" to $N\ge257$.
- **cls012_lap.**
  - The pass test `min(...) <= 0` would let a NaN through; it is now `not all(x > 0)`.
  - The decay assert used $\kappa-a_5/5$ where the Rankin term needs $\kappa-a_5/4$.
  - (Referee M1) Float drift left $V\in(11.999999999999831,12)$ uncovered. It was fixed in two independent ways. The last tile is now snapped to 12 (`if Vb-v2<1e-9: v2=Vb`), and the tail cell is run from 11.99, which overlaps the last tile. `cls012_lap_rerun.log` holds the rerun: tile margins unchanged (3.0700 / 1.6879 / 0.8339), tail margins 3.0950 / 1.7130 / 0.8589.
- **lap34.**
  - (A) Float drift left $V\in(11.99999999999983,12)$ uncovered; fixed by snapping the last tile and running the tail cell from 11.99.
  - (B) NaN pass-through in the supremum computations and in the pass test; fixed with `_fin` guards and `not all(x > 0)`. The rerun gives 0 non-finite values.
  - `lap34_cert.log` now comes from the fixed script.
- **Centre combination.**
  - (AUDIT_CENTRE.md D1) The scalar comparison $|\mathrm{OFF}|\le0.0187\le$ slack is false (0.018791 against 0.018703). It is replaced by the per-box comparison.
  - (AUDIT2_CENTRE.md D-a) The per-box comparison is now asserted in Arb inside `centre3_comb.py`, with `arb('-6.3')`, and the log was regenerated. The worst margin is 0.010043, and 0.0120 on the boxes actually needed.
- **Off-peak lemma.**
  - D1: the sliver $(\mathrm{fl}(\pi),\pi]$ was uncovered. It is fixed in the Arb version and also closed by the closure remark.
  - D2: the $h$/sign reduction for (N) was unstated. It is now stated in §7.2.
  - D3: a float key was used in the assembly. Fixed in the Arb version.
  - D4: the float checker had been validated only on zero-width cells. Superseded by Arb.
  - The Arb (F) run stalled because of an early exit in `arb_ok` that could only reject cells. The exit was removed and the run completed (ARB_F_RUN.md).
  - N1, N2: float slivers at $|y|=1.1$ and $L=0.1$. They are closed by the closure remark in §7.2.
  - N3: the constant 6.3 is now `arb('6.3')` everywhere. (E), (E0) and (N) were re-run: worst values −2.18, −0.23 and −1.8 (after also extending the lowest $L$-box to start at 0.0999).
- **Earlier centre audit (GAP_CLASS012 C.2), all fixed in piece (A).**
  - The strip next to each peak was bounded by nothing.
  - The tail was normalised with a lower curvature, 7.4× too small.
  - The tail lemma was unproved (now replaced by Lemma 7.3).
  - The EM constant 5/12·TV was invalid.
  - The Taylor constant was evaluated at the centre instead of over the window.
  - The arctan law was used without its error.
  - $L<0.05$ was not covered.
- **Documentation.**
  - The exact run for $n\le290$ is now logged.
  - `centre3_saddle_end.log`, `bandfix_analytic.log` and `bandfix_coeffs.log` now exist.
  - All the referee's minor points (m1–m9) are addressed in this version.

**Final round (19 September).** A second independent audit of the classes 0–2 certificates (AUDIT2_CLASS012.md) re-derived the representation, inventoried every term, tested monotonicity in $X$ up to $1000X_{\min}$ (including the snapped tile and the tail cell), and validated 76 hard points by an exact CRT method over 82 primes: all signs correct, true error 15–260× below the certified bound. It found no mathematical error; one latent float endpoint in `X4sup` (not triggered at any tile or grid point) and stale text in GAP_CLASS012.md. Separately, the $L$-grid of off-peak pieces (E), (E0), (N) was extended to start at 0.0999, so the exact $L=0.1$ lies in a certified box; all three re-ran certified ((E0) worst $-0.23$).

### B.3 Withdrawn components (not used anywhere above)

| component | why withdrawn | source |
|---|---|---|
| `B_uniform5.py` (primary 983-tile Laplace certificate, all classes) | NaN zone-2 cells replaced by $e^{-10^6}$, silently dropping 7168 cells; with a NaN guard, 63 of 73 tiles fail | GAP_PRIMARY.md |
| `mono_cert.py`, `zone2_cert.py`, `tail3_cell.py`, `tail3_K0.py`, `derived_tile*.py`, `repr_tile*.py`, `thin_cert.py`, the $K_0$ thin-factor lemma, $M=25$ / `Mfix*.py`, $t_{\rm eff}$ (old classes 3–4 Laplace certificate) | Centred on the kernel-only saddle $m=a_5/W^2$ and treated the $O(N)$ thin factor $e^{-5nu}$ and the twisted-product factor as amplitudes. Against the exact line integral its main term was off by factors down to $e^{-13}$, while certifying $B\le0.06$, so "$B<1$" proved nothing. The zone-2 range also passed through the fifth roots (GAP_ZONE2.md), and the main-term normalisation had been wrong earlier (PROOF.md §7). Superseded by §6.4. | PROOF.md §15, GAP_SMALL.md, GAP_MONO.md, GAP_ZONE2.md, GAP_TAIL.md |
| The folded identity with $D_0=1.12601$ (`repr_D0.py`, GAP_REPR.md) | Correct, but not needed: §6.1 is exact without folding | GAP_LAP34.md §2 |
| `centre_hybrid*.py`, `centre_final.py`, `centre_sound_full.py`, `cls012_centre.py` (old centre certificate) | The holes listed in GAP_CLASS012 C.2. Class 4 fails at the handover once the tail normalisation is corrected (0.053 + 0.082 > 0.126). Superseded by `centre3_*` and `offpeak_*`. | GAP_CLASS012.md §C, GAP_CENTRE.md |
| `band_ext_cert.py`, `band_mono.py`, `band_majorant4.py`, `band_tail5.py` | Wrong counts, sampled monotonicity, wrong-direction envelope | GAP_BAND.md |
| `tmin3.py` $t_{\min}$ table; extension to $v=1.4$ | Wrong normalisation | PROOF.md §7 |

---

## Appendix C. Standard facts relied on

- **C.1 Dedekind η / partition-function transformation.** This is Fact 3.5, in Rademacher's normalisation with $hH\equiv-1\pmod k$ (Apostol, *Modular Functions and Dirichlet Series in Number Theory*, Thm 5.1). It was confirmed numerically to $10^{-49}$ at 10 cusps.
- **C.2 Farey dissection and the Rademacher path via Ford circles.** This is Lemma 3.2 (Apostol Ch. 5 §§5.4–5.7; Rademacher, *Topics in Analytic Number Theory*, Ch. 14).
- **C.3 Farey neighbours.** Consecutive fractions of order $N$ satisfy $k+k'\ge N+1$ (Hardy–Wright, Thm 30).
- **C.4 Bessel contour integral.** DLMF 10.9.19 and Watson §6.22. The passage from the loop to the vertical line is proved in §3.5 [gap (a)].
- **C.5 Bessel facts.**
  - $I_\nu(y)>0$ for $y>0$, and $I_\nu(y)$ is decreasing in $\nu\ge0$ for fixed $y>0$ (Cochran 1967; Watson).
  - The closed forms of $I_{1/2}$, $I_{3/2}$ and $I_{5/2}$.
  - $P_5(i)=\sqrt{a_5/x}\,I_1(2\sqrt{a_5x})$ is the inverse Laplace transform of $e^{a_5/u}-1$, where $x=i-1/6$; this follows from the series of $I_1$.

  These are used in §5 (L1, $\Psi$, cls012 Part B) and §6.1.
- **C.6 Dilogarithm identities.** $\sum_{c=0}^4\mathrm{Li}_2(\zeta^cx)=\mathrm{Li}_2(x^5)/5$, with the principal branch for $|x|<1$. Also the Euler pentagonal theorem (§2).
- **C.7 Complex analysis.** Cauchy's theorem, Fubini/Tonelli, the maximum principle for subharmonic functions, and Hadamard's three-circles theorem. The last two are used in piece (F).
- **C.8 Correctness of the arithmetic software.**
  - Arb, through python-flint: every ball returned contains the true value.
  - FLINT `fmpz_poly`: exact integer arithmetic.
  - IEEE-754 correct rounding of $+,-,\times,\div$: this was relied on only in the float-interval version of piece (F), which the Arb certificate supersedes.
  - mpmath is used only in checks, never in a certificate inequality. The exception is the §3.9 constants, which were recomputed independently by the referee.


