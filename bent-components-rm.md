# The non-bent components of a quadratic function form a Reed–Muller codeword

*Bent components · quadratic functions · 2026-11-11 · answers an open problem*

## In plain words

> A component of a vectorial function is "bent" when it is as far from linear as a function can be. For quadratic functions, whether a component is bent is decided by one polynomial of low degree, so the components that are not bent always form a set of a very particular kind: the support of a word in a Reed–Muller code. Results from 1969 and 1970 on which sizes such sets can have then say how few non-bent components a function can have. That settles a question a 2026 paper left open (no function meets its bound), improves the paper's bound, and rules out in a few lines two cases in dimension 8 that the [companion note](apn8-not-sixteen.md) settled by long computer searches.

## Setting

Let $`n`$ be even and $`F:\mathbb F_2^n\to\mathbb F_2^m`$ quadratic. For $`b\in\mathbb F_2^m`$ the component $`b\cdot F`$ has an alternating matrix $`M_b=\sum_k b_kM_k`$, where $`M_k`$ belongs to the $`k`$-th coordinate of $`F`$. Let $`N_F=\{b\neq0: b\cdot F \text{ not bent}\}`$ (affine components count as not bent), $`r=n/2`$, and $`\mathrm{RM}(r,m)`$ the Reed–Muller code of Boolean functions on $`\mathbb F_2^m`$ of degree at most $`r`$, with minimum distance $`d=2^{m-r}`$. $`F`$ is *MNBC* (maximum number of bent components) when $`|N_F|=2^{m-r}-1`$, the least possible [1, after Thm III.4].

**Lemma 1** (known [2]). A quadratic Boolean function is bent exactly when its alternating form is nondegenerate, so $`g(b):=[\,b\cdot F\text{ bent}\,]=\det M_b=\mathrm{Pf}(M_b)^2=\mathrm{Pf}(M_b)`$ over $`\mathbb F_2`$. The Pfaffian has degree $`n/2`$ in the entries of $`M_b`$, which are linear in $`b`$, so $`\deg g\le n/2`$, and $`g(0)=0`$.

Carlet and Thornburgh give this formula [2, Remark 6.2], and derive $`\deg g\le n/2`$ for every plateaued $`F`$ from [1, Prop III.1], the fact that a subspace of dimension above $`n/2`$ holds an odd number of non-bent components [2, Lemma 5.3]; the two statements are equivalent. What follows adds one step: the low-weight structure of Reed–Muller codes.

## The bounds

> **Theorem 2.** Let $`n\ge4`$ be even, $`m\ge n/2+1`$, and $`F:\mathbb F_2^n\to\mathbb F_2^m`$ quadratic. Either $`F`$ is MNBC ($`N_F\cup\{0\}`$ is a subspace of dimension $`m-n/2`$), or $`|N_F|\ge 3\cdot2^{m-n/2-1}-1`$.

*Proof.* $`h:=1+g`$ lies in $`\mathrm{RM}(r,m)`$, $`h(0)=1`$, and its weight is $`|N_F|+1`$, since the zeros of $`g`$ are $`0`$ and $`N_F`$. For $`r\ge2`$, a codeword of $`\mathrm{RM}(r,m)`$ whose weight lies in $`[d,2d)`$ has weight $`2d-2^i`$ for some $`i`$ (Berlekamp and Sloane [3]; Kasami and Tokura [4] found exactly which $`i`$ occur and which codewords have them; see [5, Ch. 15, Thm 11] and [6, Thm 2]). A weight above $`d`$ of that form is at least $`2d-2^{m-r-1}=1.5d`$. Codewords of weight $`d`$ are the indicators of $`(m-r)`$-dimensional affine flats [5, Ch. 13]; as $`h(0)=1`$, the flat is a linear subspace and $`F`$ is MNBC. Otherwise $`|N_F|+1\ge1.5d=3\cdot2^{m-r-1}`$. ∎

The theorem says more than parity only from $`m\ge n/2+3`$: $`|N_F|`$ is odd [1, Prop III.1], which already gives the bound for $`m=n/2+1`$ and $`m=n/2+2`$. From $`m=n/2+3`$ on it is new. Bénéteau, Goluboff, Kölsch and Vaghasiya [1, Thm III.7] proved, for $`n\le 2m-6`$ (that is, $`m\ge n/2+3`$), that a non-MNBC $`F`$ has $`|N_F|\ge 2^{m-r}+2^{m-r-2}+2^{m-r-3}-1=11\cdot2^{m-r-3}-1`$; Theorem 2 gives $`12\cdot2^{m-r-3}-1`$.

> **Corollary 3** (answers [1, Open Problem VI.4]: no). For even $`n\le 2m-6`$, no quadratic $`F:\mathbb F_2^n\to\mathbb F_2^m`$ has exactly $`2^{m-n/2}+2^{m-n/2-2}+2^{m-n/2-3}-1`$ non-trivial non-bent components.

For $`n\ge4`$ this value lies strictly between the MNBC value and Theorem 2's bound. For $`n=2`$, $`h=1+g`$ is affine, so $`|N_F|`$ is $`2^{m-1}-1`$ or $`2^m-1`$, and the value is neither.

*A second proof, with the paper's own tools* ($`n\le2m-6`$). Equality makes $`N_F`$ a smallest nontrivial blocking set with respect to $`(n/2)`$-spaces of $`\mathrm{PG}(m-1,2)`$ [1, proof of Thm III.7]. By Govaerts and Storme [7, Thm 1.4(3)] it is then a cone: a vertex $`V`$ of vector dimension $`m-n/2-3`$ and, in a 4-dimensional $`U`$ with $`U\cap V=0`$, the base $`B=l\cup l_1\cup l_2\cup l_3`$, three concurrent lines through a point $`P`$ that span $`U`$ and a line $`l`$ skew to them (the complement of an elliptic quadric in $`\mathrm{PG}(3,2)`$). For $`Q\in l`$ the line $`PQ`$ meets $`B`$ exactly in $`P`$ and $`Q`$. Take $`C`$ complementary to $`V\oplus U`$ (dimension $`n/2-1`$) and $`W=\langle P,Q\rangle\oplus C`$, of dimension $`n/2+1`$. Since $`N_F\cup\{0\}=V+(B\cup\{0\})\subseteq V\oplus U`$ and $`W\cap(V\oplus U)=\langle P,Q\rangle`$, $`W`$ contains exactly two non-bent components, but [1, Prop III.1] says the number is odd. ∎ (The wording of [7, Thm 1.4] used here is that of a 2004 draft of [7], consistent with [1]'s quotation of the published theorem.)

> **Theorem 4.** Let $`n\ge4`$ be even and $`F:\mathbb F_2^n\to\mathbb F_2^n`$ quadratic APN. Then $`|N_F|\ge 2^{n/2+1}-3`$: $`F`$ has at most $`2^n-2^{n/2+1}+2`$ bent components.

*Proof.* $`|N_F|\equiv1\pmod 4`$ [1, Thm III.8; 8, Lemma 6.4]. (Directly: the fourth moment of the Walsh transform makes $`\sum_{b\ne0}2^{k_b}=2(2^n-1)`$, where $`2^{(n+k_b)/2}`$ is the amplitude of $`b\cdot F`$ ($`2^{k_b}`$ is the size of its radical), so $`\sum_{b\in N_F}(2^{k_b}-1)=2^n-1`$; each $`k_b\ge2`$ is even, so each term is $`3`$ modulo $`4`$, as is $`2^n-1`$.) So the weight of $`h`$ is $`\equiv2\pmod4`$. Here $`d=2^{n/2}`$, and if the weight is below $`2d`$ it is $`2d-2^i`$, which is $`\equiv2\pmod4`$ only for $`i=1`$: weight $`2^{n/2+1}-2`$. ∎

At equality, $`N_F\cup\{0\}`$ is, by Kasami and Tokura's description of these codewords [4; 6, Thm 2], the symmetric difference of two $`(n/2)`$-dimensional affine flats that meet in one point. Equality holds for $`n=4`$: $`x^3`$ on $`\mathbb F_{16}`$ has 5 non-bent components.

**What Theorem 4 adds.** [1] has three lower bounds for quadratic APN functions and compares them [1, Prop IV.4], by the dimension $`k`$ of the largest radical (linearity $`2^{(n+k)/2}`$). For $`k\neq n/2`$ their bounds (A) and (C), or simply counting points ($`k>n/2`$ forces $`|N_F|\ge2^k+1`$), already exceed $`2^{n/2+1}`$. Only for $`k=n/2`$, which needs $`4\mid n`$ since $`k`$ is even, was Thm III.9 their best: $`2^{n/2}+2^{n/2-2}+2^{n/2-3}-1`$, now $`2^{n/2+1}-3`$ (21 → 29 for $`n=8`$, 87 → 125 for $`n=12`$). So for $`4\mid n`$ Theorem 4 improves the lower bound on $`|N_F|`$ over all quadratic APN functions; for $`n\equiv2\pmod4`$ it adds nothing.

**Dimension 8.** The distributions $`[0^{238-4i},2^{5i},4^{17-i}]`$ of [1, Thm IV.5] have $`|N_F|=4i+17`$, so Theorem 4 leaves $`i\ge3`$; the remaining type $`[0^{190},2^{64},6^1]`$ has no rank-4 components. So no quadratic APN function on $`\mathbb F_2^8`$ has more than 14 rank-4 components. That is the $`i=1`$ case of [1, Open Problem VI.1] (first settled by Carlet and Thornburgh [2, Cor 5.13]) and the $`i=2`$ case (until now settled only by the searches in the [companion note](apn8-not-sixteen.md)), each in a few lines. It stops at $`i=3`$: weight 30 is a Kasami–Tokura weight. Weights alone don't go further either: $`\mathrm{RM}(4,8)`$ has words of weight 34 ($`i=4`$).

## More from the same argument

- **All plateaued functions.** Lemma 1 holds for every plateaued $`F`$ with $`n`$ even [2, Lemma 5.3], and the proofs above use nothing else, so Theorem 2 holds for plateaued $`F`$ and Theorem 4 for plateaued APN $`F`$ (the mod-4 argument above works for them unchanged).
- **The whole low range.** If $`|N_F|<2^{m-n/2+1}-1`$, then $`|N_F|=2^{m-n/2+1}-2^i-1`$ for an $`i`$ allowed by [4] (see [6, Thm 2]).
- **Divisibility.** By McEliece's theorem [5, Ch. 15], $`|N_F|+1`$ is divisible by $`2^{\lceil 2m/n\rceil-1}`$, which says more than parity once $`m>n`$.

## How tight

For $`n=4`$, Theorem 2 is the true minimum: over all subspaces of $`\mathrm{Alt}(\mathbb F_2^4)`$, the smallest non-MNBC values of $`|N_F|`$ for $`m=4,\dots,7`$ are 5, 11, 23 and 47, exactly the bound ([1]'s gives 10, 21, 43 for $`m=5,6,7`$). Appending a zero coordinate doubles both $`|N_F|+1`$ and $`d`$, so a function on the bound stays on it: $`x^3`$ on $`\mathbb F_{16}`$ gives tightness for $`n=4`$ and every $`m\ge4`$, and a quadratic $`(6,5)`$-function with 5 non-bent components (found by search; it is in the code folder) for $`n=6`$ and every $`m\ge5`$. Whether Theorem 2 is tight for $`n\ge8`$ (it would need, for instance, an $`(n,n/2+2)`$-function with exactly 5 non-bent components), and whether Theorem 4 is tight for any $`n\ge6`$, is open.

## Checks

`checks.py` in [code/bent-components-rm](code/bent-components-rm), a few minutes in all: the degree of $`g`$ on Gold functions in dimensions 6, 8 and 10 and on random quadratic $`(n,m)`$-functions (never above $`n/2`$); Kasami and Tokura's weights against an exhaustive enumeration of $`\mathrm{RM}(2,5)`$ and $`\mathrm{RM}(2,6)`$; the table for $`n=4`$; and the $`(6,5)`$ example with its paddings. Before publishing, the argument was reviewed by Claude Fable 5.2, a different model from the author, which re-derived each step, checked the second proof's base on coordinates, found an independent $`(6,5)`$ example, and corrected the comparison for Theorem 4 and the attribution of Lemma 1.

## References

1. S. H. Bénéteau, N. Goluboff, L. Kölsch, D. Vaghasiya. On the Walsh spectra of quadratic APN functions. IEEE Trans. Inf. Theory 72(7):5207–5216, 2026. [arXiv:2510.12008](https://arxiv.org/abs/2510.12008). Prop III.1, Thms III.7–III.9, IV.5–IV.6, Prop IV.4, Open Problems VI.1 and VI.4.
2. C. Carlet, D. Thornburgh. Resolving a conjecture on quadratic APN functions and a new quadratic $`(n,n)`$-function associated to crooked functions. [arXiv:2608.23888](https://arxiv.org/abs/2608.23888) (v2, September 2026), Lemma 5.3, Corollary 5.13, Remark 6.2.
3. E. R. Berlekamp, N. J. A. Sloane. Restrictions on weight distribution of Reed–Muller codes. Information and Control 14(5):442–456, 1969.
4. T. Kasami, N. Tokura. On the weight structure of Reed–Muller codes. IEEE Trans. Inf. Theory 16(6):752–759, 1970.
5. F. J. MacWilliams, N. J. A. Sloane. The Theory of Error-Correcting Codes. North-Holland, 1977.
6. C. Carlet, P. Solé. The weight spectrum of two families of Reed–Muller codes. Discrete Math., 2023. [arXiv:2301.13497](https://arxiv.org/abs/2301.13497), Theorem 2 (the statement of [4] used here).
7. P. Govaerts, L. Storme. The classification of the smallest nontrivial blocking sets in PG(n,2). J. Combin. Theory Ser. A 113(7):1543–1548, 2006.
8. L. Kölsch, A. Polujan. The combinatorial structure and value distributions of plateaued functions. [arXiv:2410.00611](https://arxiv.org/abs/2410.00611), Lemma 6.4.
