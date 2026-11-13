"""Checks for "The non-bent components of a quadratic function form a Reed-Muller codeword".

Plain Python 3 with numpy. A few minutes in all.
  python checks.py            all checks
  python checks.py 1 3        a selection
  python checks.py 1 --luts apn_8bit.txt   also check the degree on every 8-bit LUT in a file (one LUT per line)

1. The bent indicator g(b) = [b.F bent] has algebraic degree <= n/2 (Lemma 1): Gold functions in dimensions 6, 8, 10,
   random quadratic (n, m)-functions, and optionally a file of 8-bit LUTs.
2. Kasami-Tokura's weights in [d, 2d) against exhaustive enumeration of RM(2,5) and RM(2,6).
3. n = 4: every subspace T of Alt(F_2^4); |N_F| = 2^(m - dim T) |T cap D| - 1, D the degenerate forms with 0.
   Smallest non-MNBC |N_F| for m = 4..7 against Theorem 2's bound.
4. A quadratic (6,5)-function with 5 non-bent components (Theorem 2's bound), padded to (6,6) and (6,7).
5. Theorem 4 on the APN functions of item 1: |N_F| >= 2^(n/2+1) - 3.
6. The (8,6)-function (xy, Tr(alpha x^3), Tr(alpha^2 y^3)) on F_16 x F_16 with 5 non-bent components (Theorem 2's bound),
   padded to (8,7) and (8,8); and every (xy, q1(x), q2(y)) with q1, q2 quadratic forms on F_16: which have 5.
7. Theorem 5: the inequality of its steps 2 and 7 for n = 8, 12, ..., 48; for n = 8, the condition of step 5 by brute
   force over every set T of size >= 12 and every word of RM(2,4); with --luts, the radical partition (P), and (S) and
   (N) on every pair of rank-4 components of each 8-bit APN function in the file.
8. Theorem 6: which weights of Kasami-Tokura-Azumi's list for [2d, 2.5d) are 2 mod 4 when m = 2r (r = 5..12); their
   Table I entries of degree r evaluated directly for r = 5, 6, 7; and the affine stabilizer of
   x1x2x3x4 + x5x6(x7x8 + x3x4), whose orbit is all weight-34 words of RM(4,8) by their count (a few minutes).
"""
import itertools, random, re, sys
import numpy as np


def pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def rank2(rows):
    rows, r = list(rows), 0
    while rows:
        p = rows.pop()
        if p:
            r += 1
            lb = p & -p
            rows = [x ^ p if x & lb else x for x in rows]
    return r


def rows_of(c, n, P):
    rows = [0] * n
    for k, (i, j) in enumerate(P):
        if c >> k & 1:
            rows[i] |= 1 << j
            rows[j] |= 1 << i
    return rows


def bent_indicator(forms, n):
    """forms[k]: the k-th coordinate's alternating form as a bit vector over pairs(n). Returns g on F_2^m."""
    P, m = pairs(n), len(forms)
    g = np.zeros(1 << m, dtype=np.uint8)
    for b in range(1, 1 << m):
        c = 0
        for k in range(m):
            if b >> k & 1:
                c ^= forms[k]
        g[b] = rank2(rows_of(c, n, P)) == n
    return g


def anf_degree(f, m):
    a = f.astype(np.uint8).copy()
    for i in range(m):
        s = 1 << i
        a = a.reshape(-1, 2 * s)
        a[:, s:] ^= a[:, :s]
        a = a.reshape(-1)
    return max((bin(u).count('1') for u in range(1 << m) if a[u]), default=-1)


def forms_of_lut(F, n):
    beta = [[F[(1 << i) ^ (1 << j)] ^ F[1 << i] ^ F[1 << j] ^ F[0] for j in range(n)] for i in range(n)]
    out = []
    for v in range(n):
        c = 0
        for k, (i, j) in enumerate(pairs(n)):
            if beta[i][j] >> v & 1:
                c |= 1 << k
        out.append(c)
    return out


POLY = {6: 0b1000011, 8: 0b100011011, 10: 0b10000001001}


def gold(n, k):
    def mul(a, b):
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a >> n & 1:
                a ^= POLY[n]
        return r
    e = (1 << k) + 1
    out = []
    for x in range(1 << n):
        r, base, t = 1, x, e
        while t:
            if t & 1:
                r = mul(r, base)
            base = mul(base, base)
            t >>= 1
        out.append(r if x else 0)
    return out


def check1(luts):
    print('== 1. degree of the bent indicator (Lemma 1: at most n/2)')
    rng = random.Random(1)
    apn = []
    for n, k in ((6, 1), (8, 1), (8, 3), (10, 1), (10, 3)):
        g = bent_indicator(forms_of_lut(gold(n, k), n), n)
        nf = (1 << n) - 1 - int(g[1:].sum())
        apn.append((n, nf))
        print(f'  x^{(1 << k) + 1} on F_2^{n}: degree {anf_degree(g, n)}, |N_F| = {nf}')
    for n, m in ((6, 6), (6, 9), (8, 5), (8, 10), (10, 7)):
        mx = max(anf_degree(bent_indicator([rng.getrandbits(n * (n - 1) // 2) for _ in range(m)], n), m) for _ in range(20))
        print(f'  20 random quadratic ({n},{m})-functions: max degree {mx} (bound {n // 2})')
    if luts:
        degs, nfs = set(), set()
        for line in open(luts, encoding='utf-8'):
            F = [int(t) for t in re.findall(r'\d+', line)]
            if len(F) != 256:
                continue
            g = bent_indicator(forms_of_lut(F, 8), 8)
            degs.add(anf_degree(g, 8))
            nfs.add(255 - int(g[1:].sum()))
        print(f'  {luts}: degrees {sorted(degs)}, |N_F| values {sorted(nfs)}')
        apn += [(8, v) for v in nfs]
    return apn


def check2():
    print('== 2. Kasami-Tokura in small cases')
    for r, m in ((2, 5), (2, 6)):
        cols = []
        for k in range(r + 1):
            for s in itertools.combinations(range(m), k):
                cols.append(sum(1 << x for x in range(1 << m) if all(x >> i & 1 for i in s)))
        words = [0]
        for c in cols:
            words += [w ^ c for w in words]
        d = 1 << (m - r)
        seen = sorted({bin(w).count('1') for w in words} & set(range(d, 2 * d)))
        mu_max = int(max(min(m - r, r), (m - r + 2) / 2))
        pred = sorted({2 * d - (1 << (m - r + 1 - mu)) for mu in range(1, mu_max + 1)})
        print(f'  RM({r},{m}): {len(words)} codewords; weights in [{d},{2 * d}): {seen}; Kasami-Tokura: {pred}')


def check3():
    print('== 3. n = 4: smallest |N_F| over all quadratic (4,m)-functions')
    P = pairs(4)
    deg = [rank2(rows_of(c, 4, P)) < 4 for c in range(64)]
    by_dim, level = {}, {frozenset([0])}
    for t in range(7):
        by_dim[t] = sorted({sum(deg[x] for x in T) for T in level})
        level = {frozenset(list(T) + [x ^ v for x in T]) for T in level for v in range(1, 64) if v not in T}
    for m in (4, 5, 6, 7):
        vals = {(1 << (m - t)) * c - 1 for t in range(min(m, 6) + 1) for c in by_dim[t]}
        mnbc = (1 << (m - 2)) - 1
        other = min(v for v in vals if v != mnbc)
        old = f'{11 * 2 ** (m - 5) - 1}' if m >= 5 else 'n/a'
        print(f'  m = {m}: MNBC {mnbc}, smallest other {other}, Theorem 2 {3 * 2 ** (m - 3) - 1}, [1, Thm III.7] {old}')


def check4():
    print('== 4. a (6,5)-function on Theorem 2\'s bound, and its paddings')
    forms = [15103, 13438, 18684, 27620, 13103]  # alternating forms on F_2^6 over pairs(6)
    for extra in (0, 1, 2):
        f = forms + [0] * extra
        m = len(f)
        g = bent_indicator(f, 6)
        nb = [b for b in range(1, 1 << m) if not g[b]]
        print(f'  (6,{m}): |N_F| = {len(nb)}, bound 3*2^(m-4)-1 = {3 * 2 ** (m - 4) - 1}, MNBC would be {2 ** (m - 3) - 1}')


def check5(apn):
    print('== 5. Theorem 4 on the APN functions above: |N_F| >= 2^(n/2+1) - 3')
    for n, nf in sorted(set(apn)):
        print(f'  n = {n}: |N_F| = {nf} >= {2 ** (n // 2 + 1) - 3}: {nf >= 2 ** (n // 2 + 1) - 3}')


def check6():
    print("== 6. an (8,6)-function on Theorem 2's bound: (xy, Tr(alpha x^3), Tr(alpha^2 y^3)) on F_16 x F_16, and its paddings")
    def mul(a, b):  # F_16 = F_2[t]/(t^4 + t + 1); t = 2 is primitive
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a & 16:
                a ^= 0b10011
        return r
    def tr(a):
        s, x = 0, a
        for _ in range(4):
            s ^= x
            x = mul(x, x)
        return s
    cube = lambda x: mul(x, mul(x, x))
    al, al2 = 2, 4
    lut = [mul(v & 15, v >> 4) | tr(mul(al, cube(v & 15))) << 4 | tr(mul(al2, cube(v >> 4))) << 5 for v in range(256)]
    forms = forms_of_lut(lut, 8)[:6]
    P8 = pairs(8)
    for extra in (0, 1, 2):
        f = forms + [0] * extra
        m = len(f)
        g = bent_indicator(f, 8)
        nb = [b for b in range(1, 1 << m) if not g[b]]
        print(f'  (8,{m}): |N_F| = {len(nb)}, bound 3*2^(m-5)-1 = {3 * 2 ** (m - 5) - 1}, MNBC would be {2 ** (m - 4) - 1}')
        if extra == 0:
            def form(b):
                c = 0
                for k in range(m):
                    if b >> k & 1:
                        c ^= f[k]
                return c
            print(f'  its non-bent components have ranks {sorted(rank2(rows_of(form(b), 8, P8)) for b in nb)}')
            # component b = (b_xy, c1, c2): b_xy . (xy) = Tr(a xy) with b_xy[i] = Tr(a t^i)
            bvec = lambda a: sum(tr(mul(a, 1 << i)) << i for i in range(4))
            expect = {bvec(0) | 1 << 4, bvec(0) | 1 << 5} | {bvec(a) | 3 << 4 for a in range(16) if cube(a) == cube(al)}
            print(f'  non-bent are exactly Tr(alpha x^3), Tr(alpha^2 y^3) and the three Tr(axy) + Tr(alpha x^3) + Tr(alpha^2 y^3) with a^3 = alpha^3: '
                  f'{set(nb) == expect}')
    # the family (xy, q1(x), q2(y)): q1, q2 range over all 64 x 64 alternating forms on F_2^4 (linear terms don't matter)
    P4 = pairs(4)
    lift = lambda q, off: sum(1 << P8.index((i + off, j + off)) for k, (i, j) in enumerate(P4) if q >> k & 1)
    xy = forms[:4]
    def polar4(fun):
        beta = [[fun((1 << i) ^ (1 << j)) ^ fun(1 << i) ^ fun(1 << j) ^ fun(0) for j in range(4)] for i in range(4)]
        return sum(beta[i][j] << k for k, (i, j) in enumerate(P4))
    tr3 = {polar4(lambda x, l=l: tr(mul(l, cube(x)))): l for l in range(16)}
    cubes = {cube(x) for x in range(1, 16)}
    hits = []
    for q1 in range(64):
        for q2 in range(64):
            g = bent_indicator(xy + [lift(q1, 0), lift(q2, 4)], 8)
            if sum(1 for b in range(1, 64) if not g[b]) == 5:
                hits.append((q1, q2))
    ok = all(q1 in tr3 and q2 in tr3 and tr3[q1] not in cubes and tr3[q2] not in cubes and mul(tr3[q1], tr3[q2]) in cubes
             for q1, q2 in hits)
    want = sum(1 for l in range(1, 16) for u in range(1, 16) if l not in cubes and u not in cubes and mul(l, u) in cubes)
    print(f'  family (xy, q1(x), q2(y)), all 64 x 64 pairs of forms: {len(hits)} give |N_F| = 5; each is Tr(l x^3), Tr(u y^3) with '
          f'l, u not cubes and l*u a cube: {ok}; pairs of that kind: {want}')


def rm_words(m, r):
    """Every word of RM(r, m) as a bitmask over F_2^m (small m only)."""
    mons = [s for s in range(1 << m) if bin(s).count('1') <= r]
    cols = [sum(1 << x for x in range(1 << m) if x & s == s) for s in mons]
    words = [0]
    for c in cols:
        words += [w ^ c for w in words]
    return words


def check7(luts):
    print('== 7. Theorem 5: Theorem 4 is not attained for n = 4t >= 8')
    from math import comb
    print('  steps 2 and 7: a = #radicals of dimension n/2 >= 2M/3 + 11/3 - 16/(3M) (M = 2^(n/2)), |G| <= M + 2 - a,'
          ' and |G| >= dim RM(t-1, 2t) is needed')
    for t in range(2, 13):
        M = 4 ** t
        a = -(-(M * M // 2 + 11 * M // 4 - 4) * 4 // (3 * M))     # ceil((M^2/2 + 11M/4 - 4) / (3M/4))
        need = sum(comb(2 * t, i) for i in range(t))
        gmax = M + 2 - a
        if t <= 4 or gmax >= need:
            print(f'    n = {4 * t}: a >= {a}, |G| <= {gmax}, needed {need}: contradiction {gmax < need}')
        assert gmax < need
    print('    n = 20 .. 48: contradiction True (each checked)')
    # n = 8 directly: the condition of step 5 for every u in T, over all T with |T| >= 12 and all of RM(2,4)
    words = rm_words(4, 2)
    worst = 0
    for size in (12, 13, 14, 15):
        for T in itertools.combinations(range(1, 16), size):
            Tm = sum(1 << u for u in T)
            ok = sum(1 for u in T if any(w >> u & 1 and not w & (Tm ^ (1 << u)) for w in words))
            worst = max(worst, ok)
    print(f'  n = 8, all 576 sets T of size >= 12 in F_2^4 - 0: at most {worst} points of T have a word of RM(2,4)'
          f' through u and off the rest of T (all of T would be needed)')
    if not luts:
        return
    n, rk_type, pairs_seen = 8, {}, 0
    P = pairs(n)
    for line in open(luts, encoding='utf-8'):
        F = [int(t) for t in re.findall(r'\d+', line)]
        if len(F) != 256:
            continue
        forms = forms_of_lut(F, n)
        rad, cover = {}, [0] * 256
        for v in range(1, 256):
            c = 0
            for k in range(n):
                if v >> k & 1:
                    c ^= forms[k]
            rows = rows_of(c, n, P)
            R = [x for x in range(256) if all(bin(r & x).count('1') % 2 == 0 for r in rows)]
            rad[v] = (rows, R)
            for x in R[1:]:
                cover[x] += 1
        assert all(cover[x] == 1 for x in range(1, 256)), '(P) fails'
        Q = [v for v in rad if len(rad[v][1]) == 16]
        rk_type[len(Q)] = rk_type.get(len(Q), 0) + 1

        def pf_on(rows, basis):      # Pfaffian of the form with these rows, restricted to span(basis)
            def form(x, y):
                s = 0
                for i in range(n):
                    if x >> i & 1:
                        s ^= bin(rows[i] & y).count('1') & 1
                return s
            e = [[form(basis[i], basis[j]) for j in range(4)] for i in range(4)]
            return e[0][1] & e[2][3] ^ e[0][2] & e[1][3] ^ e[0][3] & e[1][2]
        for c in Q:
            basis = []
            for x in rad[c][1]:
                if x and rank2(basis + [x]) > len(basis):
                    basis.append(x)
            for c2 in Q:
                if c2 == c:
                    continue
                assert len(rad[c ^ c2][1]) == 1, '(S) fails'
                assert pf_on(rad[c2][0], basis) == 1, '(N) fails'
                assert pf_on(rad[c ^ c2][0], basis) == 1, 'B_c vanishing on R_c fails'
                pairs_seen += 1
    print(f'  {luts}: (P) on every function; (S) and (N) on all {pairs_seen} ordered pairs of rank-4 components'
          f' (functions by number of rank-4 components: {dict(sorted(rk_type.items()))})')


def table1_entries(r):
    """Kasami-Tokura-Azumi Table I entries of degree r: (number, variables used, f(X), stated weight as f(m))."""
    def pr(X, lo, hi):
        out = np.ones_like(X[1])
        for i in range(lo, hi + 1):
            out &= X[i]
        return out
    A = lambda X: pr(X, 1, r)                     # X1...Xr
    B = lambda X: pr(X, r + 1, 2 * r - 2)         # X(r+1)...X(2r-2)
    E = [(11, 2 * r - 1, lambda X: A(X) ^ ((X[1] ^ 1) & pr(X, r + 1, 2 * r - 1)), lambda m: 2 ** (m - r + 1))]
    if r >= 6:
        E.append((19, 2 * r - 2, lambda X: A(X) ^ (B(X) & ((X[r - 3] & X[r - 2]) ^ (X[r - 1] & (X[r] ^ 1)))),
                  lambda m: 2 ** (m - r + 1) + 2 ** (m - r - 1) - 2 ** (m - 2 * r + 3)))
    if r >= 5:
        E += [(20, 2 * r - 1, lambda X: A(X) ^ (B(X) & ((X[2 * r - 1] & (X[r] ^ 1)) ^ (X[r - 2] & X[r - 1]))),
               lambda m: 2 ** (m - r + 1) + 2 ** (m - r - 1) - 2 ** (m - 2 * r + 3)),
              (21, 2 * r - 1, lambda X: A(X) ^ (B(X) & ((X[2 * r - 1] & X[r - 2]) ^ (X[r - 1] & X[r]))),
               lambda m: 2 ** (m - r + 1) + 2 ** (m - r - 1) - 2 ** (m - 2 * r + 2)),
              (22, 2 * r, lambda X: A(X) ^ (B(X) & ((X[2 * r - 1] & X[r - 1]) ^ (X[2 * r] & X[r]))),
               lambda m: 2 ** (m - r + 1) + 2 ** (m - r - 1) - 2 ** (m - 2 * r + 2)),
              (24, 2 * r, lambda X: A(X) ^ (B(X) & ((X[2 * r - 1] & X[2 * r]) ^ (X[r - 1] & (X[r] ^ 1)))),
               lambda m: 2 ** (m - r + 1) + 2 ** (m - r - 1) - 2 ** (m - 2 * r + 1)),
              (26, 2 * r, lambda X: A(X) ^ (B(X) & ((X[2 * r - 1] & X[2 * r]) ^ (X[r - 1] & X[r]))),
               lambda m: 2 ** (m - r + 1) + 2 ** (m - r - 1) - 2 ** (m - 2 * r + 3) + 2 ** (m - 2 * r + 1))]
    E += [(23, 2 * r + 1, lambda X: A(X) ^ (B(X) & ((X[2 * r - 1] & X[2 * r]) ^ (X[2 * r + 1] & X[r]))),
           lambda m: 2 ** (m - r + 1) + 2 ** (m - r - 1) - 2 ** (m - 2 * r + 2)),
          (25, 2 * r + 1, lambda X: A(X) ^ (B(X) & ((X[2 * r - 1] & X[2 * r]) ^ (X[2 * r + 1] & (X[r] ^ 1)))),
           lambda m: 2 ** (m - r + 1) + 2 ** (m - r - 1) - 2 ** (m - 2 * r + 1)),
          (27, 2 * r + 2, lambda X: A(X) ^ (B(X) & ((X[2 * r - 1] & X[2 * r]) ^ (X[2 * r + 1] & X[2 * r + 2]))),
           lambda m: 2 ** (m - r + 1) + 2 ** (m - r - 1) - 2 ** (m - 2 * r + 2) + 2 ** (m - 2 * r))]
    if r == 5:
        E += [(15, 9, lambda X: pr(X, 1, 5) ^ (X[6] & X[7] & ((X[3] & X[4] & X[8]) ^ (X[2] & (X[5] ^ X[8]) & X[9]))),
               lambda m: 2 ** (m - 4) + 2 ** (m - 7) + 2 ** (m - 8)),
              (18, 10, lambda X: pr(X, 1, 5) ^ (X[6] & X[7] & ((X[3] & X[4] & X[8]) ^ ((X[5] ^ X[8]) & X[9] & X[10]))),
               lambda m: 2 ** (m - 4) + 2 ** (m - 7) + 2 ** (m - 8) + 2 ** (m - 9))]
    return E


def affine_stabilizer_order(S, m):
    """Number of affine maps of F_2^m mapping the point set S onto itself (S must span F_2^m affinely)."""
    pts = sorted(S)
    p0, vs = pts[0], []
    for p in pts[1:]:
        if rank2(vs + [p ^ p0]) > len(vs):
            vs.append(p ^ p0)
    assert len(vs) == m
    coef = {}
    for c in range(1 << m):
        y = p0
        for i in range(m):
            if c >> i & 1:
                y ^= vs[i]
        coef[y] = c
    level = [[x for x in range(1 << m) if coef[x].bit_length() == k] for k in range(m + 1)]
    count = 0

    def image(x, q0, ws):
        y, c, i = q0, coef[x], 0
        while c:
            if c & 1:
                y ^= ws[i]
            c >>= 1
            i += 1
        return y

    def rec(k, q0, ws, red):
        nonlocal count
        if k > m:
            count += 1
            return
        for q in pts:
            w = v = q ^ q0
            for b in red:                 # reduce against the images so far: independent iff something is left
                v = min(v, v ^ b)
            if not v:
                continue
            ws2 = ws + [w]
            for x in level[k]:
                if (x in S) != (image(x, q0, ws2) in S):
                    break
            else:
                rec(k + 1, q0, ws2, red + [v])

    for q0 in pts:
        rec(1, q0, [], [])
    return count


def check8():
    print('== 8. Theorem 6: weights of RM(r, 2r) in [2d, 2.5d) that are 2 mod 4')
    from fractions import Fraction as Fr
    from math import prod
    # (a) Kasami-Tokura-Azumi's list of all weights in [2d, 2.5d) (pp. 392-394), in units of d; m = 2r, ceil(m/r) = 2
    for r in range(5, 13):
        m, d = 2 * r, 2 ** r
        pats = [Fr(2), Fr(17, 8), Fr(9, 4), Fr(37, 16), Fr(19, 8), Fr(39, 16)]                        # (a)-(f)
        pats += [Fr(5, 2) - Fr(1, 2 ** (k + 1)) for k in range(4, m - r - 2 + 1)]                     # (g)
        pats += [Fr(5, 2) - Fr(1, 2 ** (k + 1)) + Fr(1, 2 ** (k + 3)) for k in range(2, m - r - 4 + 1)]  # (h)
        two = sorted({int(p * d) for p in pats if (p * d) % 4 == 2})
        print(f'  r = {r}, d = {d}: listed weights that are 2 mod 4: {two} = 2.5d - {[5 * d // 2 - w for w in two]}')
        assert [5 * d // 2 - w for w in two] == [6, 2]
    # (b) Table I's degree-r entries in at most 2r variables, evaluated directly
    for r in (5, 6, 7):
        m = 2 * r
        idx = np.arange(1 << m, dtype=np.uint32)
        X = [None] + [((idx >> (i - 1)) & 1).astype(np.uint8) for i in range(1, m + 1)]
        got = []
        for no, nv, f, stated in table1_entries(r):
            if nv > m:
                continue
            w = int(f(X).sum())
            assert w == stated(m), (r, no, w, stated(m))
            got.append((no, w))
        two = [(no, w) for no, w in got if w % 4 == 2]
        print(f'  r = {r}: Table I entries (number, weight) {got}, all as stated; 2 mod 4: {two}')
    # (c) dimension 8: KTA's count of weight-34 words of RM(4,8) (formula (b), last term) is one affine orbit
    S = set()
    for p in range(256):
        x = [None] + [p >> i & 1 for i in range(8)]
        if (x[1] & x[2] & x[3] & x[4]) ^ (x[5] & x[6] & ((x[7] & x[8]) ^ (x[3] & x[4]))):
            S.add(p)
    st = affine_stabilizer_order(S, 8)
    agl = 2 ** 8 * prod(2 ** 8 - 2 ** i for i in range(8))
    N = 2 ** 23 * 3 * 5 ** 2 * 7 ** 2 * 17 * 31 * 127
    print(f'  RM(4,8): x1x2x3x4 + x5x6(x7x8 + x3x4) has weight {len(S)}, affine stabilizer of order {st},'
          f' orbit {agl // st}; Kasami-Tokura-Azumi count of weight-34 words {N}: one orbit {agl // st == N}')


if __name__ == '__main__':
    args = sys.argv[1:]
    luts = None
    if '--luts' in args:
        i = args.index('--luts')
        luts = args[i + 1]
        del args[i:i + 2]
    which = set(args) or {'1', '2', '3', '4', '5', '6', '7', '8'}
    apn = check1(luts) if '1' in which or '5' in which else []
    if '2' in which:
        check2()
    if '3' in which:
        check3()
    if '4' in which:
        check4()
    if '5' in which:
        check5(apn)
    if '6' in which:
        check6()
    if '7' in which:
        check7(luts)
    if '8' in which:
        check8()
