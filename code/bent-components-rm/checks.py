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


if __name__ == '__main__':
    args = sys.argv[1:]
    luts = None
    if '--luts' in args:
        i = args.index('--luts')
        luts = args[i + 1]
        del args[i:i + 2]
    which = set(args) or {'1', '2', '3', '4', '5', '6'}
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
