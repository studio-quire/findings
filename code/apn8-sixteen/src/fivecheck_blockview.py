"""Rowen's block-view claims, tested on S8s in the sieve's form format (hex forms, last 8 tokens of each line).

For every ordered pair (R_i, R_j) of rank-4 radicals of S8: V = R_i + R_j, B = (P, M, Q); K = S8 cap {M = 0}, C = M(S8).
Claims: dim K in {2, 3}; K/<B_i> -> P injective with nondegenerate nonzero images; C^perp has no rank-1 matrix.
The search doesn't rest on these; this is the check the brief asks for, on what actually exists.
"""
import sys, itertools

def parity(x): return bin(x).count('1') & 1

def ev(f, a, b):
    r = 0
    for i in range(8):
        if a >> i & 1:
            r ^= (f >> (8 * i)) & 255
    return parity(r & b)

def span(vs):
    out = {0}
    for v in vs:
        out |= {x ^ v for x in out}
    return out

def rank(vs):
    piv = {}
    n = 0
    for v in vs:
        while v:
            h = v.bit_length() - 1
            if h in piv:
                v ^= piv[h]
            else:
                piv[h] = v; n += 1
                break
    return n

def basis_of(mask_elems):
    b, seen = [], {0}
    for v in sorted(mask_elems):
        if v not in seen:
            b.append(v); seen |= {x ^ v for x in seen}
    return b

def radical(f):
    return [v for v in range(1, 256) if all(ev(f, v, 1 << i) == 0 for i in range(8))]

def main():
    lines = [l.split() for l in open(sys.argv[1]) if l.strip()]
    tally = {}
    for tok in lines:
        fs = [int(x, 16) for x in tok[-8:]]
        elems = []
        for m in range(1, 256):
            f = 0
            for i in range(8):
                if m >> i & 1:
                    f ^= fs[i]
            elems.append(f)
        rads = [(f, radical(f)) for f in elems]
        solids = [(f, r) for f, r in rads if len(r) == 15]
        for (fi, ri), (fj, rj) in itertools.permutations(solids, 2):
            bi, bj = basis_of(ri), basis_of(rj)
            # M block of each basis form of S8
            def mpart(f):
                x = 0
                for p, a in enumerate(bi):
                    for q, b in enumerate(bj):
                        x |= ev(f, a, b) << (4 * p + q)
                return x
            def ppart(f):
                x, k = 0, 0
                for p in range(4):
                    for q in range(p + 1, 4):
                        x |= ev(f, bi[p], bi[q]) << k; k += 1
                return x
            ms = [mpart(f) for f in fs]
            dimC = rank(ms)
            dimK = 8 - dimC
            # K's elements: all forms with M = 0
            K = [f for f in elems if mpart(f) == 0]
            assert len(K) == 2 ** dimK - 1
            assert fi in K and fj in K
            # P injective on K / <B_i>, images nondegenerate
            for f in K:
                if f == fi:
                    continue
                p = ppart(f)
                assert p != 0
                pf = (p & 1) * (p >> 5 & 1) ^ (p >> 1 & 1) * (p >> 4 & 1) ^ (p >> 2 & 1) * (p >> 3 & 1)
                assert pf == 1, 'P degenerate'
            # C^perp has no rank-1: for each a, y nonzero some M in C has a^T M y = 1
            C = span(ms)
            for a in range(1, 16):
                for y in range(1, 16):
                    rk1 = 0
                    for p in range(4):
                        for q in range(4):
                            if (a >> p & 1) and (y >> q & 1):
                                rk1 |= 1 << (4 * p + q)
                    assert any(parity(rk1 & c) for c in ms), 'rank-1 in C perp'
            tally[dimK] = tally.get(dimK, 0) + 1
    print('ordered pairs by dim K:', tally, 'all block-view claims held')

main()
