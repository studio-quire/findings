"""For each class of r solids: all solids inside the union, and all exact covers of the union by r solids.

Tells a class of solid sets from a class of point sets: Langevin-Hou's fix-5 and Kurz's Table 9.3 count point sets, and a point
set with c covers (all in one class) has a stabiliser c times the stabiliser of one cover.
"""
import sys
from fractions import Fraction

def matvec(T, v):
    r = 0
    for i in range(4):
        r |= (bin((T >> (4 * i)) & 15 & v).count('1') & 1) << i
    return r

def solids_of(r, T4, T5):
    Ts = [0x8421, T4, T5][: r - 2]
    S = [frozenset(x for x in range(1, 16)), frozenset(y << 4 for y in range(1, 16))]
    for T in Ts:
        S.append(frozenset(x | (matvec(T, x) << 4) for x in range(1, 16)))
    return S

def solids_inside(U):
    """All 4-dim subspaces whose nonzero points lie in the set U."""
    pts = sorted(U)
    found = set()
    for a in pts:
        for b in pts:
            if b <= a or (a ^ b) not in U:
                continue
            s2 = {a, b, a ^ b}
            for c in pts:
                if c <= b or c in s2 or any((c ^ x) not in U for x in s2):
                    continue
                s3 = s2 | {c} | {c ^ x for x in s2}
                for d in pts:
                    if d <= c or d in s3 or any((d ^ x) not in U for x in s3):
                        continue
                    found.add(frozenset(s3 | {d} | {d ^ x for x in s3}))
    return found

def exact_covers(U, sol, r):
    sol = list(sol)
    out = []
    def go(rest, chosen):
        if not rest:
            out.append(frozenset(chosen)); return
        p = min(rest)
        for s in sol:
            if p in s and s <= rest:
                go(rest - s, chosen + [s])
    go(frozenset(U), [])
    return set(out)


def normal_pair(cover):
    """(T, T') of an ordered cover: coordinates with the first three as R1, R2, R3."""
    W = [sorted(s) for s in cover]
    e = []; sp = {0}
    for v in W[0]:
        if v not in sp:
            e.append(v); sp |= {x ^ v for x in sp}
    f = []
    for x in e:
        f.append(next(y for y in W[1] if (x ^ y) in cover[2]))
    # coordinates of v in the basis e, f by table
    table = {}
    for m in range(256):
        v = 0
        for i in range(4):
            if m >> i & 1: v ^= e[i]
            if m >> (4 + i) & 1: v ^= f[i]
        table[v] = m
    Ts = []
    for s in cover[3:]:
        T = 0
        for b in range(4):
            w = next(table[v] for v in s if table[v] & 15 == 1 << b)
            for i in range(4):
                if (w >> 4) >> i & 1:
                    T |= 1 << (4 * i + b)
        Ts.append(T)
    return Ts

def mul(a, b):
    r = 0
    for i in range(4):
        row, acc = (a >> (4 * i)) & 15, 0
        for k in range(4):
            if row >> k & 1: acc ^= (b >> (4 * k)) & 15
        r |= acc << (4 * i)
    return r

def main():
    G = 1
    for i in range(8): G *= 256 - 2 ** i
    rows = [l.split() for l in open('runs/check5/configs.txt') if l.strip()]
    total = {4: Fraction(0), 5: Fraction(0)}
    for i, row in enumerate(rows):
        r, T4, T5, om, gam = int(row[0]), int(row[1], 16), int(row[2], 16), int(row[3]), int(row[4])
        S = solids_of(r, T4, T5)
        U = frozenset().union(*S)
        assert len(U) == 15 * r
        sol = solids_inside(U)
        cov = exact_covers(U, sol, r)
        note = ''
        if len(cov) > 1:
            others = [c for c in cov if c != frozenset(S)]
            # is each other cover disjoint, as a set of solids, from S?
            note = ' shared solids between covers: ' + str(sorted(len(c & frozenset(S)) for c in others))
            for c in others:
                Ts = normal_pair([sorted_s for sorted_s in sorted(c, key=min)])
                I = 0x8421
                sub = [(T, mul(T, T) ^ T ^ I) for T in Ts]
                note += ' other cover T: %s, T^2+T+I: %s' % (['%04x' % t for t in Ts], ['%04x' % x for _, x in sub])
        total[r] += Fraction(G, gam * len(cov))
        print('r=%d rep %04x %04x |Gamma|=%d: solids inside %d, exact covers %d, point-set stabiliser %d%s' % (r, T4, T5, gam, len(sol), len(cov), gam * len(cov), note))
    print('point sets r=5 (each set of solids weighted 1/covers):', total[5])

main()
