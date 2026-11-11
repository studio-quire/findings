"""Aut(P) in GL(8,2) for each of five's configurations, computed from the solids alone, without pspread's stabiliser.

  python apn/five_aut.py runs/i12/p5.txt runs/i12/aut.txt

A g fixing the set of solids sends the frame (R0, R1, R2) to an ordered triple of solids (A, B, C). With V = R0 + R1 and R2 the
graph of phi: R0 -> R1, and C the graph of psi: A -> B, such g are exactly g|R0 = alpha (any isomorphism R0 -> A), g|R1 = psi alpha
phi^-1: 60 triples x 20,160 alphas, each tested on the two other solids. The count must equal the a<N> in the configuration's tag.

Output per configuration: `<tag> <|Aut|> <generators>`, a generator as the images of e0..e7, 8 hex bytes joined, chosen greedily
until they generate the whole group (each closure checked to have exactly |Aut| elements).
Rowen (Claude Opus 5.5), 2026-11-06.
"""
import itertools, re, sys
import numpy as np


def span(vs):
    s = {0}
    for v in vs:
        s |= {x ^ v for x in s}
    return s


def ordered_bases(pts):
    """Every ordered basis of the 4-space pts, as an (n, 4) array."""
    out, nz = [], sorted(pts - {0})
    for a in nz:
        for b in nz:
            if b == a:
                continue
            s2 = span([a, b])
            for c in nz:
                if c in s2:
                    continue
                s3 = span([a, b, c])
                for d in nz:
                    if d not in s3:
                        out.append((a, b, c, d))
    return np.array(out, dtype=np.uint8)


def coords(x, basis):
    """x as a bitmask over basis, or None."""
    for m in range(1 << len(basis)):
        v = 0
        for i, b in enumerate(basis):
            if m >> i & 1:
                v ^= b
        if v == x:
            return m
    return None


def split(x, b0, b1):
    """x = u + v, u in <b0>, v in <b1>: masks over b0 and b1."""
    m = coords(x, b0 + b1)
    return m & 15, m >> 4


def apply_mat(g, x):
    y = 0
    for i in range(8):
        if x >> i & 1:
            y ^= g[i]
    return y


def compose(g, h):   # g after h
    return tuple(apply_mat(g, h[i]) for i in range(8))


def closure(gens):
    ident = tuple(1 << i for i in range(8))
    seen, todo = {ident}, [ident]
    while todo:
        new = []
        for x in todo:
            for g in gens:
                y = compose(g, x)
                if y not in seen:
                    seen.add(y)
                    new.append(y)
        todo = new
    return seen


def group(bases):
    sets = [span(b) for b in bases]
    k = len(bases)
    b0, b1, b2 = bases[0], bases[1], bases[2]
    # phi: R0 -> R1 with R2 its graph: for r in R2, r = u + phi(u)
    phi_inv = {}   # mask over b1 -> mask over b0
    for r in sets[2]:
        u, v = split(r, b0, b1)
        phi_inv[v] = u
    assert len(phi_inv) == 16
    # each vector needed, as (mask over b0, mask over b0 of phi^-1 of its R1 part)
    def parts(x):
        u, v = split(x, b0, b1)
        return u, phi_inv[v]
    eparts = [parts(1 << i) for i in range(8)]
    rest = [[parts(x) for x in bases[i]] for i in range(3, k)]

    def lin(al, mask):
        y = np.zeros(len(al), dtype=np.uint8)
        for i in range(4):
            if mask >> i & 1:
                y ^= al[:, i]
        return y

    alphas = {i: ordered_bases(sets[i]) for i in range(k)}
    member = np.zeros((k, 256), dtype=bool)
    for i, s in enumerate(sets):
        member[i, list(s - {0})] = True
    found = []
    for a, b, c in itertools.permutations(range(k), 3):
        if sets[a] & sets[b] != {0}:
            continue
        psi = np.zeros(256, dtype=np.uint8)
        for r in sets[c]:   # r = x + psi(x), x in A, psi(x) in B
            x = next(p for p in sets[a] if (r ^ p) in sets[b])
            psi[x] = r ^ x
        al = alphas[a]

        def image(p):
            return lin(al, p[0]) ^ psi[lin(al, p[1])]

        others = [i for i in range(k) if i not in (a, b, c)]
        ok = np.zeros(len(al), dtype=bool)
        for perm in itertools.permutations(others):   # R3.. -> these solids, in this order
            good = np.ones(len(al), dtype=bool)
            for i, tgt in enumerate(perm):
                for p in rest[i]:
                    good &= member[tgt, image(p)]
            ok |= good
        if ok.any():
            imgs = np.stack([image(p) for p in eparts], axis=1)[ok]
            found += [tuple(int(v) for v in row) for row in imgs]
    return found


def main():
    src, out = sys.argv[1:3]
    with open(out, "w") as o:
        for l in open(src):
            t = l.split()
            if not t:
                continue
            tag, k = t[0], int(t[1])
            by = [int(x, 16) for x in t[2:2 + 4 * k]]
            bases = [by[4 * i:4 * i + 4] for i in range(k)]
            els = group(bases)
            want = int(re.search(r"-a(\d+)-", tag)[1])
            sets = {frozenset(span(b)) for b in bases}
            for g in els:
                assert len({apply_mat(g, x) for x in range(256)}) == 256
                assert {frozenset(apply_mat(g, x) for x in s) for s in sets} == sets
            assert len(set(els)) == len(els) == want, (tag, len(els), want)
            gens, have = [], {tuple(1 << i for i in range(8))}
            for g in els:
                if g not in have:
                    gens.append(g)
                    have = closure(gens)
                    if len(have) == want:
                        break
            assert len(have) == want, (tag, len(have), want)
            o.write(f"{tag} {want} " + " ".join("".join(f"{x:02x}" for x in g) for g in gens) + "\n")
            print(tag, want, "generators", len(gens), flush=True)


if __name__ == "__main__":
    main()
