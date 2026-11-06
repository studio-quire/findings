"""Positive controls for cover.cpp's hole stage (the leaves of holes_tree must reach push and on to SOLUTION): known
functions cut down to a D-dim subspace holding all their rank-4 forms, so the hole stage has to complete the rest.

  python apn/cover_ctl.py gen OUT.txt LUTFILE N D [SEED]   N functions (most rank-4 forms first, then at random) -> cover.cpp
      input lines "<tag> <k> <solid bases> seed <coefs>", and OUT.txt.spaces with each function's whole space
  python apn/cover_ctl.py check OUT.txt.spaces COVER.out   every planted space among the SOLUTION lines, nothing else
      that isn't checked to be APN with the same solids (extra spaces are printed; each is verified by its own forms)

Rowen (Claude Opus 5.5), 2026-11-04. (cover_prep.py ctl seeds only rank-4 forms and their span, too broad for this.)
"""
import random, re, sys

PAIRS = [(i, j) for i in range(8) for j in range(i + 1, 8)]
par = lambda v: bin(v).count("1") & 1


def rows_of(F, b):
    f = [par(b & F[x]) for x in range(256)]
    M = [0] * 8
    for i, j in PAIRS:
        if f[(1 << i) ^ (1 << j)] ^ f[1 << i] ^ f[1 << j] ^ f[0]:
            M[i] |= 1 << j; M[j] |= 1 << i
    return M


def coef(M):
    return sum(1 << k for k, (i, j) in enumerate(PAIRS) if M[i] >> j & 1)


def rows_from(c):
    M = [0] * 8
    for k, (i, j) in enumerate(PAIRS):
        if c >> k & 1:
            M[i] |= 1 << j; M[j] |= 1 << i
    return M


def radical(M):
    out = []
    for x in range(1, 256):
        v = 0
        for i in range(8):
            if x >> i & 1:
                v ^= M[i]
        if v == 0:
            out.append(x)
    return out


def rref(vs):
    piv = {}
    for v in vs:
        for b in range(27, -1, -1):
            if v >> b & 1:
                if b in piv:
                    v ^= piv[b]
                else:
                    piv[b] = v
                    break
    for b in sorted(piv):
        for c in piv:
            if c != b and piv[c] >> b & 1:
                piv[c] ^= piv[b]
    return tuple(sorted(piv.values()))


def span_points(basis):
    pts = [0]
    for b in basis:
        if b not in pts:
            pts += [p ^ b for p in pts]
    return pts


def indep(vs):
    """a maximal independent subfamily, in order"""
    out, r = [], ()
    for v in vs:
        r2 = rref(list(r) + [v])
        if len(r2) > len(r):
            out.append(v); r = r2
    return out


def analyse(F):
    forms = {b: rows_of(F, b) for b in range(1, 256)}
    r4, r6, r8, rads = [], [], [], {}
    for b, M in forms.items():
        R = radical(M)
        c = coef(M); rads[c] = R
        (r4 if len(R) == 15 else r6 if len(R) == 3 else r8).append(c)
    return r4, r6, r8, rads


def solid_basis(R):
    out, s = [], {0}
    for x in R:
        if x not in s:
            out.append(x); s |= {y ^ x for y in s}
    return out


def cmd_gen(out, lutfile, n, D, seed=6):
    rnd = random.Random(seed)
    Fs = []
    for line in open(lutfile):
        line = line.strip()
        if line.startswith("{"):
            Fs.append([int(t) for t in line[1:line.index("}")].split(",")])
    rnd.shuffle(Fs)
    scored = []
    for F in Fs[:400]:
        r4 = analyse(F)[0]
        scored.append((-len(r4), len(scored), F))
    scored.sort(key=lambda t: t[:2])
    with open(out, "w") as f, open(out + ".spaces", "w") as g:
        made = 0
        for _, idx, F in scored:
            if made == n:
                break
            r4, r6, r8, rads = analyse(F)
            base = indep(r4)
            if len(r4) < 2 or len(base) > D:
                continue
            rnd.shuffle(r6); rnd.shuffle(r8)
            seeds = indep(base + r6 + r8)[:D]
            if not all(c in span_points(seeds) for c in r4):
                continue
            tag = f"ctl{D}-{idx}-k{len(r4)}"
            sol = " ".join(f"{x:02x}" for c in r4 for x in solid_basis(rads[c]))
            f.write(f"{tag} {len(r4)} {sol} seed " + " ".join(f"{c:07x}" for c in seeds) + "\n")
            whole = rref(r4 + r6 + r8)
            g.write(f"{tag} " + " ".join(f"{c:07x}" for c in whole) + "\n")
            made += 1
    print(f"{made} controls at D = {D} (rank-4 counts: first has {-scored[0][0]})")


def is_apn(basis):
    seen, pts = set(), span_points(list(basis))
    for c in pts:
        if c:
            for x in radical(rows_from(c)):
                if x in seen:
                    return False
                seen.add(x)
    return len(pts) == 256 and len(seen) == 255


def cmd_check(spaces, coverout):
    want = {}
    for l in open(spaces):
        t = l.split()
        want[t[0]] = tuple(int(x, 16) for x in t[1:])
    found = {k: [] for k in want}
    bad = 0
    for l in open(coverout):
        if l.startswith("SOLUTION "):
            t = l.split()
            found[t[1]].append(rref(int(x, 16) for x in t[t.index("basis") + 1:]))
        elif l.startswith(("BADSEED", "FULL")):
            print(l.strip()); bad += 1
    miss = extra = notapn = 0
    for k, w in want.items():
        if w not in found[k]:
            miss += 1; print("MISSED", k, f"({len(found[k])} solutions)")
        for s in found[k]:
            if s != w:
                extra += 1
                if not is_apn(s):
                    notapn += 1; print("NOT APN", k, s)
    print(f"{len(want)} controls: {len(want) - miss} planted spaces found, {miss} missed; {extra} other solutions, "
          f"{notapn} of them not APN; {bad} bad seeds or full spans")
    return 1 if miss or notapn or bad else 0


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "gen":
        cmd_gen(a[1], a[2], int(a[3]), int(a[4]), int(a[5]) if len(a) > 5 else 6)
    else:
        sys.exit(cmd_check(a[1], a[2]))
