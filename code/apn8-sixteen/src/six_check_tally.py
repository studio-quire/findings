"""Tallies and comparisons for six_check (search B, 6 rank-4 components).

Subagent (Claude Opus 5.5) for Rowen, 2026-11-03. Own code (Python integers, no shared routines with six_check.cpp).
  python apn/six_check_tally.py prog <configs.txt> <prog>...             totals of R (point branching) or L (lin) lines
  python apn/six_check_tally.py cores <core file A> <core file B> [...]  core sets compared per configuration
  python apn/six_check_tally.py nodes <R prog> <L prog>                  lin nodes less its x4 nodes = point branching's
Core files: K lines (six_check search: K conf dim <rref codes, 28-bit, pair order (0,1),(0,2)..(6,7)>) and C lines
(six_check lin: C conf dim <basis forms as u64, byte i = row i>). Each core becomes the set of its nonzero elements
as 28-bit integers (pair (i, j), i < j, at bit index in the order above); lin's cores with a rank-4 element that isn't
one of the 28 forms with a member radical (a seventh solid; rank computed here) are set aside (x4).

This folder's copy (the note added by Claude Opus 5.5, the code unchanged): the files here for `prog` and `nodes` are
logs/check6-{p1asc,p1desc,asc,asc340}-R.txt (R lines) and logs/check6-{lin,lindesc}-prog.txt (L lines; the lin runs'
.prog files whole), from the bundle root, e.g.
  python src/six_check_tally.py prog data/configs6.txt logs/check6-lin-prog.txt
  python src/six_check_tally.py nodes logs/check6-asc-R.txt logs/check6-lin-prog.txt
The core files `cores` compares (p1asc.cores, p1desc.cores, lin.cores: 457 MB to 973 MB each) aren't in this folder;
its results are in logs/check6-tally-p1-output.txt and logs/check6-tally-lin-output.txt.
"""
import sys
from collections import defaultdict

PAIRS = [(i, j) for i in range(8) for j in range(i + 1, 8)]


def from_u64(f):
    v = 0
    for p, (i, j) in enumerate(PAIRS):
        if (f >> (8 * i + j)) & 1:
            v |= 1 << p
    return v


def rows(v):
    r = [0] * 8
    for p, (i, j) in enumerate(PAIRS):
        if (v >> p) & 1:
            r[i] |= 1 << j
            r[j] |= 1 << i
    return r


def rank(v):
    r = rows(v)
    rk = 0
    for c in range(8):
        piv = next((k for k in range(rk, 8) if (r[k] >> c) & 1), None)
        if piv is None:
            continue
        r[rk], r[piv] = r[piv], r[rk]
        for k in range(8):
            if k != rk and (r[k] >> c) & 1:
                r[k] ^= r[rk]
        rk += 1
    return rk


def rref(basis):
    # canonical basis: fully reduced, pivots = highest bits, sorted
    out = []
    for v in basis:
        for r in out:
            if (v >> (r.bit_length() - 1)) & 1:
                v ^= r
        if not v:
            return None
        p = v.bit_length() - 1
        out = [r ^ v if (r >> p) & 1 else r for r in out]
        out.append(v)
    return tuple(sorted(out))


def span(basis):
    s = {0}
    for b in basis:
        s |= {x ^ b for x in s}
    s.discard(0)
    return s


def stream_cores(path, lo=0, hi=10 ** 9):
    # {conf: {canonical basis}} for lo <= conf < hi (a configuration written twice, by a restarted run, is merged;
    # repeated cores counted)
    cur, cores, dup = None, set(), 0
    blocks = {}
    for line in open(path):
        t = line.split(maxsplit=3)
        if not t or t[0] not in ("K", "C"):
            continue
        conf = int(t[1])
        if not lo <= conf < hi:
            continue
        t = line.split()
        dim = int(t[2])
        basis = [int(x, 16) for x in t[3:]] if t[0] == "K" else [from_u64(int(x, 16)) for x in t[3:]]
        key = rref(basis)
        if key is None or len(key) != dim:
            raise SystemExit(f"{path}: core of conf {conf} not of dimension {dim}")
        if conf != cur:
            if cur is not None:
                blocks.setdefault(cur, set()).update(cores)
            cur, cores = conf, set()
        if key in cores or key in blocks.get(conf, ()):
            dup += 1
        cores.add(key)
    if cur is not None:
        blocks.setdefault(cur, set()).update(cores)
    return blocks, dup


def n_rank4(key):
    return sum(1 for v in span(key) if rank(v) == 4)


def cmd_cores(files):
    # the first file is the reference (point branching, K lines); lin's (C lines) cores beyond it must all hold a
    # seventh solid (more than 6 rank-4 elements), and none of the reference's may be missing
    # in chunks of 40 configurations (memory)
    tot = {f: [0, 0, 0, 0, 0, 0, 0, 0] for f in files}  # cores, configurations, repeats, missing confs, x4, other extra, ref cores (common configurations), configurations in one file only
    for lo in range(0, 400, 40):
        ref, dup0 = stream_cores(files[0], lo, lo + 40)
        tot[files[0]][0] += sum(len(v) for v in ref.values()); tot[files[0]][1] += len(ref); tot[files[0]][2] += dup0
        for f in files[1:]:
            oth, dup = stream_cores(f, lo, lo + 40)
            T = tot[f]
            T[0] += sum(len(v) for v in oth.values()); T[1] += len(oth); T[2] += dup
            T[7] += len(set(ref) ^ set(oth))
            for k in sorted(set(ref) & set(oth)):
                a, b = ref[k], oth[k]
                T[6] += len(a)
                if a - b:
                    T[3] += 1
                    if T[3] <= 5:
                        print(f"  conf {k}: {len(a - b)} of {files[0]}'s missing in {f}")
                for key in b - a:
                    if n_rank4(key) > 6:
                        T[4] += 1
                    else:
                        T[5] += 1
                        if T[5] <= 5:
                            print(f"  conf {k}: an extra core in {f} without a seventh solid")
    print(f"{files[0]}: {tot[files[0]][0]} cores on {tot[files[0]][1]} configurations, {tot[files[0]][2]} repeated")
    for f in files[1:]:
        T = tot[f]
        print(f"{f}: {T[0]} cores on {T[1]} configurations, {T[2]} repeated; {T[7]} configurations in only one of the two; on the common ones against {files[0]} ({T[6]} cores): "
              f"{T[3]} configurations with cores missing, {T[4]} extra cores with a seventh solid, {T[5]} other extra cores")


def read_prog(path, tag):
    out = {}
    for line in open(path):
        t = line.split()
        if t and t[0] == tag:
            k = int(t[1])
            if k in out:
                raise SystemExit(f"{path}: configuration {k} twice")
            out[k] = t
    return out


def cmd_nodes(rp, lp):
    R = read_prog(rp, "R")
    L = read_prog(lp, "L")
    bad = 0
    for k in sorted(set(R) & set(L)):
        r, l = R[k], L[k]
        rn = [int(x) for x in r[2:11]]
        ln = [int(x) for x in l[3:12]]
        lx = [int(x) for x in l[21:30]]
        if [a - b for a, b in zip(ln, lx)] != rn:
            bad += 1
            if bad <= 5:
                print(f"  conf {k}: R {rn} L {ln} x4 {lx}")
    print(f"nodes by dimension, point branching (R) vs lin less x4 (L): {len(set(R) & set(L))} configurations in both, {bad} differ;"
          f" only R {len(set(R) - set(L))}, only L {len(set(L) - set(R))}")


def cmd_prog(conf, progs):
    n = sum(1 for line in open(conf) if line.startswith("CONF"))
    tot = defaultdict(lambda: [0] * 60)
    seen = defaultdict(set)
    for p in progs:
        for line in open(p):
            t = line.split()
            if not t or t[0] not in ("R", "L"):
                continue
            k = int(t[1])
            if k in seen[t[0]]:
                print(f"  {p}: configuration {k} twice")
            seen[t[0]].add(k)
            vals = [int(x) for x in (t[2:-12] if t[0] == "R" else t[3:-1])]
            for i, v in enumerate(vals):
                tot[t[0]][i] += v
            if t[0] == "R" and len(t) > 58:
                tot["Rsw"][0] += int(t[-1])
    for tag in ("R", "L"):
        if not seen[tag]:
            continue
        v = tot[tag]
        print(f"{tag}: {len(seen[tag])} of {n} configurations")
        if tag == "R":
            print(f"  nodes by dim {v[0:9]}\n  tries by dim {v[9:18]}\n  p1 by dim {v[18:27]}\n  bentext {v[27]} full {v[28]}"
                  f"\n  tries in D by dim {v[29:38]}\n  cores by dim {v[38:47]}")
        else:
            print(f"  nodes by dim {v[0:9]}\n  steps by dim {v[9:18]}\n  x4 nodes by dim {v[18:27]}\n  cores by dim {v[27:36]}"
                  f"\n  x4 cores by dim {v[36:45]}\n  solutions {v[45]} dead {v[46]}")


if __name__ == "__main__":
    c, a = sys.argv[1], sys.argv[2:]
    if c == "cores":
        cmd_cores(a)
    elif c == "nodes":
        cmd_nodes(*a)
    elif c == "prog":
        cmd_prog(a[0], a[1:])
