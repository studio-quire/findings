"""A second implementation of seven_check's search (7 rank-4 radicals, 150 holes), through full depth.

Adapted for seven-check (Claude Opus 5.5, 2026-10-31) from eight_check_numpy.py, which a subagent (Claude Opus 5.5)
wrote for eight-check from the mathematics, sharing no code with the C++:
- forms are 8 row bytes of a symmetric zero-diagonal matrix (bit j of row i = B(e_i, e_j)); B w = XOR of rows j, j in w;
- the forms with a given radical R are found as the null space, inside the 28-dimensional space of alternating forms,
  of the linear conditions B r = 0 (r in a basis of R), enumerated, and kept when the radical is exactly R;
- legality and radicals of many forms at once in numpy (B w for all 256 w by doubling); disjointness of the new
  radicals as "sum of radical sizes = size of their union".
Same search as seven_check `search ... asc` (first uncovered point, off-D points ascending, then D ascending; branch
over every legal form whose radical contains it), so node counts by dimension, near-misses (nodes with every off-D
point covered) and tries at points of D must agree with the C++ exactly.

  python seven_check_numpy.py search <configs.txt> <index>...       per configuration: nodes / p1 / triesD by dim
  python seven_check_numpy.py control <lut file> <first line> <n> <seed>
      planted functions: a known quadratic APN T with k >= 2 rank-4 radicals; members = its k radicals + 7 - k random
      solids (all pairwise disjoint), T's rank-6 radicals off D allowed as an oracle (as the C++ controls); the search
      must find T's span.
Run with runs/venv's python (numpy).
"""
import random
import re
import sys
import time

import numpy as np

POP = np.array([bin(i).count("1") for i in range(256)], dtype=np.int64)


# ---------------------------------------------------------------- GF(2) linear algebra on Python ints
def span_of(vs):
    s = {0}
    for v in vs:
        if v not in s:
            s |= {x ^ v for x in s}
    return s


def basis_of(points):
    b, s = [], {0}
    for p in points:
        if p not in s:
            b.append(p)
            s |= {x ^ p for x in s}
    return b


def nullspace(rows, nvars):
    """Null space of the linear system rows (ints, bit k = coefficient of variable k) over GF(2)."""
    piv = {}  # pivot column -> reduced row
    for r in rows:
        for c, pr in piv.items():
            if r >> c & 1:
                r ^= pr
        if r:
            c = r.bit_length() - 1
            for c2 in list(piv):
                if piv[c2] >> c & 1:
                    piv[c2] ^= r
            piv[c] = r
    free = [k for k in range(nvars) if k not in piv]
    out = []
    for f in free:
        v = 1 << f
        for c, pr in piv.items():
            if pr >> f & 1:
                v |= 1 << c
        out.append(v)
    return out


PAIRS = [(i, j) for i in range(8) for j in range(i + 1, 8)]


def form_of_code(code):
    rows = [0] * 8
    for k, (i, j) in enumerate(PAIRS):
        if code >> k & 1:
            rows[i] |= 1 << j
            rows[j] |= 1 << i
    return rows


def forms_with_radical_containing(R):
    """All alternating forms whose radical contains the subspace R (as an (n, 8) uint8 array)."""
    conds = []
    for r in basis_of(sorted(R)):
        for bit in range(8):  # bit `bit` of B r = sum over j in r of B(e_bit, e_j)
            row = 0
            for k, (i, j) in enumerate(PAIRS):
                if (i == bit and r >> j & 1) or (j == bit and r >> i & 1):
                    row |= 1 << k
            conds.append(row)
    ns = nullspace(conds, 28)
    basis = np.array([form_of_code(c) for c in ns], dtype=np.uint8).reshape(len(ns), 8)
    allf = np.zeros((1, 8), dtype=np.uint8)
    for b in basis:
        allf = np.concatenate([allf, allf ^ b[None, :]])
    return allf


def mask_of(points):
    m = 0
    for p in points:
        m |= 1 << p
    return m


def packed(m):
    return np.array([(m >> (64 * k)) & ((1 << 64) - 1) for k in range(4)], dtype=np.uint64)


def to_int(row):
    return sum(int(row[k]) << (64 * k) for k in range(4))


# ---------------------------------------------------------------- radicals of many forms at once
def radicals(F):
    """F (n, 8) uint8 -> (size of radical (n,), packed radical masks (n, 4) uint64)."""
    n = F.shape[0]
    if n > 65536:
        parts = [radicals(F[k:k + 65536]) for k in range(0, n, 65536)]
        return np.concatenate([c for c, _ in parts]), np.concatenate([p for _, p in parts])
    V = np.zeros((n, 256), dtype=np.uint8)
    for j in range(8):
        V[:, 1 << j: 2 << j] = V[:, : 1 << j] ^ F[:, j: j + 1]
    Z = V == 0
    Z[:, 0] = False
    cnt = Z.sum(axis=1)
    pk = np.packbits(Z, axis=1, bitorder="little").view(np.uint64)
    return cnt, pk


class Problem:
    def __init__(self, members, oracle=()):
        self.members = [frozenset(m) for m in members]
        self.radof = {}
        for j, m in enumerate(self.members):
            assert len(m) == 15
            for v in m:
                assert v not in self.radof
                self.radof[v] = j
        self.D = frozenset(v for v in range(1, 256) if v not in self.radof)
        self.Dm = packed(mask_of(self.D))
        self.memm = np.array([packed(mask_of(m)) for m in self.members])
        self.F4 = []
        for m in self.members:
            F = forms_with_radical_containing(m)
            cnt, _ = radicals(F)
            self.F4.append(F[cnt == 15])
            assert len(self.F4[-1]) == 28
        lines = set()
        Ds = sorted(self.D)
        for i, a in enumerate(Ds):
            for b in Ds[i + 1:]:
                if (a ^ b) in self.D:
                    lines.add(frozenset((a, b, a ^ b)))
        self.lines = sorted(lines, key=sorted)
        self.lines_at = {v: [] for v in self.D}
        for li, L in enumerate(self.lines):
            for v in L:
                self.lines_at[v].append(li)
        self.lmask = [mask_of(L) for L in self.lines]
        self.lforms = {}
        self.order = sorted(self.radof) + sorted(self.D)
        # oracle (controls only): rank-6 forms whose radicals, off D, count as legal
        self.oracle = [(np.array(f, dtype=np.uint8), frozenset(r)) for f, r in oracle]
        self.orm = np.array([packed(mask_of(r)) for _, r in self.oracle]).reshape(len(self.oracle), 4)

    def line_forms(self, li):
        if li not in self.lforms:
            F = forms_with_radical_containing(self.lines[li])
            cnt, _ = radicals(F)
            self.lforms[li] = F[cnt == 3]
            assert len(self.lforms[li]) == 13888
        return self.lforms[li]

    def legal(self, G):
        cnt, pk = radicals(G)
        ok = cnt == 0
        r4 = cnt == 15
        if r4.any():
            eq = (pk[r4][:, None, :] == self.memm[None, :, :]).all(axis=2).any(axis=1)
            ok[np.flatnonzero(r4)[eq]] = True
        r6 = cnt == 3
        if r6.any():
            inD = ((pk[r6] & ~self.Dm[None, :]) == 0).all(axis=1)
            if len(self.oracle):
                inD |= (pk[r6][:, None, :] == self.orm[None, :, :]).all(axis=2).any(axis=1)
            ok[np.flatnonzero(r6)[inD]] = True
        return ok, cnt, pk


def children(P, span, cov, cand):
    """Candidates f (rows of cand) for which span + f is legal with pairwise disjoint radicals, none meeting cov.
    Returns the surviving candidates and their unions of new radicals (packed)."""
    covp = packed(cov)
    surv = np.arange(len(cand))
    acc = np.zeros((len(cand), 4), dtype=np.uint64)
    tot = np.zeros(len(cand), dtype=np.int64)
    k = 0
    # one span element at a time while there are many candidates, then the rest at once
    while k < len(span) and len(surv) > 4096:
        G = cand[surv] ^ span[k][None, :]
        ok, cnt, pk = P.legal(G)
        ok &= ((pk & covp[None, :]) == 0).all(axis=1) & ((pk & acc[surv]) == 0).all(axis=1)
        surv = surv[ok]
        acc[surv] |= pk[ok]
        tot[surv] += cnt[ok]
        k += 1
    if k < len(span) and len(surv):
        rest = span[k:]
        G = (cand[surv][:, None, :] ^ rest[None, :, :]).reshape(-1, 8)
        ok, cnt, pk = P.legal(G)
        m = len(rest)
        ok = ok.reshape(-1, m).all(axis=1)
        pk = pk.reshape(-1, m, 4)
        cnt = cnt.reshape(-1, m)
        un = acc[surv] | np.bitwise_or.reduce(pk, axis=1)
        sizes = tot[surv] + cnt.sum(axis=1)
        disjoint = np.bitwise_count(un).astype(np.int64).sum(axis=1) == sizes
        ok &= disjoint & ((un & covp[None, :]) == 0).all(axis=1)
        acc[surv] = un
        surv = surv[ok]
    return surv, acc[surv]


def search(P, record_dim8=None):
    st = {"nodes": [0] * 9, "p1": [0] * 9, "triesD": [0] * 9, "tries": [0] * 9, "full": [], "allcov_low": 0,
          "cores": [0] * 9, "fp": 0}

    def rec(span, cov, d, parent_p1=False):
        a = next((v for v in P.order if not cov >> v & 1), None)
        st["nodes"][d] += 1
        p1 = a is None or a in P.D
        if p1:
            st["p1"][d] += 1
            if not parent_p1:
                st["cores"][d] += 1
                st["fp"] = (st["fp"] + fnv(space_key(span))) % (1 << 64)
        if a is None:
            if d == 8:
                st["full"].append(span.copy())
            else:
                st["allcov_low"] += 1
            return
        if d == 8:
            return
        parts = []
        if a in P.radof:
            parts.append(P.F4[P.radof[a]])
        else:
            for li in P.lines_at[a]:
                if not P.lmask[li] & cov:
                    parts.append(P.line_forms(li))
        for f, r in P.oracle:
            if a in r and not r <= P.D and not mask_of(r) & cov:
                parts.append(f[None, :])
        if not parts:
            return
        cand = np.concatenate(parts)
        st["tries"][d] += len(cand)
        if a in P.D:
            st["triesD"][d] += len(cand)
        surv, acc = children(P, span, cov, cand)
        for f, rm in zip(cand[surv], acc):
            rec(np.concatenate([span, span ^ f[None, :]]), cov | to_int(rm), d + 1, p1)

    rec(np.zeros((1, 8), dtype=np.uint8), 0, 0)
    return st


# ---------------------------------------------------------------- core fingerprints (as seven_check's)
def code28(rows):
    """Bit p of the code = B(e_i, e_j) for the p-th pair i < j (i ascending, then j)."""
    c = 0
    for p, (i, j) in enumerate(PAIRS):
        if int(rows[i]) >> j & 1:
            c |= 1 << p
    return c


def space_key(span):
    """The span's fully reduced row-echelon basis (pivots from bit 27 down), sorted, as '%07x ' each."""
    piv = {}
    for r in (code28(f) for f in span):
        for c in sorted(piv, reverse=True):
            if r >> c & 1:
                r ^= piv[c]
        if r:
            c = r.bit_length() - 1
            for c2 in piv:
                if piv[c2] >> c & 1:
                    piv[c2] ^= r
            piv[c] = r
    return "".join(f"{x:07x} " for x in sorted(piv.values()))


def fnv(s):
    h = 1469598103934665603
    for ch in s.encode():
        h = ((h ^ ch) * 1099511628211) % (1 << 64)
    return h


# ---------------------------------------------------------------- inputs
def graph(N):
    """graph(N) = {(x, N x)}: x the low nibble, N a 4x4 matrix with row r in bits 4r..4r+3."""
    pts = set()
    for x in range(1, 16):
        y = 0
        for r in range(4):
            if bin((N >> (4 * r)) & 15 & x).count("1") & 1:
                y |= 1 << r
        pts.add(x | y << 4)
    return pts


def config_members(path, idx):
    for line in open(path):
        if line.startswith("CONF") and int(line.split()[1]) == idx:
            Ns = [int(t, 16) for t in line.split(":")[1].split()]
            assert len(Ns) == 4
            return [set(range(1, 16)), {x << 4 for x in range(1, 16)}, graph(0x8421)] + [graph(N) for N in Ns]
    raise SystemExit(f"no configuration {idx}")


def mode_search(path, idxs):
    for idx in idxs:
        t0 = time.time()
        P = Problem(config_members(path, idx))
        assert len(P.D) == 150 and len(P.lines) == 2150 and all(len(P.lines_at[v]) == 43 for v in P.D)
        st = search(P)
        print(f"config {idx}: nodes {st['nodes']} p1 {st['p1']} triesD {st['triesD']} cores {st['cores']} fp {st['fp']:016x} full {len(st['full'])} "
              f"allcovered-below-8 {st['allcov_low']} [{time.time() - t0:.0f} s]", flush=True)


def lut_forms(F):
    """Alternating forms B_i(e_j, e_k) = bit i of F(e_j + e_k) + F(e_j) + F(e_k) + F(0) (rows per component)."""
    out = []
    for i in range(8):
        rows = [0] * 8
        for j in range(8):
            for k in range(8):
                if j != k and (F[(1 << j) ^ (1 << k)] ^ F[1 << j] ^ F[1 << k] ^ F[0]) >> i & 1:
                    rows[j] |= 1 << k
        out.append(rows)
    return np.array(out, dtype=np.uint8)


def all_solids_random(rng, avoid, need):
    """`need` random solids, pairwise disjoint and disjoint from the point set `avoid`."""
    chosen = []
    used = set(avoid)
    tries = 0
    while len(chosen) < need:
        tries += 1
        if tries > 200000:
            return None
        b = []
        s = {0}
        while len(b) < 4:
            v = rng.randrange(1, 256)
            if v in s:
                continue
            ns = s | {x ^ v for x in s}
            if ns & used:
                break
            b.append(v)
            s = ns
        if len(b) < 4:
            continue
        S = s - {0}
        chosen.append(S)
        used |= S
    return chosen


def mode_control(lut, first, n, seed):
    rng = random.Random(seed)
    done = found = skipped = 0
    for ln, line in enumerate(open(lut), 1):
        if ln < first or done >= n:
            continue
        vals = [int(t) for t in re.findall(r"[0-9]+", line)]
        if len(vals) != 256:
            continue
        B = lut_forms(vals)
        T = np.zeros((1, 8), dtype=np.uint8)
        for b in B:
            T = np.concatenate([T, T ^ b[None, :]])
        cnt, pk = radicals(T[1:])
        if (cnt == 255).any() or not set(np.unique(cnt)) <= {0, 3, 15}:
            skipped += 1
            continue
        r4 = [frozenset(v for v in range(1, 256) if to_int(pk[i]) >> v & 1) for i in np.flatnonzero(cnt == 15)]
        if len(r4) < 2:
            continue
        cov4 = set().union(*r4)
        extra = all_solids_random(rng, cov4, 7 - len(r4))
        if extra is None:
            skipped += 1
            continue
        members = r4 + extra
        memset = set().union(*members)
        oracle = []
        for i in np.flatnonzero(cnt == 3):
            r = frozenset(v for v in range(1, 256) if to_int(pk[i]) >> v & 1)
            if not r.isdisjoint(memset):
                oracle.append((T[1 + i], r))
        t0 = time.time()
        P = Problem(members, oracle)
        st = search(P)
        hit = any(set(map(tuple, s)) == set(map(tuple, T)) for s in st["full"])
        done += 1
        found += hit
        print(f"control line {ln}: k = {len(r4)}, oracle {len(oracle)}, nodes {st['nodes']}, full spans {len(st['full'])}, "
              f"T re-found: {'YES' if hit else 'NO'} [{time.time() - t0:.0f} s]", flush=True)
    print(f"controls: {done} run, {found} re-found, {skipped} skipped (not quadratic APN-shaped or no room)", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "search":
        mode_search(sys.argv[2], [int(x) for x in sys.argv[3:]])
    elif sys.argv[1] == "control":
        mode_control(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]))
