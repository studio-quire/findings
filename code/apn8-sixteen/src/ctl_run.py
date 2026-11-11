"""k = 4 positive control of the orbit reduction, end to end, with solutions that exist.

  python ctl_run.py <base tag, e.g. ctl4-73-k4> <full tag in ctl.txt> [max items]

The d4 control found Sol(C) under the core C of a known function, with no symmetry involved. The reduction says Sol(r) = Sol(C) o g
for r = C o g, g in Aut(P). So: from the pipeline's own group, near list and orbit list take r, the listed representative of C's
orbit, and an explicit g. Run the job's binary with the job's flags (--only, --p1all, --reps [r], --frange) on just the frontier items
under r that contain some S o g. Everything found must be an image, and every image must be found.
"""
import re, subprocess, sys, os

BIN = r"..\..\bin"
PAIRS = [(i, j) for i in range(8) for j in range(i + 1, 8)]


def rows(c):
    M = [0] * 8
    for k, (i, j) in enumerate(PAIRS):
        if c >> k & 1:
            M[i] |= 1 << j
            M[j] |= 1 << i
    return M


def app(g, x):
    y = 0
    for i in range(8):
        if x >> i & 1:
            y ^= g[i]
    return y


def act(c, g):
    """coefficient word of B(g x, g y)"""
    M = rows(c)
    o = 0
    w = [app(M, g[i]) for i in range(8)]   # B(g e_i, .) as a functional: rows of M combined by g e_i
    for k, (i, j) in enumerate(PAIRS):
        if bin(w[i] & g[j]).count("1") & 1:
            o |= 1 << k
    return o


def key(vs):
    piv = {}
    for v in vs:
        while v:
            h = v.bit_length() - 1
            if h in piv:
                v ^= piv[h]
            else:
                piv[h] = v
                break
    for h in sorted(piv):
        for c in piv:
            if c != h and piv[c] >> h & 1:
                piv[c] ^= piv[h]
    return tuple(sorted(piv.values(), reverse=True))


def span(vs):
    s = [0]
    for v in vs:
        s += [x ^ v for x in s]
    return s


def main():
    base, tag = sys.argv[1:3]
    maxitems = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 9
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    gens = []
    for l in open("ctl_aut.txt"):
        t = l.split()
        if t[0] == tag:
            gens = [tuple(int(p[2 * i:2 * i + 2], 16) for i in range(8)) for p in t[2:]]
    core = None
    for l in open("../../cover5/ctl/d4.txt"):
        t = l.split()
        if t[0] == base:
            core = [int(x, 16) for x in t[t.index("seed") + 1:]]
    sols = []
    for l in open("../../cover5/ctl/d4-split.out", errors="replace"):
        t = l.split()
        if l.startswith("SOLUTION") and t[1] == base:
            sols.append([int(x, 16) for x in t[t.index("basis") + 1:]])
    reps = {}
    for l in open(f"ctl/reps-{tag}.txt"):
        t = l.split()
        reps[tuple(int(x, 16) for x in t[1:1 + int(t[0])])] = l
    ident = tuple(1 << i for i in range(8))
    k0 = key(core)
    seen = {k0: ident}
    todo = [k0]
    r = k0 if k0 in reps else None
    while todo and r is None:
        nxt = []
        for k in todo:
            for h in gens:
                y = key([act(c, h) for c in k])
                if y not in seen:
                    g = seen[k]
                    seen[y] = tuple(app(g, h[i]) for i in range(8))   # (B o g) o h = B o (g h)
                    nxt.append(y)
                    if y in reps:
                        r = y
        todo = nxt
    g = seen[r]
    print(tag, "core", [hex(x) for x in k0], "rep", [hex(x) for x in r], "same" if r == k0 else "different", "g", bytes(g).hex(), "|Sol(C)|", len(sols), flush=True)
    assert key([act(c, g) for c in k0]) == r
    one = f"ctl/one-{tag}.txt"
    open(one, "w").write(reps[r])
    # the frontier under r
    env = dict(os.environ, RV_FRONT=f"ctl/front-{tag}.txt")
    p = subprocess.run([BIN + r"\rv5coverF.exe", "ctl.txt", "-1", "5", "--only", tag, "--p1all", "--reps", one, "--frange", "0", "0"],
                       capture_output=True, text=True, errors="replace", env=env)
    print(re.search(r"phase 1 done.*", p.stdout)[0], flush=True)
    print(re.search(r"phase 1 counts.*", p.stdout)[0], flush=True)
    front = [[int(x, 16) for x in l.split()] for l in open(f"ctl/front-{tag}.txt")]
    images = {}
    for S in sols:
        I = [act(c, g) for c in S]
        kk = key(I)
        els = set(span(I))
        assert len(els) == 256
        hit = [i for i, f in enumerate(front) if all(x in els for x in f)]
        assert len(hit) == 1, hit
        images[kk] = hit[0]
    print("distinct images", len(images), "in", len(set(images.values())), "frontier items of", len(front), flush=True)
    found = set()
    items = sorted(set(images.values()))[:maxitems]
    want = {k for k, i in images.items() if i in set(items)}
    for it in items:
        p = subprocess.run([BIN + r"\rv5cover.exe", "ctl.txt", "-1", "5", "--only", tag, "--p1all", "--reps", one, "--frange", str(it), str(it + 1)],
                           capture_output=True, text=True, errors="replace")
        assert p.returncode == 0
        d = re.search(r"^done: .*", p.stdout, re.M)[0]
        got = [key([int(x, 16) for x in l.split("basis")[1].split()]) for l in p.stdout.splitlines() if l.startswith("SOLUTION")]
        found |= set(got)
        print("item", it, "solutions", len(got), d[:100], flush=True)
    print("images wanted", len(want), "found", len(found), "wanted and not found", len(want - found), "found and not an image", len(found - set(images)), flush=True)
    print("RESULT", tag, "PASS" if want == found and r != k0 else "FAIL or trivial", flush=True)


main()
