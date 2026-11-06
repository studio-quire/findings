"""Section 7 of six-check.md: search B against search A, configuration by configuration (comparison only, after the verdict).

  python src/six_check_vs_a.py [logs dir] [data dir]      (defaults: logs/, data/ beside src/)
Reads logs/check6-mapa-output.txt (A's configuration -> ours, from `six_check mapa`), data/configs6.txt (type, |Aut|),
logs/six-near.txt (A's near misses per configuration, from its 8 outputs, and our core sets' digests; src/six_near.py)
and our R lines in logs/check6-p1asc-R.txt and -p1desc-R.txt (cores by dimension). Where A's configuration line and ours
carry the same six solids (A's NEAR spans and our core spans are then in one frame), the sets of spans are compared.
Rowen (Claude Opus 5.5), 2026-11-04; seven_check_vs_a.py for six.

This folder's copy (made for it by Claude Opus 5.5): the original (`apn/six_check_vs_a.py [sets]`) read A's NEAR lines
from its 8 cover outputs and our P1 lines from p1asc.prog (about 790 MB and 507 MB, not included), and compared the spans
as sets in memory. Here they come through six_near.py's per-configuration counts and digests, so the set comparison is
of digests and can't tell the 116 configurations with different frames apart from a partial overlap: its line says
"others", where the original's output (logs/check6-vs-a-output.txt) has "disjoint 116, partly shared 0".
"""
import collections, os, re, sys

here = os.path.dirname(os.path.abspath(__file__))
L = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "logs")
D = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "data")

A2B, autA = {}, {}
for l in open(os.path.join(L, "check6-mapa-output.txt")):
    m = re.match(r"A (p6-\d+-a(\d+)-w(\d+)) -> (-?\d+) .* canon aut (\d+)", l)
    if m:
        A2B[m[1]] = int(m[4]); autA[m[1]] = (int(m[2]), int(m[3]), int(m[5]))
typ, autB = {}, {}
for l in open(os.path.join(D, "configs6.txt")):
    m = re.match(r"CONF (\d+) (\S+) .* aut (\d+)", l)
    if m:
        typ[int(m[1])] = m[2]; autB[int(m[1])] = int(m[3])

print(f"mapa: {len(A2B)} of A's, {len(set(A2B.values()))} distinct of ours, none unmapped: {min(A2B.values()) >= 0}")
bad_aut = [a for a, (an, w, ac) in autA.items() if not an == ac == autB[A2B[a]]]
print(f"|Aut| (A's name, our canon of A's rep, our classification) equal on {len(autA) - len(bad_aut)} of {len(autA)}")
wt = collections.Counter((w, typ[A2B[a]]) for a, (an, w, ac) in autA.items())
print("w vs type:", dict(sorted(wt.items())))
print(f"  one to one: {len({w for w, t in wt}) == len({t for w, t in wt}) == len(wt)}")

anear = collections.Counter()
adig, bdig, nspan = {}, {}, {}
for l in open(os.path.join(L, "six-near.txt")):
    t = l.split()
    if t[0] == "A":
        c = A2B[t[1]]
        adig[c], nspan[c] = t[-1], int(t[2])
        for d, n in enumerate(t[3:12]):
            if int(n):
                anear[(c, d)] += int(n)
    elif t[0] == "B-p1asc":
        bdig[int(t[1])] = t[-1]
print(f"A: {sum(anear.values())} NEAR lines, by dim "
      f"{dict(sorted(collections.Counter(d for (c, d) in anear.elements()).items()))}; on {len({c for c, d in anear})} configurations")

for order in ("asc", "desc"):
    bcore = collections.Counter()
    nrec = 0
    for l in open(os.path.join(L, f"check6-p1{order}-R.txt")):
        if l.startswith("R "):
            t = l.split(); nrec += 1
            vals = [int(x) for x in t[2:-12]]
            for d, c in enumerate(vals[38:47]):
                if c:
                    bcore[(int(t[1]), d)] += c
    keys = set(anear) | set(bcore)
    off = sorted(k for k in keys if anear[k] != bcore[k])
    print(f"== p1{order}: {nrec} R records; B {sum(bcore.values())} cores on {len({c for c, d in bcore})} configurations; "
          f"(configuration, dim) pairs {len(keys)}, unequal {len(off)} {off[:10]}")
    for d in range(3, 8):
        ca = {c for (c, dd) in anear if dd == d}; cb = {c for (c, dd) in bcore if dd == d}
        if ca or cb:
            print(f"  dim {d}: A {sum(anear[(c, d)] for c in ca)} on {len(ca)}, B {sum(bcore[(c, d)] for c in cb)} on {len(cb)}, "
                  f"same configurations {ca == cb}")

eq = [c for c in adig if adig[c] == bdig.get(c)]
print(f"sets (p1asc), by digest: A {sum(nspan.values())} spans; configurations with equal sets {len(eq)} "
      f"({sum(nspan[c] for c in eq)} spans), others {len(adig) - len(eq)}")
