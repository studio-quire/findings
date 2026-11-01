"""Section 7 of seven-check.md: search B against search A, configuration by configuration (comparison only, after the verdict).

  python src/seven_check_vs_a.py [logs dir] [data dir]      (defaults: logs/, data/ beside src/)
Reads logs/check7-mapa-output.txt (A's configuration -> ours, from `seven_check mapa`), data/configs7.txt (type, |Aut|),
A's NEAR lines in logs/seven-output-{c0..c3,c5..c7,c4x,p2222-s5}.txt (the nine outputs of the final tally; the cancelled
first run of c4 is not read), and our R lines in logs/check7-{asc,desc}-R.txt (cores by dimension, fields 39..47).
Rowen (Claude Opus 5.5), 2026-11-01; eight_check_vs_a.py for seven. This folder's copy (made for it by Claude Opus 5.5)
reads the files here (the original read Rowen's working tree: runs/check7/mapa.out and configs.txt,
runs/i10/{cover,cover-b,cover-c,split}/*.out, runs/check7/<order>-r*.prog); otherwise unchanged.
"""
import collections, glob, os, re, sys

sys.path.insert(0, os.path.dirname(__file__))
from seven_check_tally import read_progs

HERE = os.path.dirname(os.path.abspath(__file__))
LOGS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "logs")
DATA = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "..", "data")

A2B, autA, autMap = {}, {}, {}
for l in open(os.path.join(LOGS, "check7-mapa-output.txt")):
    m = re.match(r"A (p7-\d+-a(\d+)-w(\d+)) -> (-?\d+) .* canon aut (\d+)", l)
    if m:
        A2B[m[1]] = int(m[4]); autA[m[1]] = (int(m[2]), int(m[3]), int(m[5]))
typ, autB = {}, {}
for l in open(os.path.join(DATA, "configs7.txt")):
    m = re.match(r"CONF (\d+) (\S+) .* aut (\d+)", l)
    if m:
        typ[int(m[1])] = m[2]; autB[int(m[1])] = int(m[3])

print(f"mapa: {len(A2B)} of A's, {len(set(A2B.values()))} distinct of ours, none unmapped: {min(A2B.values()) >= 0}")
bad_aut = [a for a, (an, w, ac) in autA.items() if not an == ac == autB[A2B[a]]]
print(f"|Aut| (A's name, our canon of A's rep, our classification) equal on {len(autA) - len(bad_aut)} of {len(autA)}")
wt = collections.Counter((w, typ[A2B[a]]) for a, (an, w, ac) in autA.items())
print("w vs type:", dict(sorted(wt.items())))
print(f"  one to one: {len({w for w, t in wt}) == len({t for w, t in wt}) == len(wt)}")

outs = [os.path.join(LOGS, f"seven-output-{c}.txt") for c in ("c0", "c1", "c2", "c3", "c5", "c6", "c7", "c4x", "p2222-s5")]
anear, aseen = collections.Counter(), collections.Counter()
for f in outs:
    for l in open(f):
        m = re.match(r"NEAR (\S+) dim (\d+) .* basis (.*)$", l.strip())
        if m:
            anear[(A2B[m[1]], int(m[2]))] += 1
            aseen[(m[1], m[3])] += 1
print(f"A: {len(outs)} outputs, {sum(anear.values())} NEAR lines, {sum(c > 1 for c in aseen.values())} repeated "
      f"(config, basis); by dim {dict(sorted(collections.Counter(d for (c, d) in anear.elements()).items()))}; "
      f"on {len({c for c, d in anear})} configurations")

for order in ("asc", "desc"):
    recs, dup, sols = read_progs([os.path.join(LOGS, f"check7-{order}-R.txt")])
    bcore = collections.Counter({(i, d): r["cores"][d] for i, r in recs.items() for d in range(9) if r["cores"][d]})
    keys = set(anear) | set(bcore)
    off = sorted(k for k in keys if anear[k] != bcore[k])
    print(f"== {order}: {len(recs)} records ({dup} repeated, {len(sols)} solutions); B {sum(bcore.values())} cores on {len({c for c, d in bcore})} configurations; "
          f"(configuration, dim) pairs {len(keys)}, unequal {len(off)} {off[:10]}")
    for d in range(4, 8):
        ca = {c for (c, dd) in anear if dd == d}; cb = {c for (c, dd) in bcore if dd == d}
        print(f"  dim {d}: A {sum(anear[(c, d)] for c in ca)} on {len(ca)}, B {sum(bcore[(c, d)] for c in cb)} on {len(cb)}, "
              f"same configurations {ca == cb}")
