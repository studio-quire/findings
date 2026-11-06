"""Exact masses for six_check's classification (search B, 6 rank-4 components), in Python integers.

Subagent (Claude Opus 5.5) for Rowen, 2026-11-03, adapted from seven_check_mass.py.
  python src/six_check_mass.py PATH/TO/fix-6.txt [data dir] [logs dir]      (defaults: data/, logs/ beside src/)
Reads data/configs6.txt (CONF: |Aut P|; HOLE: |Aut H|, solid inside, number of covers), logs/check6-merge-output.txt
(weighted 3-cliques), logs/check6-count3-output.txt (3-cliques by plain enumeration) and Langevin-Hou's fix-6.txt
(download it from https://langevin.univ-tln.fr/project/spread/psp.html; not included), and checks:
  (a) hole classes 341, of them 0 containing no solid (Kurz, Table 9.3, r = 6);
  (b) the number of hole classes per |Aut H| = Langevin-Hou's per fixator order, and per order their qte
      (= classes * |GL(8,2)| / order), and the total sum_H |GL(8,2)| / |Aut H| = their Total. Since the classes found
      are pairwise inequivalent, equality of the total means no class is missing;
  (c) labelled 6-partial spreads three ways: sum_H |GL|/|Aut H| * covers(H) (covers counted from the hole set alone),
      sum_P |GL|/|Aut P| over the configurations, and W3 * |GL| / (120 * 20160) with W3 = the weighted 3-cliques of
      Part A and = count3's plain count.

This folder's copy (made for it by Claude Opus 5.5) reads the files here; the original, `six_check_mass.py
<runs/check6>`, read configs.txt, merge.out and count3.out from that folder and fix-6.txt from runs/papers/langevin/.
Otherwise unchanged.
"""
import os
import re
import sys

GL8 = 1
for k in range(8):
    GL8 *= 2 ** 8 - 2 ** k


def main(fix, data, logs):
    confs, holes = [], {}
    for line in open(os.path.join(data, "configs6.txt")):
        if line.startswith("CONF"):
            m = re.match(r"CONF (\d+) (\S+) hole (\d+) aut (\d+)", line)
            confs.append((int(m[1]), m[2], int(m[3]), int(m[4])))
        elif line.startswith("HOLE"):
            m = re.match(r"HOLE (\d+) aut (\d+) solid (\d) covers (\d+)", line)
            holes[int(m[1])] = (int(m[2]), int(m[3]), int(m[4]))
    print(f"|GL(8,2)| = {GL8}")
    nmax = sum(1 for a, s, c in holes.values() if not s)
    print(f"(a) hole classes {len(holes)}, containing no solid {nmax}   [Kurz r = 6: 341, 0]  "
          f"{'OK' if (len(holes), nmax) == (341, 0) else 'MISMATCH'}")
    lh, lhq, lhtot = {}, {}, None
    txt = open(fix).read()
    for m in re.finditer(r"fix=(\d+)\s+number of class (\d+)\s+qte=(\d+)", txt):
        lh[int(m[1])] = int(m[2])
        lhq[int(m[1])] = int(m[3])
    lhtot = int(re.search(r"Total=(\d+)", txt)[1])
    assert int(re.search(r"Order of GL\(2,8\) is (\d+)", txt)[1]) == GL8
    ours = {}
    for a, s, c in holes.values():
        ours[a] = ours.get(a, 0) + 1
    print(f"    Langevin-Hou: {sum(lh.values())} classes over {len(lh)} fixator orders; ours: {sum(ours.values())} over {len(ours)}")
    bad = [(o, ours.get(o, 0), lh.get(o, 0)) for o in sorted(set(lh) | set(ours)) if ours.get(o, 0) != lh.get(o, 0)]
    qbad = [o for o in lh if lhq[o] != lh[o] * GL8 // o or (GL8 % o)]
    print(f"(b) classes per |Aut H| against Langevin-Hou: {len(bad)} orders differ {bad[:10]}; their qte = classes * |GL| / order: "
          f"{'OK' if not qbad else 'differs for ' + str(qbad)}")
    tot = sum(GL8 // a for a, s, c in holes.values())
    print(f"    sum_H |GL|/|Aut H| = {tot}; Langevin-Hou Total = {lhtot}  {'OK' if tot == lhtot else 'MISMATCH'}")
    by_cov = sum((GL8 // a) * c for a, s, c in holes.values())
    by_conf = sum(GL8 // aut for _, _, _, aut in confs)
    W4 = int(re.search(r"weighted 3-cliques (\d+)", open(os.path.join(logs, "check6-merge-output.txt")).read())[1])
    c4 = open(os.path.join(logs, "check6-count3-output.txt")).read()
    W4p = int(re.search(r" 3:(\d+)", c4)[1])
    by_w4 = W4 * GL8 // (120 * 20160)
    print(f"(c) labelled 6-partial spreads: from covers {by_cov}; from configurations {by_conf}; from Part A's W3 "
          f"({W4}) {by_w4}; count3's plain 3-clique count {W4p} {'= W3' if W4p == W4 else '!= W3 MISMATCH'}  "
          f"{'OK' if by_cov == by_conf == by_w4 and W4 * GL8 % (120 * 20160) == 0 and W4p == W4 else 'MISMATCH'}")
    types = {}
    for _, t, _, _ in confs:
        types[t] = types.get(t, 0) + 1
    print(f"configurations {len(confs)}: {types}")


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "..", "data"),
         sys.argv[3] if len(sys.argv) > 3 else os.path.join(here, "..", "logs"))
