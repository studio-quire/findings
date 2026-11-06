"""Search A's tree against search B's, per configuration and dimension (review-six.md item 5). From the bundle's folder:

  python src/cover_percfg_cmp.py logs
A's nodes come from logs/cover-percfg-output.txt (cover_percfg.py: the published cover5t, one configuration at a time), B's from
lin less x4 (spans holding a seventh solid, which A and the point branching reject outright). A and B share a frame on a
configuration where their near-miss digests (six-near.txt, labels A and B-p1asc) are equal. B's first phase alone is p1asc, so
lin less x4 less p1asc is B's hole stage; R files are the point branching's whole searches, which use no filter.
Totals must be the published ones. Rowen (Claude Opus 5.5), 2026-11-05.
"""
import re
import sys

logs = sys.argv[1]
DIMS = range(1, 9)


def counts(path, tag):
    """First nine numbers after the configuration id, on the R lines of a .prog-style file."""
    out = {}
    for l in open(f"{logs}/{path}"):
        if l.startswith(tag + " "):
            t = l.split()
            out[int(t[1])] = list(map(int, t[2:11]))
    return out


def lin(path):
    out = {}
    for l in open(f"{logs}/{path}"):
        if l.startswith("L "):
            t = l.split()
            v = list(map(int, t[3:48]))
            out[int(t[1])] = [a - b for a, b in zip(v[0:9], v[18:27])]
    return out


A2B = {}
for l in open(f"{logs}/check6-mapa-output.txt"):
    m = re.match(r"A (p6-\S+) -> (\d+) ", l)
    if m:
        A2B[m[1]] = int(m[2])

digest = {"A": {}, "B-p1asc": {}}
near = {}
for l in open(f"{logs}/six-near.txt"):
    t = l.split()
    if t[0] in digest:
        digest[t[0]][t[1]] = t[-1]
    if t[0] == "A":
        near[A2B[t[1]]] = list(map(int, t[3:12]))
shared = {A2B[a] for a, d in digest["A"].items() if digest["B-p1asc"].get(str(A2B[a])) == d}

A, cand, sent, sols = {}, 0, 0, 0
for part in open(f"{logs}/cover-percfg-output.txt").read().split("== ")[1:]:
    tag = part.split()[0]
    done = re.search(r"^done: .*solutions (\d+), full \(dim < 8\) (\d+), dim 8 non-APN (\d+), bad seeds (\d+)", part, re.M)
    sols += sum(map(int, done.groups()))
    A[A2B[tag]] = [0] + list(map(int, re.search(r"^nodes by span dimension 1\.\.8: (.*)", part, re.M)[1].split()))
    m = re.search(r"^hole-stage candidates (\d+), passing the filter \(sent to push\) (\d+)", part, re.M)
    cand += int(m[1]); sent += int(m[2])

L, Ldesc = lin("check6-lin-prog.txt"), lin("check6-lindesc-prog.txt")
P1, P1desc = counts("check6-p1asc-R.txt", "R"), counts("check6-p1desc-R.txt", "R")
R = {**counts("check6-asc-R.txt", "R"), **counts("check6-asc340-R.txt", "R")}
hole = {c: [a - b for a, b in zip(L[c], P1[c])] for c in L}

print(f"configurations {len(A)} of {len(A2B)}, {len(shared)} with a shared frame")
print(f"totals: nodes {sum(sum(v) for v in A.values())}, by dimension 1..8 {[sum(v[d] for v in A.values()) for d in DIMS]}, "
      f"solutions and full and non-APN and bad seeds {sols}, candidates {cand}, sent to push {sent}")

print("\nA against B (lin less x4):")
for name, group in (("shared frames", shared), ("frames differ", set(A) - shared)):
    print(f"{name}: {len(group)} configurations")
    for d in DIMS:
        a, b = sum(A[c][d] for c in group), sum(L[c][d] for c in group)
        print(f"  dim {d}: equal on {sum(A[c][d] == L[c][d] for c in group)}; A {a}, B {b}, difference {a - b}")

print("\nwhere the phases sit (B, all 341):")
print(f"  first phase: nodes above dim 6 on {sum(any(v[7:]) for v in P1.values())} configurations; "
      f"dim 6 equal to the 6-dimensional near misses on {sum(P1[c][6] == near[c][6] for c in P1)}")
print(f"  hole stage by dimension 1..8: {[sum(v[d] for v in hole.values()) for d in DIMS]}; "
      f"at dim 5 only on {[(c, v[5]) for c, v in hole.items() if v[5]]}, the configurations with 4-dimensional near misses {[c for c in near if near[c][4]]}")
print("  point order, B ascending against B descending:")
for d in DIMS:
    print(f"    dim {d}: first phase equal on {sum(P1[c][d] == P1desc[c][d] for c in P1)} of 341; "
          f"whole tree on {sum(L[c][d] == Ldesc[c][d] for c in L)} of 341, {sum(L[c][d] == Ldesc[c][d] for c in shared)} of the shared")

below = [sum(A[c][6] - near[c][6] for c in g) for g in (shared, A)]
print(f"\nnodes below the near misses on shared frames: dim 6 {below[0]} of {below[1]} ({below[0] / below[1]:.1%}), "
      f"dim 7 {sum(A[c][7] for c in shared)} of {sum(A[c][7] for c in A)} ({sum(A[c][7] for c in shared) / sum(A[c][7] for c in A):.1%})")

print(f"\nA against the point branching, {len(R)} configurations (lin less x4 equal to it on all: {all(R[c] == L[c] for c in R)}): A less it, dims 1..8")
for c in sorted(R):
    print(f"  {c}{' shared' if c in shared else ''}{f' ({near[c][4]} near misses of dim 4)' if near[c][4] else ''}: {[A[c][d] - R[c][d] for d in DIMS]}")
