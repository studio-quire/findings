"""Compare cover_sample.cmd's runs: per frontier item, the tree's counts against the filter-free hole stage's.

  python apn/cover_sample_cmp.py runs/cover5/sample [splits: s6 s5 s4 by default]
Each ITEM line: configuration, the item's generators, nodes by dimension 1..8, solutions, full, non-APN. Rowen, 2026-11-04.
"""
import collections, re, sys

d = sys.argv[1]
bad = 0
for split in sys.argv[2:] or ("s6", "s5", "s4"):
    a = sorted(l.strip() for l in open(f"{d}/tree-{split}.out") if l.startswith("ITEM "))
    b = sorted(l.strip() for l in open(f"{d}/old-{split}.out") if l.startswith("ITEM "))
    sa, sb = set(a), set(b)
    only_a, only_b = sa - sb, sb - sa
    confs = {l.split()[1] for l in a}
    nodes = collections.Counter()
    for l in a:
        v = l.split(" nodes ")[1].split(" sols ")[0].split()
        for k, x in enumerate(v):
            nodes[k + 1] += int(x)
    sols = sum(int(re.search(r" sols (\d+)", l)[1]) for l in a)
    done = [re.search(r"^done: .*", open(f"{d}/{m}-{split}.out").read(), re.M) for m in ("tree", "old")]
    print(f"{split}: tree {len(a)} items, old {len(b)} items, on {len(confs)} configurations; lines only in tree {len(only_a)}, only in old {len(only_b)}; "
          f"nodes by dimension {dict(sorted(nodes.items()))}, solutions {sols}; both finished: {all(done)}")
    for l in sorted(only_a)[:3] + sorted(only_b)[:3]:
        print("  DIFF", l[:200])
    bad += len(only_a) + len(only_b) + (not all(done)) + (len(a) != len(b))
print("identical" if not bad else f"{bad} differences")
sys.exit(1 if bad else 0)
