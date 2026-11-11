"""Search A's 22 configurations for five against search B's 22 classes, by contents: five_vs_b's B-rule count of admissible
cores on each of A's configurations should be exactly one of B's per-class counts (five-check.md §5), with the same number of
dimension 4 and |Gamma| (fivecheck_classes gamma 5) equal to A's |Aut P| tag. Rowen (Claude Opus 5.5), 2026-11-09.

python apn/five_classmatch.py runs/i12/vs-b.txt apn/five-check.md runs/check5/bundle/gamma5.txt
"""
import re
import sys

vsb, check, gamma = sys.argv[1:4]
B = {}
for line in open(check, encoding='utf-8'):
    m = re.match(r'\| (\d+) \| (\d+) \| (\d+) \|', line)
    if m:
        B[int(m.group(2))] = (int(m.group(1)), int(m.group(3)))
gam = {}
for line in open(gamma, encoding='utf-8'):
    m = re.match(r'cfg (\d+) \|Gamma\| = (\d+)', line)
    if m:
        gam[int(m.group(1))] = int(m.group(2))
assert len(B) == 22 and len(gam) == 22
hit, bad = set(), 0
for line in open(vsb, encoding='utf-8'):
    m = re.match(r'(p5-\d+-a(\d+)-w\d+) B-rule (\d+) dim4 (\d+)', line)
    if not m:
        continue
    tag, a, n, d4 = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
    c, bd4 = B.get(n, (None, None))
    ok = c is not None and bd4 == d4 and gam.get(c) == a and c not in hit
    hit.add(c)
    bad += not ok
    print(f'{tag} -> class {c} |Gamma| {gam.get(c)} {"OK" if ok else "MISMATCH"}')
print(f'B counts distinct {len(B)}; classes hit {len(hit)}; mismatches {bad}')
sys.exit(1 if bad or len(hit) != 22 else 0)
