"""Invariance check of the Gamma reduction on the dimension-4 sets: within each orbit, the sieve's survivor counts must agree.

Reads runs/check5/k4all (every tuple, sieve only) and k4orbit_cNN.txt. Also checks each orbit's representative has a finished
assembly in runs/check5/k4.
"""
import glob, re, os

total = orbits = mixed = 0
for path in sorted(glob.glob('runs/check5/k4orbit_c*.txt')):
    c = int(path[-6:-4])
    orb = [int(x) for x in open(path).read().split()]
    if not orb:
        continue
    by = {}
    for i, o in enumerate(orb):
        p = 'runs/check5/k4all/l_r5_c%02d_%06d.txt' % (c, i)
        t = open(p).read()
        assert t.rstrip().endswith('DONE'), p
        m = re.search(r'processed (\d)', t)
        assert m.group(1) == '1', p
        m = re.search(r' f0 (\d+) f1 (\d+) ', t)
        by.setdefault(o, set()).add((int(m.group(1)), int(m.group(2))))
        total += 1
    for o, v in by.items():
        orbits += 1
        if len(v) != 1:
            mixed += 1
            print('class', c, 'orbit', o, 'not constant:', sorted(v)[:4])
    nreps = len(by)
    for r in range(nreps):
        p = 'runs/check5/k4/l_r5_c%02d_%06d.txt' % (c, r)
        if not (os.path.exists(p) and open(p).read().rstrip().endswith('DONE')):
            print('representative unfinished:', p)
print('tuples', total, 'orbits', orbits, 'orbits with differing counts', mixed)
