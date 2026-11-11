"""The S8s of the solids-only runs (runs/check5/review/nou), by the plain definition, in Python, sharing nothing with the C++.

For each: all 255 forms, radicals by evaluation, how many points are covered 0, 1, 2+ times, and whether the 75 solid points are covered
exactly once. Also rebuilds the quadratic function from the forms' dual? No: APN of a quadratic F is the partition condition, re-derived.
Writes the tuples of the sets S that carry a T, for the reference.
"""
import sys, glob, re, collections

def rows(f):
    return [(f >> (8 * i)) & 255 for i in range(8)]

def radical(f):
    R = rows(f)
    comb = [0] * 256
    out = []
    for v in range(1, 256):
        b = (v & -v).bit_length() - 1
        comb[v] = comb[v & (v - 1)] ^ R[b]
        if comb[v] == 0:
            out.append(v)
    return out

forms = {}
for l in open('runs/check5/review/forms.txt'):
    c, i, d, x = l.split()
    forms[(int(c), int(x, 16))] = (int(i), int(d))

allrows = 0
apn = 0
stats = collections.Counter()
tuples = collections.Counter()
for p in sorted(glob.glob('runs/check5/review/nou/r5_c*.txt')):
    c = int(re.search(r'_c(\d+)_', p)[1])
    nfile = 0
    for l in open(p):
        if not l.startswith('S8 '):
            continue
        nfile += 1
        allrows += 1
        b = [int(x, 16) for x in l.split(':')[1].split()]
        assert len(b) == 8
        cover = [0] * 256
        dist = collections.Counter()
        seen = set()
        for m in range(1, 256):
            f = 0
            for i in range(8):
                if m >> i & 1:
                    f ^= b[i]
            assert f and f not in seen   # the eight forms are independent
            seen.add(f)
            r = radical(f)
            dist[len(r)] += 1
            for v in r:
                cover[v] += 1
        # solid points: union of the first five forms' radicals
        sol = set()
        for f in b[:5]:
            r = radical(f)
            assert len(r) == 15
            sol |= set(r)
        assert len(sol) == 75
        good = all(cover[v] == 1 for v in sol)
        unc = sum(1 for v in range(1, 256) if cover[v] == 0)
        multi = sum(1 for v in range(1, 256) if cover[v] > 1)
        stats[(c, 'solids once' if good else 'SOLIDS NOT ONCE', 'uncovered %d' % unc, 'multi %d' % multi, 'rad sizes ' + str(sorted(dist.items())))] += 1
        if unc == 0 and multi == 0:
            apn += 1
            print('APN!', l)
        tup = [None] * 5
        for f in b[:5]:
            i, d = forms[(c, f)]
            tup[i] = d
        tuples[(c,) + tuple(tup)] += 1
    print(p, nfile, 'S8')
for k, v in sorted(stats.items(), key=str):
    print(v, k)
print('S8 checked', allrows, 'truly APN', apn, 'distinct sets S carrying a T', len(tuples))
with open('runs/check5/review/nou_tuples.txt', 'w') as o:
    for t, n in sorted(tuples.items()):
        o.write(' '.join(map(str, t)) + '\n')
with open('runs/check5/review/nou_tuple_counts.txt', 'w') as o:
    for t, n in sorted(tuples.items()):
        o.write(' '.join(map(str, t)) + ' %d\n' % n)
