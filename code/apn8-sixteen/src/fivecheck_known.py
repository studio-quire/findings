"""Known quadratic APN functions -> S8 bases, in the sieve's form format (8 row bytes in a u64), by rank-4 count.

Reads the LUTs in runs/a/seeds-known-rare.txt, derives the 8 forms straight from F, checks APN by the definition
(all derivative images of size 128) and writes runs/check5/known.txt: one line per function, `n4 n2 n0 forms...`.
"""
import re, sys

def forms_of(F):
    out = []
    for b in range(8):
        f = 0
        for i in range(8):
            row = 0
            for j in range(8):
                x, y = 1 << i, 1 << j
                v = F[x ^ y] ^ F[x] ^ F[y] ^ F[0] if i != j else 0
                if v >> b & 1:
                    row |= 1 << j
            f |= row << (8 * i)
        out.append(f)
    return out

def rad_size(f):
    rows = [(f >> (8 * i)) & 255 for i in range(8)]
    n = 0
    for v in range(256):
        if all(bin(r & v).count('1') % 2 == 0 for r in rows):
            n += 1
    return n

def is_apn(F):
    for a in range(1, 256):
        if len({F[x] ^ F[x ^ a] for x in range(256)}) != 128:
            return False
    return True

def main():
    lines = open(sys.argv[1] if len(sys.argv) > 1 else 'runs/a/seeds-known-rare.txt').read().splitlines()
    tally, out = {}, []
    for ln in lines:
        vals = [int(x) for x in re.findall(r'\d+', ln)]
        if len(vals) != 256:
            print('skipped line of length', len(vals)); continue
        F = vals
        if not is_apn(F):
            print('not APN'); continue
        fs = forms_of(F)
        dist = {}
        for m in range(1, 256):
            f = 0
            for i in range(8):
                if m >> i & 1:
                    f ^= fs[i]
            s = rad_size(f)
            dist[s] = dist.get(s, 0) + 1
        key = (dist.get(16, 0), dist.get(4, 0), dist.get(1, 0), dist.get(64, 0))
        tally[key] = tally.get(key, 0) + 1
        out.append('%d %d %d %s' % (key[0], key[1], key[2], ' '.join('%016x' % x for x in fs)))
    open('runs/check5/known.txt', 'w').write('\n'.join(out) + '\n')
    for k in sorted(tally):
        print('rank4 %d rank6 %d rank8 %d rank2 %d: %d functions' % (k[0], k[1], k[2], k[3], tally[k]))

main()
