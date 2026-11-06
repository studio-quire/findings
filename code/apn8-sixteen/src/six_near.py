"""Near misses for i = 11, per configuration: counts by dimension and a digest of the set of spans.

Written for this folder (Claude Opus 5.5, 2026-11-04). Search A's NEAR lines (7,941,871 for six, about 790 MB of
cover output) and search B's cores (the P1 lines of its .prog files, as large) are too big to include, so this reduces
each run to one line per configuration:
  <label> <configuration> <spans> <by dim 0..8> <digest>
where the digest is the first 16 hex digits of SHA-256 over the configuration's spans, each taken as its canonical
basis (the reduced echelon basis of its 28-bit coefficient vectors, pivots at the highest bits, sorted) packed into one
integer (element i at bit 28 i) and written in hex, one per line, the spans sorted: equal digests mean equal sets.
Configurations are search A's names for NEAR lines and search B's indices for P1 lines (B's .prog files give a
configuration's P1 lines after its R line). A span repeated within a configuration is counted once in <spans> and
reported.
  python src/six_near.py OUT LABEL=FILE[,FILE...] [LABEL=...]
  python src/six_near.py compare logs/six-near.txt      labels whose configurations should agree, compared
"""
import collections
import hashlib
import sys


def canon(vs):
    out = []
    for v in vs:
        for r in out:
            if (v >> (r.bit_length() - 1)) & 1:
                v ^= r
        if not v:
            raise ValueError("dependent basis")
        p = v.bit_length() - 1
        out = [r ^ v if (r >> p) & 1 else r for r in out]
        out.append(v)
    return sum(v << (28 * i) for i, v in enumerate(sorted(out)))


def spans(files):
    per = collections.defaultdict(set)
    rep = 0
    for f in files:
        cur = None
        for line in open(f):
            if line.startswith("NEAR "):
                t = line.split()
                key = t[1]
            elif line.startswith("R "):
                cur = int(line.split()[1])
                continue
            elif line.startswith("P1 "):
                t = line.split()
                key = cur
            else:
                continue
            b = canon(int(x, 16) for x in t[t.index("basis") + 1:])
            if b in per[key]:
                rep += 1
            per[key].add(b)
    return per, rep


def main(out, specs):
    with open(out, "w", newline="\n") as fo:
        for spec in specs:
            label, files = spec.split("=", 1)
            per, rep = spans(files.split(","))
            for c in sorted(per, key=lambda k: (str(type(k)), k)):
                s = per[c]
                dims = collections.Counter((b.bit_length() + 27) // 28 for b in s)
                h = hashlib.sha256("".join(f"{b:x}\n" for b in sorted(s)).encode()).hexdigest()[:16]
                fo.write(f"{label} {c} {len(s)} {' '.join(str(dims[d]) for d in range(9))} {h}\n")
            print(f"{label}: {sum(len(s) for s in per.values())} spans on {len(per)} configurations, {rep} repeated",
                  flush=True)


def compare(path, pairs=(("A-cover4", "A"), ("B-p1desc", "B-p1asc"), ("A-prof-gray40", "A-prof-tree40"), ("A-prof-old40", "A-prof-tree40"))):
    runs = collections.defaultdict(dict)
    for line in open(path):
        t = line.split()
        runs[t[0]][t[1]] = t[2:]
    for a, b in pairs:
        same = sum(runs[a][c] == runs[b].get(c) for c in runs[a])
        print(f"{a} against {b}: {len(runs[a])} configurations, counts and digest equal on {same}")


if __name__ == "__main__":
    if sys.argv[1] == "compare":
        compare(sys.argv[2])
    else:
        main(sys.argv[1], sys.argv[2:])
