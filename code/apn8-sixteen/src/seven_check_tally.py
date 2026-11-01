"""Tallies and comparisons for seven_check's search outputs (search B, 7 rank-4 components).

Subagent (Claude Opus 5.5) for Rowen, 2026-10-31.
  python apn/seven_check_tally.py prog <configs.txt|-> <n> <prog files...>
      every configuration 0..n-1 exactly once (no duplicates, none missing); totals by dimension, by type
  python apn/seven_check_tally.py orders <runs/check7>
      per configuration: cores by dimension and fingerprint equal in asc, desc, sweep (and rand7 if present); the asc and
      desc core files equal as sets, configuration by configuration
  python apn/seven_check_tally.py numpy <runs/check7>
      numpy-<i>.out against asc's R line for i: nodes, p1, triesD by dimension, cores and fingerprint
R line: R id, nodes[9], tries[9], p1[9], bentext, full, triesD[9], cores[9], fp (hex), swFound[9], swChecks, swSols.

This folder's copy (made for it by Claude Opus 5.5) reads the files here, from the bundle root:
  python src/seven_check_tally.py prog data/configs7.txt 3731 logs/check7-asc-R.txt      (desc, sweep, phase1)
  python src/seven_check_tally.py orders logs [CORES]
  python src/seven_check_tally.py numpy logs
  python src/seven_check_tally.py reps logs data/configs7.txt
`prog` is unchanged; the R files here are the R lines of the runs' .prog files (their P1 lines aren't included).
`orders` and `numpy` read logs/check7-<run>-R.txt and logs/check7-numpy-output.txt instead of runs/check7/'s .prog and
numpy-<i>.out files; `orders` compares the core files (asc-r*.cores, desc-r*.cores: ~13 MB, not in this folder) only
when given a folder holding them. `reps` reads logs/check7-ctl-*.txt and checks its list against
logs/check7-reps-unreached.txt instead of writing it. Otherwise unchanged.
"""
import glob
import os
import re
import sys

FIELDS = ["nodes", "tries", "p1"]


def parse_r(line):
    t = line.split()
    v = t[1:]
    r = {"id": int(v[0])}
    x = [int(a) for a in v[1:39]]
    r["nodes"], r["tries"], r["p1"] = x[0:9], x[9:18], x[18:27]
    r["bentext"], r["full"] = x[27], x[28]
    r["triesD"] = x[29:38]
    r["cores"] = [int(a) for a in v[39:48]]
    r["fp"] = v[48]
    w = [int(a) for a in v[49:60]]
    r["swFound"], r["swChecks"], r["swSols"] = w[0:9], w[9], w[10]
    return r


def read_progs(files):
    recs, dup, sols = {}, 0, []
    for f in files:
        for line in open(f):
            if line.startswith("R "):
                r = parse_r(line)
                if r["id"] in recs:
                    dup += 1
                    if recs[r["id"]] != r:
                        print(f"  config {r['id']}: two different records")
                    continue
                recs[r["id"]] = r
            elif line.startswith("SOLUTION"):
                sols.append(line.strip())
    return recs, dup, sols


def types_of(conf):
    ty = {}
    if conf != "-":
        for line in open(conf):
            if line.startswith("CONF"):
                p = line.split()
                ty[int(p[1])] = p[2]
    return ty


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def prog(conf, n, files):
    n = int(n)
    recs, dup, sols = read_progs(files)
    missing = [i for i in range(n) if i not in recs]
    extra = [i for i in recs if not 0 <= i < n]
    print(f"{len(recs)} configurations, {dup} duplicate records, {len(missing)} missing {missing[:10]}, {len(extra)} out of range")
    ty = types_of(conf)
    groups = {"all": list(recs.values())}
    for r in recs.values():
        groups.setdefault(ty.get(r["id"], "?"), []).append(r)
    for g, rs in groups.items():
        tot = {k: [0] * 9 for k in ("nodes", "p1", "triesD", "cores", "swFound")}
        full = sw = chk = bent = 0
        for r in rs:
            for k in tot:
                tot[k] = add(tot[k], r[k])
            full += r["full"]; sw += r["swSols"]; chk += r["swChecks"]; bent += r["bentext"]
        print(f"[{g}] {len(rs)} configurations: nodes {tot['nodes']}; p1 {tot['p1']}; cores {tot['cores']}; "
              f"tries in D {tot['triesD']} (total {sum(tot['triesD']):,}); sweep legal spans {tot['swFound']}, "
              f"checks {chk:,}; bent extensions {bent}; full {full}; sweep solutions {sw}")
    print(f"SOLUTION lines: {len(sols)}")
    for s in sols[:20]:
        print("  " + s)


def read_cores(files):
    # K <conf> <dim> <dim codes>. A run killed mid-configuration can leave a truncated line that the rerun's first line
    # is glued onto: scan the text for complete records (exactly dim codes) instead of splitting lines.
    by, bad = {}, 0
    rec = re.compile(r"K (\d+) (\d+) ((?:[0-9a-f]{7} )+)")
    for f in files:
        for m in rec.finditer(open(f).read()):
            d, codes = int(m[2]), m[3].split()
            if len(codes) != d:
                bad += 1
                continue
            by.setdefault(int(m[1]), set()).add((d, " ".join(codes)))
    if bad:
        print(f"  ({bad} truncated core records skipped)")
    return by


def orders(d, cores=None):
    """Cores (by dimension, fingerprint) of every run against the reference (phase1 if present, else asc), each on the
    configurations both have; the sweep's spans against asc; the asc and desc core files as sets (if given)."""
    runs = {}
    for name in ("asc", "desc", "sweep", "rand7", "phase1"):
        fs = sorted(glob.glob(os.path.join(d, f"check7-{name}-R.txt")))
        if fs:
            runs[name] = read_progs(fs)[0]
    refn = "phase1" if "phase1" in runs else "asc"
    ref = runs[refn]
    print(f"runs: " + ", ".join(f"{k} {len(v)}" for k, v in runs.items()) + f"; reference {refn}")
    for name, r in runs.items():
        if name == refn:
            continue
        ids = sorted(set(r) & set(ref))
        bad = [i for i in ids if r[i]["cores"] != ref[i]["cores"] or r[i]["fp"] != ref[i]["fp"]]
        print(f"  {name}: cores by dimension and fingerprint on {len(ids)} configurations: {len(bad)} disagreements {bad[:10]}")
    if "sweep" in runs and "asc" in runs:  # the sweep follows asc's branching point with K_a
        ids = sorted(set(runs["sweep"]) & set(runs["asc"]))
        sb = 0
        for i in ids:
            a, w = runs["asc"][i], runs["sweep"][i]
            if a["nodes"] != add(w["nodes"], w["swFound"]) or a["full"] != w["swSols"]:
                sb += 1
                if sb <= 10:
                    print(f"  config {i}: asc nodes {a['nodes']} vs sweep {w['nodes']} + {w['swFound']}")
        print(f"  asc nodes by dimension = sweep's phase-1 nodes + sweep's legal spans: {sb} disagreements over {len(ids)}")
    if cores is None:
        print("  core files: not in this folder (give the folder holding asc-r*.cores and desc-r*.cores)")
        return
    ca, cd = read_cores(glob.glob(os.path.join(cores, "asc-r*.cores"))), read_cores(glob.glob(os.path.join(cores, "desc-r*.cores")))
    ids = set(runs.get("asc", {})) & set(runs.get("desc", {}))
    diff = [i for i in ids if ca.get(i, set()) != cd.get(i, set())]
    ncore = sum(len(ca.get(i, ())) for i in ids)
    cnt_bad = [i for i in ids if len(ca.get(i, ())) != sum(runs["asc"][i]["cores"])]
    print(f"  core files: {ncore} distinct cores (asc) over the {len(ids)} configurations both finished; asc and desc sets "
          f"differ on {len(diff)} {sorted(diff)[:10]}; asc file count != R-line count on {len(cnt_bad)}")


def numpy(d):
    asc = read_progs([os.path.join(d, "check7-asc-R.txt")])[0]
    n = bad = pend = 0
    tot = 0
    for f in [l for l in open(os.path.join(d, "check7-numpy-output.txt")) if not l.startswith("#")]:
        txt = f
        m = re.match(r"config (\d+): nodes \[([^]]*)\] p1 \[([^]]*)\] triesD \[([^]]*)\] cores \[([^]]*)\] fp (\w+) full (\d+)", txt)
        if not m:
            print(f"  {f}: unfinished")
            continue
        i = int(m[1])
        got = {"nodes": [int(x) for x in m[2].split(",")], "p1": [int(x) for x in m[3].split(",")],
               "triesD": [int(x) for x in m[4].split(",")], "cores": [int(x) for x in m[5].split(",")], "fp": m[6]}
        r = asc.get(i)
        if r is None:
            pend += 1
            continue
        n += 1
        tot += sum(got["triesD"])
        diffs = [k for k in got if got[k] != r[k]]
        # configurations from 1527 on ran with the candidate caches (cores of dimension <= 5), which change the tries in
        # D but not the spans visited
        if diffs == ["triesD"] and i >= 1527 and sum(r["cores"][:6]):
            print(f"  config {i}: tries in D differ by design (cache), everything else equal")
            diffs = []
        if diffs or int(m[7]):
            bad += 1
            print(f"  config {i}: differs in {diffs}, full {m[7]}")
    print(f"numpy: {n} configurations compared, {bad} differ, {pend} waiting for asc; tries in D over those compared {tot:,}")


def reps(d, conf):
    """Configurations reached by the ctlbatch controls (REPS lines of ctl-*.txt); checks the list against
    check7-reps-unreached.txt (count, then the unreached indices), which ctlrot read."""
    n = sum(1 for line in open(conf) if line.startswith("CONF"))
    ty = types_of(conf)
    hit = set()
    for f in glob.glob(os.path.join(d, "check7-ctl-*.txt")):
        if f.endswith("-output.txt"):
            continue
        for line in open(f):
            if line.startswith("REPS"):
                hit |= {int(x.split(":")[0]) for x in line.split()[2:]}
    un = [i for i in range(n) if i not in hit]
    by = {}
    for i in un:
        by[ty[i]] = by.get(ty[i], 0) + 1
    print(f"configurations reached by the control batches: {len(hit)} of {n}; unreached {len(un)} {by}")
    same = open(os.path.join(d, "check7-reps-unreached.txt")).read() == f"{len(un)}\n{' '.join(map(str, un))}\n"
    print(f"the same list as check7-reps-unreached.txt: {same}")


if __name__ == "__main__":
    {"prog": lambda: prog(sys.argv[2], sys.argv[3], sys.argv[4:]), "orders": lambda: orders(*sys.argv[2:4]),
     "numpy": lambda: numpy(sys.argv[2]), "reps": lambda: reps(sys.argv[2], sys.argv[3])}[sys.argv[1]]()
