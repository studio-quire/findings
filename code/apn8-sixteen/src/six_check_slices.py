"""Resumable driver for six_check (search B for 6 rank-4 components), for the share runner.

Subagent (Claude Opus 5.5) for Rowen, 2026-11-03; adapted from seven_check_slices.py. Every stage writes per-slice or
per-range outputs and skips finished ones, so an eviction loses at most the slice or range in progress.
  python apn/six_check_slices.py count3     3-cliques (and 1..4) of V5824, two ways -> count3.out
  python apn/six_check_slices.py parta      Part A: one slice per class
  python apn/six_check_slices.py merge      merge of all Part A slices -> orbits.txt, merge.out
  python apn/six_check_slices.py configs    hole sets, covers, source families -> configs.txt, configs.out
  python apn/six_check_slices.py selftest
  python apn/six_check_slices.py search <name> [ranges]   see SEARCHES
  python apn/six_check_slices.py lin <name> [ranges]      the linear search (six_check lin), see LINS
Outputs in runs/check6.
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "runs", "bin")
EXE = os.path.join(BIN, "six_check.exe")
OUT = os.path.join(ROOT, "runs", "check6")
R = os.path.join(ROOT, "runs")
SPLIT = {0: 1, 1: 1, 2: 1, 3: 1, 4: 1}


def done(path):
    try:
        with open(path, "rb") as f:
            f.seek(0, 2)
            n = f.tell()
            f.seek(max(0, n - 16))
            return f.read().strip().endswith(b"END")
    except OSError:
        return False


def run_logged(log, args, mode="w", exe=EXE):
    with open(os.path.join(OUT, log), mode) as lf:
        r = subprocess.run([exe] + args, stdout=lf, stderr=subprocess.STDOUT)
    print(f"{args[0]} exit {r.returncode}", flush=True)
    return r.returncode


def count3():
    sys.exit(run_logged("count3.out", ["count3"]))


def parts():
    return [os.path.join(OUT, f"parta-c{c}-s{s}.txt") for c, ns in SPLIT.items() for s in range(ns)]


def parta():
    for c, ns in SPLIT.items():
        for s in range(ns):
            txt = os.path.join(OUT, f"parta-c{c}-s{s}.txt")
            if done(txt):
                continue
            mx = os.path.join(OUT, f"maxc-c{c}-s{s}.txt")
            rc = run_logged(f"parta-c{c}-s{s}.out", ["parta", str(c), txt, mx, str(s), str(ns)])
            if rc != 0 or not done(txt):
                print(f"class {c} slice {s}: exit {rc}", flush=True)
                sys.exit(1)
    print("parta: all slices done", flush=True)


def merge():
    assert all(done(p) for p in parts()), "Part A incomplete"
    sys.exit(run_logged("merge.out", ["merge", os.path.join(OUT, "orbits.txt"), os.path.join(R, "check7", "configs.txt")] + parts()))


def configs():
    sys.exit(run_logged("configs.out", ["configs", os.path.join(OUT, "orbits.txt"), os.path.join(R, "check7", "configs.txt"),
                                        os.path.join(R, "check8", "configs.txt"), os.path.join(R, "check14", "spreads.txt"),
                                        os.path.join(OUT, "configs.txt")]))


def selftest():
    sys.exit(run_logged("selftest.out", ["selftest", os.path.join(OUT, "configs.txt")]))


def prog_ids(path, tag="R "):
    ids = set()
    if os.path.exists(path):
        for line in open(path):
            if line.startswith(tag):
                ids.add(int(line.split()[1]))
    return ids


def nconfs():
    return sum(1 for line in open(os.path.join(OUT, "configs.txt")) if line.startswith("CONF"))


def by_ids(tag, mode, inp, ids, extra, post=(), exe=EXE, rtag="R "):
    # one prog file per run name; configurations in the given order, one process per configuration (so an eviction
    # loses at most one configuration, and the order of `ids` decides what is done first)
    pf = os.path.join(OUT, f"{tag}.prog")
    for i in ids:
        if i in prog_ids(pf, rtag):
            continue
        p = [a.replace("{t}", tag) for a in post]
        rc = run_logged(f"{tag}.out", [mode, inp, pf] + extra + [str(i), str(i)] + p, "a", exe)
        if rc or i not in prog_ids(pf, rtag):
            print(f"{tag}: configuration {i}: exit {rc}, incomplete", flush=True)
            sys.exit(1)
    print(f"{tag}: all {len(ids)} done", flush=True)


SEARCHES = {
    # the off-D phase alone in two orders (cores by dimension and fingerprint, written to <t>.cores)
    "p1asc": (["asc"], ["phase1", "cores", os.path.join(OUT, "{t}.cores")]),
    "p1desc": (["desc"], ["phase1", "cores", os.path.join(OUT, "{t}.cores")]),
    "asc": (["asc"], []),
}


def id_list(spec):
    # "all", "a-b", or a file of ids
    n = nconfs()
    if spec in (None, "all"):
        return list(range(n))
    if "-" in spec and spec.replace("-", "").isdigit():
        a, b = spec.split("-")
        return list(range(int(a), int(b) + 1))
    return [int(x) for x in open(os.path.join(OUT, spec)).read().split()]


def search(name, spec=None, exe=None):
    extra, post = SEARCHES[name]
    by_ids(name, "search", os.path.join(OUT, "configs.txt"), id_list(spec), extra, post, exe or EXE)


LINS = {
    "lin": ([], ["cores", os.path.join(OUT, "{t}.cores")]),
    "linp1": ([], ["phase1", "cores", os.path.join(OUT, "{t}.cores")]),
    "lindesc": ([], ["desc", "cores", os.path.join(OUT, "{t}.cores")]),
}


EXE_LIN = os.path.join(BIN, "six_check_lin.exe")


def lin(name, spec=None):
    extra, post = LINS[name]
    by_ids(name, "lin", os.path.join(OUT, "configs.txt"), id_list(spec), extra, post,
           os.path.join(BIN, "six_check_lin2.exe") if name == "lindesc" else EXE_LIN, "L ")


# controls (six_check_ctl.exe): name, LUT file, choices per path, ctlbatch options (tokens joined by '+'), max functions.
# lin = through the linear search; nond / nond2 / nond3 = in-D controls (the D phase from a subspace of T of
# codimension 1 / 2 / 3); no nond = unguided: the representative's whole search (every core with its D phase).
EXE_CTL = os.path.join(BIN, "six_check_ctl.exe")
# from the fix for cores that are already a whole space (the first lfull batch, kept as ctl-lfull-bl22-nofix.*, failed on it)
EXE_CTL2 = os.path.join(BIN, "six_check_ctl2.exe")
CTL = [("lind-bl22", "apn_8bit.txt", 3, "nond+lin", 10 ** 6), ("lind2-bl22", "apn_8bit.txt", 3, "nond2+lin", 10 ** 6),
       ("lind3-bl22", "apn_8bit.txt", 3, "nond3+lin", 10 ** 6), ("ind3-bl22", "apn_8bit.txt", 3, "nond3", 10 ** 6),
       ("ind2-bl22", "apn_8bit.txt", 3, "nond2", 10 ** 6),
       ("lind3-new", "new_apns.txt", 1, "nond3+lin", 10 ** 6), ("ind3-new", "new_apns.txt", 1, "nond3", 10 ** 6),
       ("lfull-bl22", "apn_8bit.txt", 1, "lin", 20), ("full-bl22", "apn_8bit.txt", 1, "", 5)]


def controls(*names):
    # a batch isn't resumable: one with fewer than two SUMMARY lines is deleted and rerun
    worst = 0
    for name, lut, nch, opt, maxn in CTL:
        if names and name not in names:
            continue
        txt = os.path.join(OUT, f"ctl-{name}.txt")
        if os.path.exists(txt) and sum(1 for line in open(txt) if line.startswith("SUMMARY")) == 2:
            continue
        if os.path.exists(txt):
            os.remove(txt)
        rc = run_logged(f"ctl-{name}.out", ["ctlbatch", os.path.join(OUT, "configs.txt"), os.path.join(R, "known", lut),
                                            "2", str(maxn), str(nch), txt] + ([opt] if opt else []), exe=EXE_CTL2 if name.startswith("lfull") else EXE_CTL)
        worst = max(worst, rc)
    sys.exit(worst)


# planted controls: plantings file, LUT, search (lin or pb), how many (from the start of the file), one process per
# planting so an eviction loses one; results appended to ctlplant-<name>.txt (FOUND / NOTFOUND / FAIL lines)
PLANTS = {"c2new-lin": ("plants-new-c2.txt", "new_apns.txt", "lin", 1000), "c2bl-lin": ("plants-bl22-c2.txt", "apn_8bit.txt", "lin", 1000),
          "c1bl-lin": ("plants-bl22-c1.txt", "apn_8bit.txt", "lin", 40),
          "c1new-lin": ("plants-new-c1.txt", "new_apns.txt", "lin", 1000),
          "c2new-pb": ("plants-new-c2.txt", "new_apns.txt", "pb", 1000), "c2bl-pb": ("plants-bl22-c2.txt", "apn_8bit.txt", "pb", 1000),
          "c1bl-pb": ("plants-bl22-c1.txt", "apn_8bit.txt", "pb", 10)}


def ctlplant(name, count=None):
    pf, lut, how, n = PLANTS[name]
    n = int(count) if count else n
    plants = [l for l in open(os.path.join(OUT, pf)) if l.startswith("PLANT")]
    txt = os.path.join(OUT, f"ctlplant-{name}.txt")
    done_ = set()
    if os.path.exists(txt):
        for l in open(txt):
            t = l.split()
            if len(t) > 2 and t[0] in ("FOUND", "NOTFOUND", "FAIL") and t[1] == "plant":
                done_.add(int(t[2].rstrip(":")))
    worst = 0
    for i in range(min(n, len(plants))):
        if i in done_:
            continue
        rc = run_logged(f"ctlplant-{name}.out", ["ctlplant", os.path.join(OUT, "configs.txt"), os.path.join(R, "known", lut),
                                                 os.path.join(OUT, pf), txt, str(i), "1"] + (["lin"] if how == "lin" else []), "a", EXE_CTL)
        worst = max(worst, rc)
    sys.exit(worst)


def ctlrot_lin(pool="40"):
    # one unguided control per configuration, all 341, through lin: BL22's first `pool` functions with k >= 3 in
    # rotation, T's radicals embedded into the representative's members (seven's ctlrot, unguided), chunks of 31
    n = nconfs()
    exe = os.path.join(BIN, "six_check_ctl3.exe")
    worst = 0
    for k in range(0, n, 31):
        p = os.path.join(OUT, f"ctlrot-lin-{k // 31}.out")
        reps = [str(i) for i in range(k, min(n, k + 31))]
        if os.path.exists(p) and sum(1 for line in open(p) if line.startswith("rep ")) == len(reps):
            continue
        with open(p, "w") as lf:
            r = subprocess.run([exe, "ctlrot", os.path.join(OUT, "configs.txt"), os.path.join(R, "known", "apn_8bit.txt"), "1", pool] + reps + ["lin"],
                               stdout=lf, stderr=subprocess.STDOUT)
        print(f"ctlrot-lin chunk {k // 31}: exit {r.returncode}", flush=True)
        worst = max(worst, r.returncode)
    sys.exit(worst)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    a = sys.argv[2:]
    {"count3": count3, "parta": parta, "merge": merge, "configs": configs, "selftest": selftest,
     "search": lambda: search(*a), "lin": lambda: lin(*a), "controls": lambda: controls(*a),
     "ctlplant": lambda: ctlplant(*a), "ctlrotlin": lambda: ctlrot_lin(*a)}[sys.argv[1]]()
