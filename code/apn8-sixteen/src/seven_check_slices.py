"""Resumable driver for seven_check (search B for 7 rank-4 components), for the share runner.

Subagent (Claude Opus 5.5) for Rowen, 2026-10-31; adapted from eight_check_slices.py. Every stage writes per-slice or
per-range outputs and skips finished ones, so an eviction loses at most the slice or range in progress.
  python apn/seven_check_slices.py count4     4-cliques (and 1..5) of V5824, two ways -> count4.out
  python apn/seven_check_slices.py parta      Part A: class 0 in 16 slices, 1 and 2 in 8, 3 in 4, 4 in 1
  python apn/seven_check_slices.py merge      merge of all Part A slices -> orbits.txt, merge.out
  python apn/seven_check_slices.py configs    hole sets, covers, source families -> configs.txt, configs.out
  python apn/seven_check_slices.py selftest
  python apn/seven_check_slices.py search [names]   asc, desc (cores logged), sweep (coset sweep D phase), rand7
  python apn/seven_check_slices.py cross <name>     cross-runs on inputs that don't use this classification (CROSS)
  python apn/seven_check_slices.py controls / ctlrot [pool] [dsweep] / numpy [nproc]
Outputs in runs/check7.
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE = os.path.join(ROOT, "runs", "bin", "seven_check.exe")
OUT = os.path.join(ROOT, "runs", "check7")
R = os.path.join(ROOT, "runs")
SPLIT = {0: 16, 1: 8, 2: 8, 3: 4, 4: 1}


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


def count4():
    sys.exit(run_logged("count4.out", ["count4"]))


def parts():
    return [os.path.join(OUT, f"parta-c{c}-s{s}.txt") for c, ns in SPLIT.items() for s in range(ns)]


def parta():
    for c, ns in SPLIT.items():
        for s in range(ns):
            txt = os.path.join(OUT, f"parta-c{c}-s{s}.txt")
            if done(txt):
                continue
            mx = os.path.join(OUT, f"maxc-c{c}-s{s}.txt")
            print(f"class {c} slice {s} of {ns}", flush=True)
            rc = run_logged(f"parta-c{c}-s{s}.out", ["parta", str(c), txt, mx, str(s), str(ns)])
            if rc != 0 or not done(txt):
                print(f"class {c} slice {s}: exit {rc}", flush=True)
                sys.exit(1)
    print("parta: all slices done", flush=True)


def merge():
    assert all(done(p) for p in parts()), "Part A incomplete"
    sys.exit(run_logged("merge.out", ["merge", os.path.join(OUT, "orbits.txt"), os.path.join(R, "check8", "configs.txt")] + parts()))


def configs():
    sys.exit(run_logged("configs.out", ["configs", os.path.join(OUT, "orbits.txt"), os.path.join(R, "check8", "configs.txt"),
                                        os.path.join(R, "check9", "configs.txt"), os.path.join(R, "check14", "spreads.txt"),
                                        os.path.join(OUT, "configs.txt")]))


def selftest():
    sys.exit(run_logged("selftest.out", ["selftest", os.path.join(OUT, "configs.txt")]))


def prog_ids(path):
    ids = set()
    if os.path.exists(path):
        for line in open(path):
            if line.startswith("R "):
                ids.add(int(line.split()[1]))
    return ids


def ranged(tag, mode, inp, n, nr, extra, post=(), only=None, exe=EXE):
    size = -(-n // nr)
    for k in range(nr):
        if only is not None and k not in only:
            continue
        lo, hi = k * size, min(n, (k + 1) * size) - 1
        if lo > hi:
            continue
        t = tag if nr == 1 else f"{tag}-r{k}"
        pf = os.path.join(OUT, f"{t}.prog")
        if len({i for i in prog_ids(pf) if lo <= i <= hi}) == hi - lo + 1:
            continue
        print(f"{t}: configurations {lo}..{hi}", flush=True)
        p = [a.replace("{t}", t) for a in post]
        rc = run_logged(f"{t}.out", [mode, inp, pf] + extra + [str(lo), str(hi)] + p, "a", exe)
        if rc or len({i for i in prog_ids(pf) if lo <= i <= hi}) != hi - lo + 1:
            print(f"{t}: exit {rc}, incomplete", flush=True)
            sys.exit(1)
    print(f"{tag}: all {n} done", flush=True)


# the sweep runs from seven_check_sw.exe (the same source, built while seven_check.exe was busy with asc and desc)
EXE_SW = os.path.join(ROOT, "runs", "bin", "seven_check_sw.exe")
SEARCHES = {
    "asc": (["asc"], ["cores", os.path.join(OUT, "{t}.cores")], EXE),
    "desc": (["desc"], ["cores", os.path.join(OUT, "{t}.cores")], EXE),
    "sweep": (["asc"], ["dsweep"], EXE_SW),
    "rand7": (["rand", "7"], [], EXE),
    # the off-D phase alone (cores by dimension and fingerprint, no D phase): cost estimates, and a fourth core count
    "phase1": (["desc"], ["phase1"], os.path.join(ROOT, "runs", "bin", "seven_check_p1.exe")),
}


def search(*names):
    C = os.path.join(OUT, "configs.txt")
    n = sum(1 for line in open(C) if line.startswith("CONF"))
    for name in names or ("asc", "desc", "sweep"):
        extra, post, exe = SEARCHES[name]
        ranged(name, "search", C, n, 8, extra, post, exe=exe)


# cross-runs on inputs that don't use this classification: mode, input, count, ranges, step (every step-th index)
CROSS = {
    "s8": ("search8", os.path.join("check8", "configs.txt"), 9372 * 8, 16),
}


def cross(name, *rs):
    mode, inp, n, nr = CROSS[name]
    ranged(name, mode, os.path.join(R, inp), n, nr, ["asc"], only={int(r) for r in rs} if rs else None)


def cross_sample(nblocks="10", seed="7"):
    # eight-check's 9372 x 8 deletions (no classification involved): nblocks random blocks of 24 consecutive family
    # indices (3 eight-configurations x their 8 deletions each), one prog file per block (s8-b<start>.prog)
    import random
    rng = random.Random(int(seed))
    starts = sorted(24 * rng.randrange(9372 * 8 // 24) for _ in range(int(nblocks)))
    inp = os.path.join(R, "check8", "configs.txt")
    for st in starts:
        pf = os.path.join(OUT, f"s8-b{st}.prog")
        if len({i for i in prog_ids(pf) if st <= i < st + 24}) == 24:
            continue
        rc = run_logged(f"s8-b{st}.out", ["search8", inp, pf, "asc", str(st), str(st + 23)], "a")
        if rc:
            sys.exit(rc)
    print("cross sample done:", starts, flush=True)


# controls: name, LUT file, choices per path, extra flag
# controls: name, LUT file, choices per path, options (ctlbatch tokens joined by '+'), max functions. All from
# seven_check_ctl.exe (the source with guided and two-step controls, built while the other exes were busy). "guided": only
# cores inside the planted T get a D phase; in these controls T's off-D points need T's rank-6 (oracle) forms, so the
# core is usually T itself and the D phase is tested by the in-D controls (nond: one D step from a hyperplane, nond2:
# two from codimension 2, with the point branching or the sweep). ind-bl22 ran first, unguided, from seven_check.exe
# (with a start at a hyperplane the only core is the start, so guided and unguided are the same search).
EXE_CTL = os.path.join(ROOT, "runs", "bin", "seven_check_ctl.exe")
CTL = [("ind-bl22", "apn_8bit.txt", 3, "nond", 1000000), ("ind2-bl22", "apn_8bit.txt", 3, "nond2+guided", 1000000),
       ("bl22", "apn_8bit.txt", 3, "guided", 1000000), ("blp22", "blp22_8bit.txt", 6, "guided", 1000000),
       ("yuperrin", "yuperrin_qam.txt", 6, "guided", 1000000), ("yuwangli", "yuwangli_qam.txt", 6, "guided", 1000000),
       ("new", "new_apns.txt", 1, "guided", 1000000), ("ind2-new", "new_apns.txt", 1, "nond2+guided", 1000000),
       ("swind-bl22", "apn_8bit.txt", 1, "nond+dsweep+guided", 1000000), ("swind2-bl22", "apn_8bit.txt", 1, "nond2+dsweep+guided", 1000000),
       ("full-bl22", "apn_8bit.txt", 1, "", 20)]


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
                                            "2", str(maxn), str(nch), txt] + ([opt] if opt else []), exe=EXE_CTL)
        worst = max(worst, rc)
    sys.exit(worst)


def ctlrot(pool="40", *flags):
    # one control per configuration (all of them, or those in reps-unreached.txt if it exists), BL22's first `pool`
    # functions with k >= 3 in rotation, parallel over the configurations, chunks of 400 (ctlrot-<k>.out)
    C = os.path.join(OUT, "configs.txt")
    n = sum(1 for line in open(C) if line.startswith("CONF"))
    reps = [str(i) for i in range(n)]
    un = os.path.join(OUT, "reps-unreached.txt")
    if os.path.exists(un) and "all" not in flags:
        reps = open(un).read().splitlines()[1].split()
    flags = [f for f in flags if f != "all"]
    tag = "ctlrot" + ("-sw" if "dsweep" in flags else "")
    worst = 0
    for k in range(0, len(reps), 400):
        p = os.path.join(OUT, f"{tag}-{k // 400}.out")
        if os.path.exists(p) and sum(1 for line in open(p) if line.startswith("rep ")) == len(reps[k:k + 400]):
            continue
        with open(p, "w") as lf:
            r = subprocess.run([EXE_CTL, "ctlrot", C, os.path.join(R, "known", "apn_8bit.txt"), "1", pool] + reps[k:k + 400] + list(flags),
                               stdout=lf, stderr=subprocess.STDOUT)
        print(f"{tag} chunk {k // 400}: exit {r.returncode}", flush=True)
        worst = max(worst, r.returncode)
    sys.exit(worst)


def numpy_sample(nproc="8", *ids):
    py = os.path.join(R, "venv", "Scripts", "python.exe")
    ids = [int(i) for i in ids] or [int(x) for x in open(os.path.join(OUT, "numpy-sample.txt")).read().split()]
    todo = [i for i in ids if not (os.path.exists(os.path.join(OUT, f"numpy-{i}.out")) and
                                    open(os.path.join(OUT, f"numpy-{i}.out")).read().startswith(f"config {i}:"))]
    print(f"numpy sample: {len(ids)} configurations, {len(todo)} to do", flush=True)
    import time
    running = []
    while todo or running:
        while todo and len(running) < int(nproc):
            i = todo.pop()
            f = open(os.path.join(OUT, f"numpy-{i}.out"), "w")
            running.append((subprocess.Popen([py, os.path.join(ROOT, "apn", "seven_check_numpy.py"), "search",
                                              os.path.join(OUT, "configs.txt"), str(i)], stdout=f, stderr=subprocess.STDOUT), f))
        for p, f in list(running):
            if p.poll() is not None:
                f.close()
                running.remove((p, f))
        time.sleep(5)
    print("numpy sample done", flush=True)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    a = sys.argv[2:]
    {"count4": count4, "parta": parta, "merge": merge, "configs": configs, "selftest": selftest,
     "search": lambda: search(*a), "cross": lambda: cross(*a), "crosssample": lambda: cross_sample(*a), "controls": lambda: controls(*a),
     "ctlrot": lambda: ctlrot(*a), "numpy": lambda: numpy_sample(*a)}[sys.argv[1]]()
