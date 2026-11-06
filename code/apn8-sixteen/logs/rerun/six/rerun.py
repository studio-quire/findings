"""Re-runs for six from the bundle's own files (cwd = bundle root), with fresh builds of src/cover.cpp, pspread3.cpp and
six_check.cpp (Build line). Each step writes OUT/<name>.log ending in '# exit <code>'; finished steps are skipped, so
the script resumes. Times in OUT/times.txt.
  python logs/rerun/six/rerun.py [step names]      (groups: a = search A's steps, b = search B's)
Environment: BIN (the builds, named cover.exe, pspread3.exe, six_check.exe; default this folder), OUT (default
BIN/out), KNOWN (folder with apn_8bit.txt).
Claude Opus 5.5, for this folder, 2026-11-04."""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
BUNDLE = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BIN = os.environ.get("BIN", HERE)
O = os.environ.get("OUT", os.path.join(BIN, "out"))
K = os.environ.get("KNOWN", "known")
E, P, C = (os.path.join(BIN, f"{n}.exe") for n in ("six_check", "pspread3", "cover"))
os.makedirs(O, exist_ok=True)
os.chdir(BUNDLE)
o = lambda f: os.path.join(O, f)
CONFIGS = "data/configs6.txt"
BL22 = os.path.join(K, "apn_8bit.txt")
PROF = ["-1", "5", "--only", "p6-101-a24-w11", "--frange", "0", "40"]

A = [
    ("p3-enum6", [P, "enum", "6", o("p6.txt"), "--ckpt", o("p6.ckpt"), "--wcheck", "1000"]),
    ("p3-selftest6", [P, "selftest", "6", "data/p6.txt"]),
    ("p3-holes6", [P, "holes", "6", "data/p6.txt", o("holes6.txt")]),
    ("cover-p11", [C, "data/p11.txt"]),
    ("cover-p10", [C, "data/p10.txt"]),
    ("cover-prof-tree40", [C, o("c2.txt")] + PROF),
    ("cover-prof-gray40", [C, o("c2.txt")] + PROF + ["--gray"]),
    ("cover-c0", [C, o("c0.txt"), "-1", "5"]),
    ("cover-c1", [C, o("c1.txt"), "-1", "5"]),
    ("cover-c2", [C, o("c2.txt"), "-1", "5"]),
    # the positive controls through the rewritten hole stage (cover has had k4, k3, k2 since thirteen, which leave the
    # hole stage little to do); k1 leaves it far more, and was cancelled as too slow on the old hole stage
    ("cover-ctl-k4", [C, "data/ctl-12631-k4.txt"]),
    ("cover-ctl-k3", [C, "data/ctl-12631-k3.txt"]),
    ("cover-ctl-k2", [C, "data/ctl-12631-k2.txt"]),
    ("prep-ctl-k1", [sys.executable, "src/cover_prep.py", "ctl", o("ctl-12631-k1.txt"), BL22, "12631", "1", "1"]),
    ("cover-ctl-k1", [C, o("ctl-12631-k1.txt")]),
]
parts = [o(f"parta-c{c}-s0.txt") for c in range(5)]
B = [("s6-count3", [E, "count3"])]
B += [(f"s6-parta-c{c}-s0", [E, "parta", str(c), o(f"parta-c{c}-s0.txt"), o(f"maxc-c{c}-s0.txt"), "0", "1"]) for c in range(5)]
B += [
    ("s6-merge", [E, "merge", o("orbits.txt"), "data/configs7.txt"] + parts),
    ("s6-configs", [E, "configs", o("orbits.txt"), "data/configs7.txt", "data/configs8.txt", "data/spreads.txt", o("configs6.txt")]),
    ("s6-mapa", [E, "mapa", CONFIGS, "data/p6.txt"]),
    ("s6-selftest", [E, "selftest", CONFIGS]),
    ("s6-linselftest", [E, "linselftest", CONFIGS, "16"]),
    ("s6-plantfind-c2", [E, "plantfind", BL22, "4", "100000", "2", o("plants-bl22-c2.txt")]),
    ("s6-plantfind-c1", [E, "plantfind", BL22, "2", "100000", "1", o("plants-bl22-c1.txt")]),
    # the list search B's controls read, logs/check6-plants-bl22-c2.txt, has 5 plantings with k >= 3; the logged job
    # (the step above, k >= 4) wrote 3 and a later run, not logged, overwrote them: this is that run, as far as can be told
    ("s6-plantfind-c2k3", [E, "plantfind", BL22, "3", "100000000", "2", o("plants-bl22-c2k3.txt")]),
]
# the off-D phase on two stretches in each order, one process per range (the originals ran one per configuration)
for name, order, lo, hi in (("p1asc", "asc", 0, 40), ("p1asc", "asc", 300, 340), ("p1desc", "desc", 100, 140)):
    B.append((f"s6-{name}-{lo}-{hi}", [E, "search", CONFIGS, o(f"{name}.prog"), order, str(lo), str(hi), "phase1", "cores", o(f"{name}.cores")]))
# lin on a sample: the first of its type, a dimension-4 core, the point-branching re-run's configuration, the lone ext10
for i in (0, 24, 57, 190, 248, 340):
    B.append((f"s6-lin-{i}", [E, "lin", CONFIGS, o("lin.prog"), str(i), str(i), "cores", o("lin.cores")]))
for i in (0, 340):
    B.append((f"s6-lindesc-{i}", [E, "lin", CONFIGS, o("lindesc.prog"), str(i), str(i), "desc", "cores", o("lindesc.cores")]))
B += [
    ("s6-asc-248", [E, "search", CONFIGS, o("asc.prog"), "asc", "248", "248"]),
    ("s6-ctl-lind-bl22", [E, "ctlbatch", CONFIGS, BL22, "2", "1000000", "3", o("ctl-lind-bl22.txt"), "nond+lin"]),
    ("s6-ctlplant-c2bl-0", [E, "ctlplant", CONFIGS, BL22, "logs/check6-plants-bl22-c2.txt", o("ctlplant-c2bl-lin.txt"), "0", "1", "lin"]),
    ("s6-ctlplant-c2bl-3", [E, "ctlplant", CONFIGS, BL22, "logs/check6-plants-bl22-c2.txt", o("ctlplant-c2bl-lin.txt"), "3", "1", "lin"]),
    ("s6-ctlrot-lin", [E, "ctlrot", CONFIGS, BL22, "1", "40", "0", "340", "lin"]),
]
groups = {"a": [s[0] for s in A], "b": [s[0] for s in B]}

# search A's chunks, cut from data/p6.txt as the run cut them (chunks of 43 lines)
lines = open("data/p6.txt", newline="").read().splitlines(True)
for k in range(3):
    with open(o(f"c{k}.txt"), "w", newline="") as f:
        f.writelines(lines[k * 43:(k + 1) * 43])

only = set()
for a in sys.argv[1:]:
    only |= set(groups.get(a, [a]))
for name, cmd in A + B:
    if only and name not in only:
        continue
    log = o(name + ".log")
    if os.path.exists(log):
        tail = open(log, errors="replace").read().rstrip().splitlines()[-1:]
        if tail and tail[0].startswith("# exit"):
            continue
    t0 = time.time()
    with open(log, "w") as lf:
        rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode
    with open(log, "a") as lf:
        lf.write(f"# exit {rc}\n")
    with open(o("times.txt"), "a") as tf:
        tf.write(f"{name}: exit {rc}, {time.time() - t0:.0f} s\n")
    print(f"{name}: exit {rc}, {time.time() - t0:.0f} s", flush=True)
