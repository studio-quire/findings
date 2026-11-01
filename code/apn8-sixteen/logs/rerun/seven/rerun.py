"""Re-runs for seven from the bundle's own files (cwd = bundle root), with fresh builds of src/cover.cpp, pspread3.cpp and
seven_check.cpp (Build line). Each step writes OUT/<name>.log ending in '# exit <code>'; finished steps are skipped, so
the script resumes. Times in OUT/times.txt.
  python logs/rerun/seven/rerun.py [step names]      (groups: a = search A's steps, b = search B's)
Environment: BIN (the builds, named cover<SUFFIX>.exe, pspread3<SUFFIX>.exe, seven_check<SUFFIX>.exe; default this
folder), OUT (default BIN/out), KNOWN (folder with apn_8bit.txt), VPY (a Python with numpy).
Claude Opus 5.5, for this folder, 2026-11-01."""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
BUNDLE = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BIN = os.environ.get("BIN", HERE)
SUF = os.environ.get("SUFFIX", "")
O = os.environ.get("OUT", os.path.join(BIN, "out"))
K = os.environ.get("KNOWN", "known")
VPY = os.environ.get("VPY", "python")
E, P, C = (os.path.join(BIN, f"{n}{SUF}.exe") for n in ("seven_check", "pspread3", "cover"))
os.makedirs(O, exist_ok=True)
os.chdir(BUNDLE)
o = lambda f: os.path.join(O, f)
NOCACHE = {"SEVEN_NOCACHE": "1"}  # asc and desc ran without the candidate caches up to configuration 1526

SPLIT = {0: 16, 1: 8, 2: 8, 3: 4, 4: 1}
parts = [o(f"parta-c{c}-s{s}.txt") for c, ns in SPLIT.items() for s in range(ns)]
A = [
    ("p3-enum7", [P, "enum", "7", o("p7.txt"), "--ckpt", o("p7.ckpt"), "--wcheck", "1000"], None),
    ("p3-selftest7", [P, "selftest", "7", "data/p7.txt"], None),
    ("p3-holes7", [P, "holes", "7", "data/p7.txt", o("holes7.txt")], None),
    ("cover-c2", [C, o("c2.txt")], None),
    ("cover-c2-split5", [C, o("c2.txt"), "-1", "5"], None),
]
B = [("s7-count4", [E, "count4"], None)]
B += [(f"s7-parta-c{c}-s{s}", [E, "parta", str(c), o(f"parta-c{c}-s{s}.txt"), o(f"maxc-c{c}-s{s}.txt"), str(s), str(ns)], None)
      for c, ns in SPLIT.items() for s in range(ns)]
B += [
    ("s7-merge", [E, "merge", o("orbits.txt"), "data/configs8.txt"] + parts, None),
    ("s7-configs", [E, "configs", o("orbits.txt"), "data/configs8.txt", "data/configs9.txt", "data/spreads.txt",
                    o("configs7.txt")], None),
    ("s7-mapa", [E, "mapa", "data/configs7.txt", "data/p7.txt"], None),
    ("s7-asc-0-10", [E, "search", "data/configs7.txt", o("asc.prog"), "asc", "0", "10", "cores", o("asc.cores")], NOCACHE),
    ("s7-asc-1000-1099", [E, "search", "data/configs7.txt", o("asc.prog"), "asc", "1000", "1099", "cores", o("asc.cores")], NOCACHE),
    ("s7-asc-3631-3730", [E, "search", "data/configs7.txt", o("asc.prog"), "asc", "3631", "3730", "cores", o("asc.cores")], None),
    ("s7-desc-1000-1099", [E, "search", "data/configs7.txt", o("desc.prog"), "desc", "1000", "1099", "cores", o("desc.cores")], NOCACHE),
    ("s7-sweep-2000-2099", [E, "search", "data/configs7.txt", o("sweep.prog"), "asc", "2000", "2099", "dsweep"], None),
    # the sweep's first ranges ran without its cache (third build): tries and checks change where a core has dimension
    # <= 5 (configuration 1), the spans visited don't
    ("s7-sweep-0-8", [E, "search", "data/configs7.txt", o("sweep.prog"), "asc", "0", "8", "dsweep"], None),
    ("s7-phase1-1000-1099", [E, "search", "data/configs7.txt", o("phase1.prog"), "desc", "1000", "1099", "phase1"], None),
    ("s7-s8-b4728", [E, "search8", "data/configs8.txt", o("s8-b4728.prog"), "asc", "4728", "4751"], None),
    ("s7-ctl-ind-bl22", [E, "ctlbatch", "data/configs7.txt", os.path.join(K, "apn_8bit.txt"), "2", "1000000", "3",
                         o("ctl-ind-bl22.txt"), "nond"], None),
    ("s7-ctlrot", [E, "ctlrot", "data/configs7.txt", os.path.join(K, "apn_8bit.txt"), "1", "40", "99", "102", "3707", "3729",
                   "guided"], None),
    ("s7-numpy-3", [VPY, "src/seven_check_numpy.py", "search", "data/configs7.txt", "3"], None),
    ("s7-numpy-ctl", [VPY, "src/seven_check_numpy.py", "control", os.path.join(K, "apn_8bit.txt"), "1", "1", "7"], None),
    ("s7-selftest", [E, "selftest", "data/configs7.txt"], None),
]
groups = {"a": [s[0] for s in A], "b": [s[0] for s in B]}

# search A's chunk c2, cut from data/p7.txt as the run cut it (chunks of 467 lines)
lines = open("data/p7.txt", newline="").read().splitlines(True)
with open(o("c2.txt"), "w", newline="") as f:
    f.writelines(lines[2 * 467:3 * 467])

only = set()
for a in sys.argv[1:]:
    only |= set(groups.get(a, [a]))
for name, cmd, env in A + B:
    if only and name not in only:
        continue
    log = o(name + ".log")
    if os.path.exists(log):
        tail = open(log, errors="replace").read().rstrip().splitlines()[-1:]
        if tail and tail[0].startswith("# exit"):
            continue
    t0 = time.time()
    with open(log, "w") as lf:
        rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, env=dict(os.environ, **(env or {}))).returncode
    with open(log, "a") as lf:
        lf.write(f"# exit {rc}\n")
    with open(o("times.txt"), "a") as tf:
        tf.write(f"{name}: exit {rc}, {time.time() - t0:.0f} s\n")
    print(f"{name}: exit {rc}, {time.time() - t0:.0f} s", flush=True)
