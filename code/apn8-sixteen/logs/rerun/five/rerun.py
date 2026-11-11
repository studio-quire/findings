"""Re-runs for five from the bundle's own files, with fresh builds of its sources (Build line). Each step writes OUT/<name>.log
ending in '# exit <code>'; finished steps are skipped, so the script resumes. Times in OUT/times.txt.
  python logs/rerun/five/rerun.py [step or group names]     (groups: a = search A's steps, b = search B's, c = the critical-problem
                                                              solver's, r = the reviewers')
Environment: BIN (the builds: <name>.exe for each src/<name>.cpp, and r5/fivecheck_sieve.exe from src/r5/), OUT (default
BIN/out), FIX5 (Langevin and Hou's fix-5.txt), VENVPY (a python with numpy, for five_aut.py; default this one), OMP_NUM_THREADS.
The drivers and some programs read and write fixed paths of Rowen's working tree (runs\\i12\\..., runs\\check5\\...,
runs\\bin\\<build>.exe), so they run unchanged in OUT/work, a folder laid out like it: the inputs from data/ and logs/, the builds
copied under the names the drivers call. Claude Opus 5.5, for this folder, 2026-11-11."""
import concurrent.futures as cf
import filecmp
import glob
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
BUNDLE = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BIN = os.path.abspath(os.environ.get("BIN", HERE))
O = os.path.abspath(os.environ.get("OUT", os.path.join(BIN, "out")))
FIX5 = os.environ.get("FIX5", "fix-5.txt")
PY, VPY = sys.executable, os.environ.get("VENVPY", sys.executable)
W = os.path.join(O, "work")
os.makedirs(O, exist_ok=True)
os.chdir(BUNDLE)
o = lambda f: os.path.join(O, f)
b = lambda n: os.path.join(BIN, n + ".exe")
wb = lambda n: os.path.join(W, "runs", "bin", n + ".exe")


def cp(src, dst):
    # left alone when already there: a second invocation's setup mustn't touch a build the first is running
    if os.path.exists(dst) and filecmp.cmp(src, dst, shallow=False):
        return
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(src, dst)


def setup():
    for d in ("apn", "runs/bin", "runs/i12/reg", "runs/i11", "runs/check5/review"):
        os.makedirs(os.path.join(W, d), exist_ok=True)
    for d in ("r5", "r5head", "k4sel"):
        os.makedirs(o(d), exist_ok=True)
    for f in glob.glob("src/*.cmd") + glob.glob("src/*.py"):
        cp(f, os.path.join(W, "apn", os.path.basename(f)))
    for src, dst in (("data/p5.txt", "runs/i12/p5.txt"), ("data/aut5.txt", "runs/i12/aut.txt"), ("data/p6.txt", "runs/i11/p6.txt"),
                     ("data/configs5.txt", "runs/check5/configs.txt"), ("data/known4.txt", "runs/check5/known4.txt"),
                     ("logs/check5-review/forms.txt", "runs/check5/review/forms.txt")):
        cp(src, os.path.join(W, dst))
    # the names the drivers call: cover5v ran search A and cover5u the regression, both now src/cover.cpp; fivecheck_r5 is
    # src/r5/; the controls' fivecheck_ctl, the checks' _b builds and the k4list builds are now src/'s fivecheck_sieve.cpp,
    # fivecheck_classes.cpp and fivecheck_configs.cpp (Build says what changed since the builds that ran)
    for src, names in (("cover", ("cover5v", "cover5u")), ("five_orbits", ("five_orbits",)), ("five_critical", ("five_critical",)),
                       ("r5/fivecheck_sieve", ("fivecheck_r5",)), ("fivecheck_sieve", ("fivecheck_ctl", "fivecheck_sieve_b")),
                       ("fivecheck_classes", ("fivecheck_cl2", "fivecheck_cl3", "fivecheck_classes_b")),
                       ("fivecheck_configs", ("fivecheck_configs_b",))):
        for n in names:
            if not os.path.exists(wb(n)):   # a rebuild differs in its link timestamp, and a running exe can't be overwritten
                cp(b(src), wb(n))
    for src, dst in (("logs/five-a", "runs/i12/cover"), ("logs/check5-k4", "runs/check5/k4"),
                     ("logs/check5-review/nou", "runs/check5/review/nou")):
        shutil.copytree(src, os.path.join(W, dst), dirs_exist_ok=True, copy_function=cp)
    # the last review's definition check: the first 10 of its 494 items
    with open(o("defcheck10.txt"), "w") as f:
        f.writelines(open("logs/review-five-page-defcheck/items.txt").readlines()[:10])
    # the reference's sample: its first 300 sets
    with open(o("ref_sample300.txt"), "w") as f:
        f.writelines(open("logs/check5-review/ref_sample.txt").readlines()[:300])
    # three of the cheapest dimension-4 representatives (72 to 75 s each in the run), as lists for 'list'
    for c, idx in ((12, (3, 12)), (13, (2,))):
        with open(o(f"k4sel_c{c}.txt"), "w") as f:
            for i in idx:
                f.write(" ".join(open(f"logs/check5-k4/l_r5_c{c}_{i:06d}.txt").readline().split()[1:6]) + "\n")


def known_forms(line):
    return open("data/known4.txt").read().split("\n")[line - 1].split()[3:8]


def hyp_sample():
    """Cores of logs/five-hyp/ run again as five_hyp.py's run_one runs them (the whole problem, then the 255 hyperplanes), into
    OUT/hyp/: the eight whose hyperplanes aren't all SAT, and the first two of each configuration."""
    files = sorted(glob.glob("logs/five-hyp/p5-*.txt"))
    pick = [f for f in files if not re.search(r"^sat 255 ", open(f).read(), re.M)]
    for tag in sorted({os.path.basename(f).rsplit("_", 1)[0] for f in files}):
        pick += [f for f in files if os.path.basename(f).rsplit("_", 1)[0] == tag][:2]
    os.makedirs(o("hyp"), exist_ok=True)

    def one(f):
        forms = re.search(r"forms (.*)", open(f).readline()).group(1).split()
        full = subprocess.run([b("five_critical")] + forms + ["1", "0"], capture_output=True, text=True).stdout
        hyp = subprocess.run([b("five_critical")] + forms + ["1", "0", "-1", "1"], capture_output=True, text=True).stdout
        keep = [l for l in full.splitlines() if l.startswith(("G5", "iter", "SAT", "cap", "core rank", "no 3-dim"))]
        keep += [l for l in hyp.splitlines() if l.startswith(("sat ", "H ")) and ("UNSAT" in l or "cap" in l or l.startswith("sat "))]
        with open(o("hyp/" + os.path.basename(f)), "w") as g:
            g.write(open(f).readline() + "\n".join(keep) + "\n")

    with cf.ThreadPoolExecutor(int(os.environ.get("OMP_NUM_THREADS", "4"))) as ex:
        list(ex.map(one, sorted(set(pick))))
    print(f"{len(set(pick))} cores run again")
    return 0


A = [
    ("a-p3-enum5", [b("pspread3"), "enum", "5", o("p5.txt"), "--ckpt", o("p5.ckpt"), "--wcheck", "1000"], None),
    ("a-p3-holes5", [b("pspread3"), "holes", "5", "data/p5.txt", o("holes5.txt")], None),
    ("a-holes-tally", [PY, "src/holes_tally.py", o("holes5.txt"), "data/p5.txt", FIX5, "22", "0"], None),
    ("a-five-aut", [VPY, "src/five_aut.py", "data/p5.txt", o("aut.txt")], None),
    ("a-five-cores", ["cmd", "/c", os.path.join(W, "apn", "five_cores.cmd")], W),
    ("a-rv5-indep", [b("rv5_indep"), r"runs\i12\p5.txt", r"runs\i12\aut.txt", r"runs\i12\near", r"runs\i12\reps"], W),
    ("a-five-tally", [PY, r"apn\five_slices.py", "tally"], W),
    ("a-five-regress", ["cmd", "/c", os.path.join(W, "apn", "five_regress.cmd")], W),
    ("a-vs-b", [b("five_vs_b"), "data/p5.txt", o("vs-b.txt")], None),
    ("a-classmatch", [PY, "src/five_classmatch.py", o("vs-b.txt"), "five-check.md", "logs/check5-bundle/gamma5.txt"], None),
    ("a-cover-p11", [b("cover"), "data/p11.txt"], None),
]
B = [
    ("b-configs", [b("fivecheck_configs"), o("configs5.txt")], None),
    ("b-covers", [PY, r"apn\fivecheck_covers.py"], W),
    ("b-bundle", ["cmd", "/c", os.path.join(W, "apn", "fivecheck_bundle.cmd")], W),
    ("b-bundle2", ["cmd", "/c", os.path.join(W, "apn", "fivecheck_bundle2.cmd")], W),
    ("b-k4list", ["cmd", "/c", os.path.join(W, "apn", "fivecheck_k4list.cmd")], W),
    ("b-k4all", ["cmd", "/c", os.path.join(W, "apn", "fivecheck_k4all.cmd")], W),
    ("b-k4orbits", [PY, r"apn\fivecheck_k4orbits.py"], W),
    ("b-k4reps-c12", [wb("fivecheck_ctl"), "list", "5", "12", o("k4sel_c12.txt"), o("k4sel"), "pivot"], W),
    ("b-k4reps-c13", [wb("fivecheck_ctl"), "list", "5", "13", o("k4sel_c13.txt"), o("k4sel"), "pivot"], W),
    ("b-r5-c0q0", [b("r5/fivecheck_sieve"), "run", "5", "0", "0", "1", o("r5"), "onlyk=5"], W),
    ("b-r5head-c0q0", [b("fivecheck_sieve"), "run", "5", "0", "0", "1", o("r5head"), "onlyk=5"], W),
]
C = [(f"c-pos-{n}", [b("five_critical")] + known_forms(n) + ["1", "0"], None) for n in (1, 7, 20)]
C += [("c-hyp-sample", hyp_sample, None), ("c-hyp-tally", [PY, "src/five_hyp.py", "logs/five-hyp", "--tally"], None)]
R = [
    ("r-rvb5-count", [b("rvb5_ref"), "count"], W),
    ("r-rvb5-known", [b("rvb5_ref"), "known", "runs/check5/known4.txt"], W),
    ("r-rvb5-sample300", [b("rvb5_ref"), "list", o("ref_sample300.txt")], W),
    ("r-rvb5-gamma", [b("rvb5_gamma"), "gamma"], W),
    ("r-closure-c21", [b("rvb5_gamma5"), "21", "0", r"runs\check5\review\nou\r5_c21_q00.txt"], W),
    ("r-closure-c20", [b("rvb5_gamma5"), "20", "3", r"runs\check5\review\nou\r5_c20_q03.txt"], W),
    ("r-nou-check", [PY, r"apn\nou_check.py"], W),
    ("r-defcheck10", [b("defcheck"), o("defcheck10.txt")], None),
    ("r-defcheck-k4", [b("defcheck2"), "logs/review-five-page-defcheck/k4.txt", "5"], None),
]
PRE = {"a-five-cores": lambda: [os.remove(p) for p in glob.glob(os.path.join(W, "runs", "i12", "reps", "summary.txt"))]}
groups = {"a": [s[0] for s in A], "b": [s[0] for s in B], "c": [s[0] for s in C], "r": [s[0] for s in R]}

setup()
only = set()
for a in sys.argv[1:]:
    only |= set(groups.get(a, [a]))
for name, cmd, cwd in A + B + C + R:
    if only and name not in only:
        continue
    log = o(name + ".log")
    if os.path.exists(log):
        tail = open(log, errors="replace").read().rstrip().splitlines()[-1:]
        if tail and tail[0].startswith("# exit"):
            continue
    PRE.get(name, lambda: None)()
    t0 = time.time()
    with open(log, "w") as lf:
        if callable(cmd):
            sys.stdout, keep = lf, sys.stdout
            try:
                rc = cmd()
            finally:
                sys.stdout = keep
        else:
            rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, cwd=cwd).returncode
    with open(log, "a") as lf:
        lf.write(f"# exit {rc}\n")
    with open(o("times.txt"), "a") as tf:
        tf.write(f"{name}: exit {rc}, {time.time() - t0:.0f} s\n")
    print(f"{name}: exit {rc}, {time.time() - t0:.0f} s", flush=True)
