"""Compare five's bundle re-runs (rerun.py's OUT folder) with the original runs. Ran from Rowen's working tree (the gleaner
root): the originals are this folder's logs/ and data/, and for what isn't here, search A's runs/i12/ (orbit lists) and search B's
runs/check5/ (the dimension-4 tuple lists, the sieve-only outputs of k4all).
  python findings/code/apn8-sixteen/logs/rerun/five/compare.py OUT > compare.log
Claude Opus 5.5, for this folder, 2026-11-11."""
import glob
import os
import re
import sys

RR = sys.argv[1]
W = os.path.join(RR, "work")
B = "findings/code/apn8-sixteen/"
L = B + "logs/"


def lines(p):
    return open(p, errors="replace").read().splitlines()


def body(name):
    return [l for l in lines(os.path.join(RR, name)) if not l.startswith("# exit")]


def have(name):
    p = os.path.join(RR, name + ".log")
    ok = os.path.exists(p) and lines(p)[-1:] and lines(p)[-1].startswith("# exit")
    if not ok:
        print(f"{name}: not run")
    return ok


def exit_of(name):
    return lines(os.path.join(RR, name + ".log"))[-1]


def notime(ls):
    # timings ("31 s", "1.2 s", "(0 s)", "seconds 3.1", "0.5 s per item", ", 535 s" at a done line's end) and paths
    out = []
    for l in ls:
        l = re.sub(r"\b\d+(\.\d+)? s\b|\(\d+(\.\d+)? s\)|seconds [\d.]+|[\d.]+ s per item per thread", "T", l)
        l = re.sub(r"[A-Za-z]:[\\/]\S*[\\/]", "", l)
        out.append(l.rstrip())
    return out


def same(label, a, b, extra=""):
    print(f"{label}: {'IDENTICAL' if a == b else 'DIFFER'}{(' ' + extra) if extra else ''}")
    return a == b


def filebytes(p):
    return open(p, "rb").read().replace(b"\r\n", b"\n")


print("== search A")
if have("a-p3-enum5"):
    same("pspread3 enum 5: p5.txt vs data/p5.txt", filebytes(os.path.join(RR, "p5.txt")), filebytes(B + "data/p5.txt"))
    keep = lambda ls: [re.sub(r", w \d+\)", ")", l) for l in notime(ls) if re.match(r"\s*w=|orbits \d|w spot|GL\(4,2\)", l)]
    same("pspread3 enum 5: log vs pspread5-output.txt (w lines, orbits, spot checks)", keep(body("a-p3-enum5.log")), keep(lines(L + "pspread5-output.txt")), exit_of("a-p3-enum5"))
if have("a-p3-holes5"):
    same("pspread3 holes 5: holes5.txt vs logs/holes5.txt, as sets of lines", set(lines(os.path.join(RR, "holes5.txt"))), set(lines(L + "holes5.txt")))
    same("pspread3 holes 5: log vs holes5-output.txt (timings aside)", notime(body("a-p3-holes5.log")), notime(lines(L + "holes5-output.txt")))
if have("a-holes-tally"):
    same("holes_tally.py: log vs holes5-tally.txt", body("a-holes-tally.log"), lines(L + "holes5-tally.txt"))
if have("a-five-aut"):
    same("five_aut.py: aut.txt vs data/aut5.txt", filebytes(os.path.join(RR, "aut.txt")), filebytes(B + "data/aut5.txt"))
    same("five_aut.py: log vs five-aut-output.txt", body("a-five-aut.log"), lines(L + "five-aut-output.txt"))
if have("a-five-cores"):
    same("five_cores.cmd: summary vs five-orbits-summary.txt", lines(os.path.join(W, "runs/i12/reps/summary.txt")), lines(L + "five-orbits-summary.txt"))
    tags = [l.split()[0] for l in lines(B + "data/p5.txt") if l.strip()]
    bad = [t for t in tags if filebytes(os.path.join(W, f"runs/i12/reps/{t}.txt")) != filebytes(f"runs/i12/reps/{t}.txt")]
    print(f"five_cores.cmd: the 22 orbit lists vs the run's (runs/i12/reps/): {22 - len(bad)} identical byte for byte, differing {bad}")
if have("a-rv5-indep"):
    same("rv5_indep: log vs review-five-orbits-indep.txt, as sets of lines", set(body("a-rv5-indep.log")), set(lines(L + "review-five-orbits-indep.txt")),
         exit_of("a-rv5-indep") + "; " + [l for l in body("a-rv5-indep.log") if l.startswith("total")][0])
if have("a-five-tally"):
    same("five_slices.py tally over logs/five-a/: vs five-tally-output.txt", body("a-five-tally.log"), lines(L + "five-tally-output.txt"), exit_of("a-five-tally"))
if have("a-five-regress"):
    strip = lambda ls: [re.sub(r"threads \d+", "threads N", l) for l in notime(ls) if not re.match(r"frontier \d+/\d+ ", l)]
    for k in range(3):
        same(f"five_regress.cmd: slice s{k} vs five-regress-s{k}.txt (timings and the thread count aside)", strip(lines(os.path.join(W, f"runs/i12/reg/s{k}.out"))), strip(lines(L + f"five-regress-s{k}.txt")))
if have("a-vs-b"):
    same("five_vs_b: vs-b.txt vs logs/five-vs-b.txt", filebytes(os.path.join(RR, "vs-b.txt")), filebytes(L + "five-vs-b.txt"))
if have("a-classmatch"):
    same("five_classmatch.py: log vs check5-bundle/classmatch.txt", body("a-classmatch.log"), lines(L + "check5-bundle/classmatch.txt"), exit_of("a-classmatch"))
if have("a-cover-p11"):
    a, o = body("a-cover-p11.log"), lines(L + "eleven-output.txt")
    pick = lambda ls, p: [l for l in ls if l.startswith(p)]
    done = lambda ls: [re.sub(r", \d+ s$", "", l) for l in pick(ls, "done:")]
    print(f"cover data/p11.txt (src/cover.cpp, no five flags) vs eleven-output.txt (cover3): nodes line equal "
          f"{pick(a, 'nodes by span') == pick(o, 'nodes by span')}; done line (sans time) equal {done(a) == done(o)}; NEAR "
          f"{len(pick(a, 'NEAR'))} vs {len(pick(o, 'NEAR'))}, equal as sets {set(pick(a, 'NEAR')) == set(pick(o, 'NEAR'))}; {done(a)}")

print("== search B")
if have("b-configs"):
    same("fivecheck_configs: configs5.txt vs data/configs5.txt", filebytes(os.path.join(RR, "configs5.txt")), filebytes(B + "data/configs5.txt"))
    same("fivecheck_configs: log vs check5-bundle/configs.txt (paths aside)", notime(body("b-configs.log")), notime(lines(L + "check5-bundle/configs.txt")))
if have("b-covers"):
    same("fivecheck_covers.py: log vs check5-bundle/covers.txt", body("b-covers.log"), lines(L + "check5-bundle/covers.txt"))
for step, names in (("b-bundle", ["gamma4", "slevel5", "slevel4", "gamma5"] + [f"flat5_c{c}" for c in (0, 2, 4, 8, 10, 13, 14, 19, 20, 21)] + [f"flat4_c{c}" for c in (0, 1, 2)]),
                    ("b-bundle2", ["configs", "configs-rerun", "configs-fc", "covers", "blockview", "plant", "plant-pivot"])):
    if have(step):
        diff = [n for n in names if notime(lines(os.path.join(W, f"runs/check5/bundle/{n}.txt"))) != notime(lines(L + f"check5-bundle/{n}.txt"))]
        print(f"{step}: {len(names) - len(diff)} of {len(names)} outputs identical to logs/check5-bundle/ (timings and paths aside); differing {diff}; {exit_of(step)}")
if have("b-k4list"):
    orig = [l for l in lines(L + "20261106-031934-check5-k4list.log") if l.startswith("cfg ")]
    same("k4list: per-class lines vs the k4list job's", [l for l in body("b-k4list.log") if l.startswith("cfg ")], orig)
    bad = [p for p in sorted(glob.glob("runs/check5/k4_c*.txt")) if filebytes(p) != filebytes(os.path.join(W, p))]
    print(f"k4list: the 22 lists k4_c*.txt vs the run's (runs/check5/): differing {bad}")
if have("b-k4all"):
    orig = sorted(glob.glob("runs/check5/k4all/*.txt"))
    st = lambda p: [l for l in lines(p) if not l.startswith("seconds")]
    bad = [p for p in orig if st(p) != st(os.path.join(W, p))]
    print(f"k4all: {len(orig) - len(bad)} of {len(orig)} sieve-only outputs identical to the run's (runs/check5/k4all/; the seconds line aside), "
          f"{len(glob.glob(os.path.join(W, 'runs/check5/k4all/*.txt')))} written; differing {bad[:5]}; {exit_of('b-k4all')}")
    for pat in ("k4all_c*.txt", "k4orbit_c*.txt"):
        bad = [p for p in sorted(glob.glob("runs/check5/" + pat)) if filebytes(p) != filebytes(os.path.join(W, p))]
        print(f"k4all: lists {pat} vs the run's: differing {bad}")
if have("b-k4orbits"):
    print("fivecheck_k4orbits.py:", " | ".join(body("b-k4orbits.log")[-3:]))
for c, pairs in ((12, ((0, 3), (1, 12))), (13, ((0, 2),))):
    if have(f"b-k4reps-c{c}"):
        for new, old in pairs:
            st = lambda p: [l for l in lines(p) if not l.startswith("seconds")]
            same(f"k4 representative class {c} list line {old} (list mode, src/fivecheck_sieve.cpp) vs logs/check5-k4/ (the seconds line aside)",
                 st(os.path.join(RR, f"k4sel/l_r5_c{c}_{new:06d}.txt")), st(L + f"check5-k4/l_r5_c{c}_{old:06d}.txt"))
if have("b-r5-c0q0"):
    st = lambda p: [l for l in lines(p) if not l.startswith("seconds")]
    same("main-run slice c0 q0 on src/r5/ vs logs/check5-r5/r5_c00_q00.txt (the seconds line aside)", st(os.path.join(RR, "r5/r5_c00_q00.txt")), st(L + "check5-r5/r5_c00_q00.txt"))
if have("b-r5head-c0q0"):
    print("the same slice on src/fivecheck_sieve.cpp:", lines(os.path.join(RR, "r5head/r5_c00_q00.txt"))[1])
    print("  the run's                              :", lines(L + "check5-r5/r5_c00_q00.txt")[1])

print("== the critical-problem solver")
for n in (1, 7, 20):
    if have(f"c-pos-{n}"):
        print(f"positive control, data/known4.txt line {n}, its first five forms:", " | ".join(l for l in body(f"c-pos-{n}.log") if l.startswith(("G5", "SAT", "cap", "no 3"))), exit_of(f"c-pos-{n}"))
if have("c-hyp-sample"):
    got = sorted(glob.glob(os.path.join(RR, "hyp", "p5-*.txt")))
    bad = [os.path.basename(p) for p in got if lines(p) != lines(L + "five-hyp/" + os.path.basename(p))]
    print(f"five_critical on {len(got)} of the 2,200 cores again (the 8 with UNSAT hyperplanes and the first 2 of each configuration): "
          f"{len(got) - len(bad)} files identical to logs/five-hyp/, differing {bad}")
if have("c-hyp-tally"):
    norm = lambda ls: [re.sub(r"\S*[\\/](p5-\S+\.txt)", r"\1", l) for l in ls]
    same("five_hyp.py --tally over logs/five-hyp/ vs five-hyp-tally.txt (paths aside)", norm(body("c-hyp-tally.log")), norm(lines(L + "five-hyp-tally.txt")))

print("== the reviewers' code")
if have("r-rvb5-count"):
    orig = [l for l in lines(L + "20261106-172106-rvb5-count.log") if not l.startswith("[runner") and l.strip()]
    same("rvb5_ref count (src/): output vs the rvb5-count job's", body("r-rvb5-count.log"), orig, exit_of("r-rvb5-count"))
    for f in ("ref_count.txt", "ref_sample.txt", "ref_k4.txt"):   # the last two are written in the order the threads finish
        same(f"rvb5_ref count: {f} vs logs/check5-review/{f}, as sets of lines", set(lines(os.path.join(W, "runs/check5/review", f))),
             set(lines(L + "check5-review/" + f)))
if have("r-rvb5-known"):
    orig = [l for l in lines(L + "20261106-172106-rvb5-known.log") if not l.startswith("[runner") and l.strip()]
    same("rvb5_ref known (src/): vs the rvb5-known job's", body("r-rvb5-known.log"), orig)
if have("r-rvb5-sample300"):
    a, o = set(body("r-rvb5-sample300.log")) - {"REFLISTDONE"}, set(lines(L + "check5-review/ref_sample_out.txt"))
    print(f"rvb5_ref list (src/) on the sample's first 300 sets: {len(a)} lines, {len(a & o)} of them in ref_sample_out.txt, {len(a - o)} not; {exit_of('r-rvb5-sample300')}")
if have("r-rvb5-gamma"):
    print("rvb5_gamma gamma:", body("r-rvb5-gamma.log")[-1], "| lines 'same':", sum(" same" in l for l in body("r-rvb5-gamma.log")))
for c in (21, 20):
    if have(f"r-closure-c{c}"):
        print(f"rvb5_gamma5 closure, class {c}:", " ".join(body(f"r-closure-c{c}.log")))
if have("r-nou-check"):
    same("nou_check.py: output vs logs/check5-review/nou_check_out.txt", body("r-nou-check.log"), lines(L + "check5-review/nou_check_out.txt"))
if have("r-defcheck10"):
    a = [l for l in body("r-defcheck10.log") if l.startswith("ITEM")]
    o = set(lines(L + "review-five-page-defcheck/items.out"))
    print(f"defcheck (src/) on the first 10 of the last review's 494 items: {len(a)} ITEM lines, {sum(l in o for l in a)} of them in its items.out; {exit_of('r-defcheck10')}")
if have("r-defcheck-k4"):
    print("defcheck2 (src/) on the four dimension-4 representatives, stopping at dimension 5:", " | ".join(body("r-defcheck-k4.log")))
