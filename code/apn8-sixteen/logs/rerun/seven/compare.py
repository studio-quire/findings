"""Compare seven's bundle re-runs (rerun.py's OUT folder) with the original runs. Ran from Rowen's working tree: the
originals are search A's runs/i10/ and search B's runs/check7/ (their outputs are in logs/ here under other names).
  python findings/code/apn8-sixteen/logs/rerun/seven/compare.py OUT > compare.log
Claude Opus 5.5, for this folder, 2026-11-01."""
import glob
import os
import re
import sys

RR = sys.argv[1]
B = "findings/code/apn8-sixteen/"
A = "runs/i10/"
C7 = "runs/check7/"


def lines(p):
    return open(p, errors="replace").read().splitlines()


def body(p):
    return [l for l in lines(p) if not l.startswith("# exit")]


def notime(ls):
    # timings ("31 s", "[4.8 s, total 5 s]", "(0 s)", "secs 31", "thread-s canon 13313") and the re-run's paths
    out = []
    for l in ls:
        l = re.sub(r"\[[\d.]+ s(, total \d+ s)?\]|\(\d+(\.\d+)? s\)|\b\d+(\.\d+)? s\b|secs \d+|thread-s canon \d+", "T", l)
        l = re.sub(r"[A-Za-z]:\\\S*\\", "", l)
        out.append(l)
    return out


def report(label, ok, extra=""):
    print(f"{label}: {'IDENTICAL' if ok else 'DIFFER'}{(' ' + extra) if extra else ''}")


def rlines(files):
    out, dup = {}, 0
    for f in files:
        for l in lines(f):
            if l.startswith("R "):
                v = l.split()
                dup += int(v[1]) in out
                out[int(v[1])] = v[2:]
    return out, dup


def cores(files):
    by = {}
    rec = re.compile(r"K (\d+) (\d+) ((?:[0-9a-f]{7} )+)")
    for f in files:
        for m in rec.finditer(open(f).read()):
            if len(m[3].split()) == int(m[2]):
                by.setdefault(int(m[1]), set()).add((int(m[2]), m[3].strip()))
    return by


def have(name):
    p = os.path.join(RR, name)
    ok = os.path.exists(p) and lines(p)[-1:] and lines(p)[-1].startswith("# exit")
    if not ok:
        print(f"{name}: not run (yet)")
    return ok


print("== search A")
if have("p3-enum7.log"):
    report("pspread3 enum 7: p7.txt vs data/p7.txt", open(os.path.join(RR, "p7.txt"), "rb").read() == open(B + "data/p7.txt", "rb").read())
    # the spot checks are every N-th clique of each thread, so their number depends on the thread count; "w <s>" on the
    # class lines is time spent on w
    keep = lambda ls: [re.sub(r"spot checks: \d+", "spot checks: N", re.sub(r", w \d+\)", ")", l)) for l in notime(ls)
                       if re.match(r"\s*w=|orbits \d|w spot|wrote|reduced cliques|GL\(4,2\)|class \d done", l)]
    a, b = keep(body(os.path.join(RR, "p3-enum7.log"))), keep(lines(A + "p7.out"))
    spot = [l for l in body(os.path.join(RR, "p3-enum7.log")) + lines(A + "p7.out") if l.startswith("w spot")]
    report("pspread3 enum 7: log vs p7.out (class totals, w lines, orbits; timings and the spot checks' number aside)",
           a == b, f"({len(a)} lines; {spot})")
if have("p3-selftest7.log"):
    print("pspread3 selftest 7:", [l for l in body(os.path.join(RR, "p3-selftest7.log")) if "selftest" in l])
if have("p3-holes7.log"):
    report("pspread3 holes 7: holes7.txt vs logs/holes7.txt as sets of lines",
           sorted(lines(os.path.join(RR, "holes7.txt"))) == sorted(lines(B + "logs/holes7.txt")))
    a, b = notime(body(os.path.join(RR, "p3-holes7.log"))), notime(lines(A + "holes/holes7.out"))
    report("pspread3 holes 7: log vs holes7.out (timings and paths aside)", a == b, f"({len(a)} vs {len(b)} lines)")
for name, split in (("cover-c2.log", False), ("cover-c2-split5.log", True)):
    if not have(name):
        continue
    a, b = body(os.path.join(RR, name)), lines(A + "cover/c2.out")
    pick = lambda ls, p: [l for l in ls if l.startswith(p)]
    done = lambda ls: [re.sub(r", \d+ s$", "", l) for l in pick(ls, "done:")]
    na, nb = sorted(pick(a, "NEAR")), sorted(pick(b, "NEAR"))
    print(f"{name} vs cover3's c2.out: nodes line equal {pick(a, 'nodes by span') == pick(b, 'nodes by span')}; done line "
          f"(sans time) equal {done(a) == done(b)}; near-miss line equal {pick(a, 'near misses') == pick(b, 'near misses')}; "
          f"NEAR {len(na)} vs {len(nb)}, equal as sets {na == nb}" + (f"; {pick(a, 'phase 1 done')}" if split else ""))

print("== search B")
if have("s7-count4.log"):
    report("count4 vs count4.out (timings aside)", notime(body(os.path.join(RR, "s7-count4.log"))) == notime(lines(C7 + "count4.out")))
parts = [p for p in sorted(glob.glob(os.path.join(RR, "parta-c*-s*.txt")))
         if (lines(os.path.join(RR, "s7-" + os.path.basename(p).replace(".txt", ".log")))[-1:] or [""])[0].startswith("# exit")]
if parts:
    nbad = nlog = 0
    for p in parts:
        n = os.path.basename(p)
        a = [re.sub(r" secs \d+$", "", l) for l in lines(p)]
        b = [re.sub(r" secs \d+$", "", l) for l in lines(C7 + n)]
        nbad += a != b
        lg = os.path.join(RR, "s7-" + n.replace(".txt", ".log"))
        if os.path.exists(lg):
            nlog += notime(body(lg)) != notime(lines(C7 + n.replace(".txt", ".out")))
        mx = os.path.join(RR, n.replace("parta", "maxc"))
        nbad += os.path.getsize(mx) != 0
    print(f"parta: {len(parts)} slices re-made; slice files differing from the originals (the CLASS line's secs aside) or "
          f"with a non-empty maxc file: {nbad}; outputs differing (timings aside): {nlog}")
if have("s7-merge.log"):
    report("merge: orbits.txt vs data/orbits7.txt", lines(os.path.join(RR, "orbits.txt")) == lines(B + "data/orbits7.txt"))
    report("merge: output vs merge.out (timings aside)", notime(body(os.path.join(RR, "s7-merge.log"))) == notime(lines(C7 + "merge.out")))
if have("s7-configs.log"):
    report("configs: configs7.txt vs data/configs7.txt", lines(os.path.join(RR, "configs7.txt")) == lines(B + "data/configs7.txt"))
    a, b = notime(body(os.path.join(RR, "s7-configs.log"))), notime(lines(C7 + "configs.out"))
    d = [(x, y) for x, y in zip(a, b) if x != y]
    report("configs: output vs configs.out (timings and paths aside)", a == b,
           f"({len(a)} vs {len(b)} lines; {len(d)} differ: {d[:3]})")
if have("s7-mapa.log"):
    report("mapa: output vs mapa.out", body(os.path.join(RR, "s7-mapa.log")) == lines(C7 + "mapa.out"))
for run, orig in (("asc", "asc-r*.prog"), ("desc", "desc-r*.prog"), ("sweep", "sweep-r*.prog"), ("phase1", "phase1-r*.prog"),
                  ("s8-b4728", "s8-b4728.prog")):
    pf = os.path.join(RR, run + ".prog")
    if not os.path.exists(pf):
        continue
    a, dup = rlines([pf])
    b, _ = rlines(glob.glob(C7 + orig))
    bad = [i for i in a if a[i] != b.get(i)]
    rng = f"{min(a)}..{max(a)}" if a else "-"
    print(f"{run}: re-run R lines {len(a)} (configurations {rng}, {dup} repeated), unequal to the original's {len(bad)} {bad[:8]}")
    for i in bad[:3]:
        names = ([f"nodes[{d}]" for d in range(9)] + [f"tries[{d}]" for d in range(9)] + [f"p1[{d}]" for d in range(9)]
                 + ["bentext", "full"] + [f"triesD[{d}]" for d in range(9)] + [f"cores[{d}]" for d in range(9)] + ["fp"]
                 + [f"swFound[{d}]" for d in range(9)] + ["swChecks", "swSols"])
        d = [names[k] for k, (x, y) in enumerate(zip(a[i], b[i])) if x != y]
        print(f"   config {i}: fields differing {d}" + (" (the sweep's cache, for its dimension-5 cores: by design)"
                                                         if run == "sweep" and d == ["swChecks"] else ""))
    cf = os.path.join(RR, run + ".cores")
    if os.path.exists(cf):
        ca, cb = cores([cf]), cores(glob.glob(C7 + orig.replace(".prog", ".cores")))
        diff = [i for i in a if ca.get(i, set()) != cb.get(i, set())]
        print(f"   core files: {sum(len(ca.get(i, ())) for i in a)} cores on those configurations, sets differing from the original's on {len(diff)} {diff[:8]}")
if have("s7-ctl-ind-bl22.log"):
    report("ctl ind-bl22: ctl-ind-bl22.txt vs the original", lines(os.path.join(RR, "ctl-ind-bl22.txt")) == lines(C7 + "ctl-ind-bl22.txt"))
    a, b = notime(body(os.path.join(RR, "s7-ctl-ind-bl22.log"))), notime(lines(C7 + "ctl-ind-bl22.out"))
    report("ctl ind-bl22: output vs ctl-ind-bl22.out (timings and paths aside)", a == b,
           f"({len(a)} vs {len(b)} lines; only in the re-run {[x for x in a if x not in b]}, only in the original {[x for x in b if x not in a]})")
if have("s7-ctlrot.log"):
    a = [l for l in body(os.path.join(RR, "s7-ctlrot.log")) if l.startswith("rep ")]
    keep = {l.split()[1] for l in a}
    b = [l for l in lines(C7 + "ctlrot-0.out") if l.startswith("rep ") and l.split()[1] in keep]
    report("ctlrot: rep lines vs ctlrot-0.out", sorted(a) == sorted(b), str(a))
strip = lambda x: re.sub(r" \[\d+ s\]", "", x)
if have("s7-numpy-3.log"):
    a = [strip(l) for l in body(os.path.join(RR, "s7-numpy-3.log")) if l.startswith("config 3:")]
    report("numpy config 3 vs numpy-3.out (time aside)", a == [strip(l) for l in lines(C7 + "numpy-3.out")], str(a))
if have("s7-numpy-ctl.log"):
    a = [strip(l) for l in body(os.path.join(RR, "s7-numpy-ctl.log")) if l.startswith("control line")]
    b = [strip(l) for l in lines(C7 + "numpy-ctl-bl22.out") if l.startswith("control line")][:len(a)]
    report("numpy control: first control vs numpy-ctl-bl22.out's first (time aside)", a == b and len(a) == 1, str(a))
if have("s7-selftest.log"):
    report("selftest vs selftest.out", notime(body(os.path.join(RR, "s7-selftest.log"))) == notime(lines(C7 + "selftest.out")))
