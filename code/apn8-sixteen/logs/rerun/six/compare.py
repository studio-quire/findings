"""Compare six's bundle re-runs (rerun.py's OUT folder) with the original runs. Ran from Rowen's working tree: the
originals are search A's runs/i11/ and runs/cover5/, and search B's runs/check6/ (their outputs are in logs/ here under
other names, without the NEAR lines and cores, which are too large; logs/six-near.txt has those per configuration).
  python findings/code/apn8-sixteen/logs/rerun/six/compare.py OUT > compare.log
Claude Opus 5.5, for this folder, 2026-11-04."""
import glob
import os
import re
import sys

RR = sys.argv[1]
B = "findings/code/apn8-sixteen/"
A = "runs/i11/"
A5 = "runs/cover5/"
C6 = "runs/check6/"


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


def have(name):
    p = os.path.join(RR, name)
    ok = os.path.exists(p) and lines(p)[-1:] and lines(p)[-1].startswith("# exit")
    if not ok:
        print(f"{name}: not run (yet)")
    return ok


def pick(ls, p):
    return [l for l in ls if l.startswith(p)]


def cover_vs(name, orig, label, hole=True):
    # a cover run against an original: its per-dimension nodes, done line (time aside), near-miss count, the hole
    # stage's filter counts when both have them, and the NEAR lines as sets
    if not have(name):
        return
    a, b = body(os.path.join(RR, name)), lines(orig)
    done = lambda ls: [re.sub(r", \d+ s$", "", l) for l in pick(ls, "done:")]
    na, nb = set(pick(a, "NEAR")), set(pick(b, "NEAR"))
    s = (f"{name} vs {label}: nodes line equal {pick(a, 'nodes by span') == pick(b, 'nodes by span')}; done line (sans time) "
         f"equal {done(a) == done(b)}; near-miss line equal {pick(a, 'near misses') == pick(b, 'near misses')}; "
         f"NEAR {len(pick(a, 'NEAR'))} vs {len(pick(b, 'NEAR'))}, equal as sets {na == nb}")
    if hole:
        s += f"; filter line equal {pick(a, 'hole-stage') == pick(b, 'hole-stage')}"
    print(s + f"; {pick(a, 'done:')[-1:]}")


print("== search A")
if have("p3-enum6.log"):
    report("pspread3 enum 6: p6.txt vs data/p6.txt", open(os.path.join(RR, "p6.txt"), "rb").read() == open(B + "data/p6.txt", "rb").read())
    # the spot checks are every N-th clique of each thread, so their number depends on the thread count; "w <s>" on the
    # class lines is time spent on w
    keep = lambda ls: [re.sub(r"spot checks: \d+", "spot checks: N", re.sub(r", w \d+\)", ")", l)) for l in notime(ls)
                       if re.match(r"\s*w=|orbits \d|w spot|wrote|reduced cliques|GL\(4,2\)|class \d done", l)]
    a, b = keep(body(os.path.join(RR, "p3-enum6.log"))), keep(lines(A + "p6.out"))
    spot = [l for l in body(os.path.join(RR, "p3-enum6.log")) + lines(A + "p6.out") if l.startswith("w spot")]
    report("pspread3 enum 6: log vs p6.out (class totals, w lines, orbits; timings and the spot checks' number aside)",
           a == b, f"({len(a)} lines; {spot})")
if have("p3-selftest6.log"):
    print("pspread3 selftest 6:", [l for l in body(os.path.join(RR, "p3-selftest6.log")) if "selftest" in l])
if have("p3-holes6.log"):
    report("pspread3 holes 6: holes6.txt vs logs/holes6.txt as sets of lines",
           sorted(lines(os.path.join(RR, "holes6.txt"))) == sorted(lines(B + "logs/holes6.txt")))
    a, b = notime(body(os.path.join(RR, "p3-holes6.log"))), notime(lines(A + "holes/h6.out"))
    report("pspread3 holes 6: log vs h6.out (timings and paths aside)", a == b, f"({len(a)} vs {len(b)} lines)")
for k in range(3):
    cover_vs(f"cover-c{k}.log", A + f"cover5/c{k}.out", f"cover5t's c{k}.out")
for k in range(2):
    cover_vs(f"cover-c{k}.log", A + f"cover/c{k}.out", f"cover4's c{k}.out (the old hole stage)", hole=False)
if have("cover-c2.log"):
    a, b = body(os.path.join(RR, "cover-c2.log")), lines(A + "cover/c2.out")
    na, nb = set(pick(a, "NEAR")), set(pick(b, "NEAR"))
    print(f"cover-c2.log vs cover4's cancelled c2.out: NEAR {len(pick(a, 'NEAR'))} vs {len(pick(b, 'NEAR'))}, equal as sets {na == nb}; "
          f"the cancelled run's last line {b[-1:]}")
cover_vs("cover-p11.log", B + "logs/eleven-output.txt", "eleven's output (cover3)", hole=False)
cover_vs("cover-p10.log", B + "logs/ten-output.txt", "ten's output (cover3)", hole=False)
cover_vs("cover-p11.log", A5 + "reg6t.out", "the regression reg6t.out (cover5t)")
cover_vs("cover-p10.log", A5 + "reg7t.out", "the regression reg7t.out (cover5t)")
cover_vs("cover-prof-tree40.log", A5 + "prof-tree40.out", "prof-tree40.out (cover5t)")
cover_vs("cover-prof-gray40.log", A5 + "prof-new40.out", "prof-new40.out (the Gray build cover5)")
cover_vs("cover-prof-tree40.log", A5 + "prof-old40.out", "prof-old40.out (cover5 --old)", hole=False)
for k in (4, 3, 2):
    if have(f"cover-ctl-k{k}.log"):
        a, b = body(os.path.join(RR, f"cover-ctl-k{k}.log")), lines(B + f"logs/rerun/cover-ctl-k{k}.log")
        done = lambda ls: [re.sub(r", \d+ s$", "", l) for l in pick(ls, "done:")]
        print(f"cover-ctl-k{k}: SOLUTION lines equal the old hole stage's {pick(a, 'SOLUTION') == pick(b, 'SOLUTION')}; done line (sans time) "
              f"equal {done(a) == done(b)}; {pick(a, 'done:')}")
if have("cover-ctl-k1.log"):
    a = body(os.path.join(RR, "cover-ctl-k1.log"))
    print("cover-ctl-k1:", pick(a, "SOLUTION"), pick(a, "done:"), pick(a, "nodes by"), pick(a, "near misses"), pick(a, "hole-stage"))

print("== search B")
for n, orig in (("s6-count3", "count3.out"), ("s6-selftest", "selftest.out"), ("s6-linselftest", "linselftest.out"), ("s6-merge", "merge.out")):
    if have(n + ".log"):
        a, b = notime(body(os.path.join(RR, n + ".log"))), notime(lines(C6 + orig))
        d = [(x, y) for x, y in zip(a, b) if x != y]
        report(f"{n} vs {orig} (timings aside)", a == b, f"({len(a)} vs {len(b)} lines; {len(d)} differ: {d[:3]})")
parts = [p for p in sorted(glob.glob(os.path.join(RR, "parta-c*-s0.txt"))) if have("s6-" + os.path.basename(p).replace(".txt", ".log"))]
if parts:
    nbad = nlog = 0
    for p in parts:
        n = os.path.basename(p)
        a = [re.sub(r" secs \d+$", "", l) for l in lines(p)]
        b = [re.sub(r" secs \d+$", "", l) for l in lines(C6 + n)]
        nbad += a != b
        lg = os.path.join(RR, "s6-" + n.replace(".txt", ".log"))
        # (progress lines dropped: the running clique counts in them depend on the threads' timing)
        keep = lambda ls: [l for l in notime(ls) if " first vertices, " not in l]
        nlog += keep(body(lg)) != keep(lines(C6 + n.replace(".txt", ".out")))
        nbad += os.path.getsize(os.path.join(RR, n.replace("parta", "maxc"))) != 0
    print(f"parta: {len(parts)} slices re-made; slice files differing from the originals (the CLASS line's secs aside) or "
          f"with a non-empty maxc file: {nbad}; outputs differing (timings aside): {nlog}")
if have("s6-merge.log"):
    report("merge: orbits.txt vs data/orbits6.txt", lines(os.path.join(RR, "orbits.txt")) == lines(B + "data/orbits6.txt"))
if have("s6-configs.log"):
    report("configs: configs6.txt vs data/configs6.txt", lines(os.path.join(RR, "configs6.txt")) == lines(B + "data/configs6.txt"))
    a, b = notime(body(os.path.join(RR, "s6-configs.log"))), notime(lines(C6 + "configs.out"))
    a = [l.replace("to configs6.txt", "to configs.txt") for l in a]   # the file it writes has another name here
    d = [(x, y) for x, y in zip(a, b) if x != y]
    report("configs: output vs configs.out (timings, paths and the output's name aside)", a == b, f"({len(a)} vs {len(b)} lines; {len(d)} differ: {d[:3]})")
if have("s6-mapa.log"):
    report("mapa: output vs mapa.out", body(os.path.join(RR, "s6-mapa.log")) == lines(C6 + "mapa.out"))
nopath = lambda ls: [re.sub(r"^plantfind \S*apn_8bit", "plantfind apn_8bit", l) for l in ls]
if have("s6-plantfind-c1.log"):
    report("plantfind c1: plants-bl22-c1.txt vs the original", lines(os.path.join(RR, "plants-bl22-c1.txt")) == lines(C6 + "plants-bl22-c1.txt"))
    report("plantfind c1: output vs plantfind-bl22-c1.out (paths aside)", nopath(body(os.path.join(RR, "s6-plantfind-c1.log"))) == nopath(lines(C6 + "plantfind-bl22-c1.out")))
if have("s6-plantfind-c2.log"):
    report("plantfind c2 (k >= 4): output vs plantfind-bl22-c2.out (paths aside)", nopath(body(os.path.join(RR, "s6-plantfind-c2.log"))) == nopath(lines(C6 + "plantfind-bl22-c2.out")))
    print("   plantings, k >= 4:", [l.split()[1] for l in lines(os.path.join(RR, "plants-bl22-c2.txt"))])
if have("s6-plantfind-c2k3.log"):
    report("plantfind c2 (k >= 3): plants-bl22-c2k3.txt vs the original plants-bl22-c2.txt", lines(os.path.join(RR, "plants-bl22-c2k3.txt")) == lines(C6 + "plants-bl22-c2.txt"))
    print("   ", [l for l in body(os.path.join(RR, "s6-plantfind-c2k3.log")) if l.startswith("plantfind")])


def records(path, tag):
    out, dup = {}, 0
    for l in open(path):
        if l.startswith(tag + " "):
            v = l.split()
            dup += int(v[1]) in out
            out[int(v[1])] = v[2:]
    return out, dup


def coreset(path, tag, confs):
    by = {}
    for l in open(path):
        if l.startswith(tag + " "):
            v = l.split()
            if int(v[1]) in confs:
                by.setdefault(int(v[1]), set()).add(" ".join(v[2:]))
    return by


for run, tag, ctag, orig in (("p1asc", "R", "K", "p1asc"), ("p1desc", "R", "K", "p1desc"), ("lin", "L", "C", "lin"),
                             ("lindesc", "L", "C", "lindesc"), ("asc", "R", None, "asc")):
    pf = os.path.join(RR, run + ".prog")
    if not os.path.exists(pf):
        continue
    a, dup = records(pf, tag)
    b, _ = records(C6 + orig + ".prog", tag)
    bad = [i for i in a if a[i] != b.get(i)]
    print(f"{run}: re-run {tag} lines {len(a)} (configurations {sorted(a)[:6]}{' ...' if len(a) > 6 else ''}, {dup} repeated), unequal to the original's {len(bad)} {bad[:8]}")
    for i in bad[:3]:
        print(f"   config {i}: fields differing {[k for k, (x, y) in enumerate(zip(a[i], b.get(i, []))) if x != y]}")
    if ctag:
        cf = os.path.join(RR, run + ".cores")
        if os.path.exists(cf):
            ca, cb = coreset(cf, ctag, set(a)), coreset(C6 + orig + ".cores", ctag, set(a))
            diff = [i for i in a if ca.get(i, set()) != cb.get(i, set())]
            print(f"   core files: {sum(len(v) for v in ca.values())} cores on those configurations, sets differing from the original's on {len(diff)} {diff[:8]}")
if have("s6-ctl-lind-bl22.log"):
    report("ctl lind-bl22: ctl-lind-bl22.txt vs the original (paths aside)",
           notime(lines(os.path.join(RR, "ctl-lind-bl22.txt"))) == notime(lines(C6 + "ctl-lind-bl22.txt")))
    a, b = notime(body(os.path.join(RR, "s6-ctl-lind-bl22.log"))), notime(lines(C6 + "ctl-lind-bl22.out"))
    report("ctl lind-bl22: output vs ctl-lind-bl22.out (timings and paths aside)", a == b,
           f"({len(a)} vs {len(b)} lines; only in the re-run {[x for x in a if x not in b]}, only in the original {[x for x in b if x not in a]})")
for k in (0, 3):
    if have(f"s6-ctlplant-c2bl-{k}.log"):
        strip = lambda ls: [re.sub(r"\s*\[\d+ s\]$", "", l) for l in ls if l.startswith(f"FOUND plant {k}") or l.startswith(f"NOTFOUND plant {k}")]
        a = strip(lines(os.path.join(RR, "ctlplant-c2bl-lin.txt")))
        b = strip(lines(C6 + "ctlplant-c2bl-lin.txt"))
        report(f"ctlplant c2bl-lin plant {k}", a == b and len(a) == 1, str(a))
if have("s6-ctlrot-lin.log"):
    a = [l for l in body(os.path.join(RR, "s6-ctlrot-lin.log")) if l.startswith("rep ")]
    keep = {l.split()[1] for l in a}
    b = [l for f in glob.glob(C6 + "ctlrot-lin-*.out") for l in lines(f) if l.startswith("rep ") and l.split()[1] in keep]
    report("ctlrot lin: rep lines vs ctlrot-lin-*.out", sorted(a) == sorted(b), str(a))
