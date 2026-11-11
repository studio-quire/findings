"""Search A on five (i = 12): cover5v over runs/i12/p5.txt, one process per slice of one configuration's frontier.

  python apn/five_slices.py run [SLICE]    run every slice not finished yet; safe to kill and start again
  python apn/five_slices.py tally          sum finished slices; check they tile each frontier exactly

Days of work on a box that evicts without warning, so a slice (default 10,000 items, ~15 min in a shared room) is the
unit of loss: its output is written to .tmp and renamed only after exit code 0 and a "done: covers" line. The frontier is
sorted before --frange cuts it, so slices of separate processes partition it; each slice checks the frontier size it saw
against the first slice's. Every process repeats phase 1 and prints the configuration's kept cores as NEAR lines: those
are kept once, gzipped, from the slice at item 0, and dropped from the rest. NEAR lines after phase 1 are
the slice's own and stay in its output.

--p1all is not optional. Plain --frange skips the hole stage in phase 1, so a core below dimension 5 (five forms spanning 4
dimensions) never reaches the frontier and its whole subtree is silently lost; on six seven such cores held 43% of the
candidates. --frange also resets the counters after phase 1, so the tally adds each configuration's "phase 1 counts" once.

--reps searches one core per Aut(P)-orbit (five_aut.py, five_orbits.cpp, five_cores.cmd): each S8 over P has exactly one core, the
span of its rank-4 forms, and g in Aut(P) carries a solution on one core to a solution on its image, so an orbit is settled by any
one member. Counts here are therefore of the reduced search, not comparable to a whole run's.

Anything but a clean empty slice stops the run with exit code 3, so the runner reports a failed job and wakes me: a
SOLUTION here would answer Open Problem VI.1, and FULL, BADSEED or non-APN lines mean something to look at, not to pass.
Rowen (Claude Opus 5.5), 2026-11-06.
"""
import glob, gzip, os, re, subprocess, sys

EXE = r"runs\bin\cover5v.exe"   # cover5t plus --p1all, --reps and the phase 1 counts line; the search below a core is unchanged
CONFIGS = r"runs\i12\p5.txt"
OUT = "runs/i12/cover"
SPLIT = "5"
STOP = "runs/i12/cover/STOPPED.txt"
ALARM = ("SOLUTION", "FULL", "BADSEED")


def tags():
    return [l.split()[0] for l in open(CONFIGS) if l.strip()]


def reps_path(tag):
    return rf"runs\i12\reps\{tag}.txt"


def path(tag, a, b):
    return f"{OUT}/{tag}/{a:08d}-{b:08d}.out"


def clean(summary):
    m = re.search(r"solutions (\d+), full \(dim < 8\) (\d+), dim 8 non-APN (\d+), bad seeds (\d+)", summary)
    return m and not any(int(x) for x in m.groups())


def run_slice(tag, a, b, size_seen):
    """Returns the frontier size this process saw."""
    os.makedirs(f"{OUT}/{tag}", exist_ok=True)
    final = path(tag, a, b)
    tmp, near_tmp = final + ".tmp", f"{OUT}/{tag}/near-phase1.txt.gz.tmp"
    keep_near = a == 0
    p = subprocess.Popen([EXE, CONFIGS, "-1", SPLIT, "--only", tag, "--p1all", "--reps", reps_path(tag), "--frange", str(a), str(b)],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace")
    size, phase1, alarm, done, p1near = None, True, [], "", None
    near = gzip.open(near_tmp, "wt", compresslevel=6) if keep_near else None
    with open(tmp, "w") as o:
        for l in p.stdout:
            if phase1 and l.startswith("NEAR "):
                if near:
                    near.write(l)
                continue
            if phase1:
                m = re.match(r"phase 1 done \d+ s: frontier (\d+) items", l)
                if m:
                    size, phase1 = int(m[1]), False
            m = re.match(r"phase 1 counts: .*near misses (\d+),", l)
            if m:
                p1near = int(m[1])
            if l.startswith(ALARM):
                alarm.append(l)
            if l.startswith("done: covers"):
                done = l
            o.write(l)
    if near:
        near.close()
    rc = p.wait()
    if rc or not done or size is None:
        sys.exit(f"{tag} [{a}, {b}): exit {rc}, no complete output (kept {tmp})")
    # a key mismatch in --reps would drop every core and look like an empty search: phase 1 must keep exactly one core per orbit
    want = sum(1 for _ in open(reps_path(tag)))
    if p1near != want:
        stop(f"{tag} [{a}, {b}): phase 1 kept {p1near} cores, {want} orbits\n")
    if size_seen is not None and size != size_seen:
        stop(f"{tag} [{a}, {b}): frontier {size}, first slice saw {size_seen}\n")
    # before the rename: a kill between the two must not leave an alarmed slice final with no STOPPED.txt (the .tmp stays as evidence)
    if alarm or not clean(done):
        stop(f"{tag} [{a}, {b}): {done}" + "".join(alarm))
    os.replace(tmp, final)
    if keep_near:
        os.replace(near_tmp, f"{OUT}/{tag}/near-phase1.txt.gz")
    return size


def stop(text):
    with open(STOP, "a") as f:
        f.write(text)
    print("STOPPED:", text, flush=True)
    sys.exit(3)


def frontier_size(tag):
    for f in sorted(glob.glob(f"{OUT}/{tag}/*.out")):
        for l in open(f, errors="replace"):
            m = re.match(r"phase 1 done \d+ s: frontier (\d+) items", l)
            if m:
                return int(m[1])
    return None


def cmd_run(step):
    os.makedirs(OUT, exist_ok=True)
    if os.path.exists(STOP):
        sys.exit(f"{STOP} exists: read it first")
    for tag in tags():
        size = frontier_size(tag)
        if size is None:
            size = run_slice(tag, 0, step, None)
            print(tag, "frontier", size, flush=True)
        for a in range(0, size, step):
            b = min(a + step, size)
            if a == 0:
                b = step   # the first slice's name, fixed before the size was known
            if not os.path.exists(path(tag, a, b)):
                run_slice(tag, a, b, size)
        print(tag, "done", flush=True)


def cmd_tally():
    tot = dict(nodes=0, cand=0, passed=0, sols=0, near2=0)
    dims = [0] * 8
    missing, items, whole = 0, 0, 0
    for tag in tags():
        size, covered, first = None, 0, True
        orbits = sum(1 for _ in open(reps_path(tag)))
        for f in sorted(glob.glob(f"{OUT}/{tag}/*.out")):
            a, b = map(int, re.search(r"(\d+)-(\d+)\.out$", f).groups())
            t = open(f, errors="replace").read()
            s = int(re.search(r"phase 1 done \d+ s: frontier (\d+) items", t)[1])
            assert size in (None, s), (f, size, s)
            size = s
            assert a == covered, (f, "expected start", covered)
            covered = min(b, s)
            assert not re.search(r"^(SOLUTION|FULL|BADSEED)", t, re.M), f
            m = re.search(r"phase 1 counts: nodes (\d+), by dimension ([\d ]+), near misses (\d+), candidates (\d+), passed (\d+)", t)
            assert int(m[3]) == orbits, (f, m[3], orbits)
            if first:   # phase 1 is repeated by every slice and identical: count it once
                first = False
                tot["nodes"] += int(m[1])
                dims = [x + int(y) for x, y in zip(dims, m[2].split())]
                tot["cand"] += int(m[4])
                tot["passed"] += int(m[5])
            assert clean(re.search(r"^done: covers.*$", t, re.M)[0]), f
            tot["nodes"] += int(re.search(r"^done: covers \d+, tasks \d+, nodes (\d+)", t, re.M)[1])
            d = list(map(int, re.search(r"nodes by span dimension 1\.\.8: (.*)", t)[1].split()))
            dims = [x + y for x, y in zip(dims, d)]
            c = re.search(r"hole-stage candidates (\d+), passing the filter \(sent to push\) (\d+)", t)
            tot["cand"] += int(c[1])
            tot["passed"] += int(c[2])
            tot["near2"] += len(re.findall(r"^NEAR ", t, re.M))
        items += covered
        if size is not None and covered == size:
            whole += 1
        else:
            missing += 1
    print(f"configurations finished {whole} of {whole + missing}; frontier items searched {items}")
    print("nodes by span dimension 1..8:", *dims)
    print(tot)
    if os.path.exists(STOP):
        print(f"{STOP} exists: not clean")
    if missing or os.path.exists(STOP):
        sys.exit(1)   # the verdict rests on this exit code, not on the run job's: only a complete, clean tally exits 0
    print("complete: every frontier tiled, every slice clean, every slice kept exactly the orbit representatives")


if __name__ == "__main__":
    if sys.argv[1] == "run":
        cmd_run(int(sys.argv[2]) if len(sys.argv) > 2 else 10000)
    else:
        cmd_tally()
