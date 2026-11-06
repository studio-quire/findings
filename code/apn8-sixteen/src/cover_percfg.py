"""Search A for six again, one configuration at a time, for per-configuration counts: the published cover5t binary, same
split, with --only TAG, its output kept without the NEAR lines. review-six.md item 5: on the configurations where A and B
share a frame, the two final hole stages should be one tree, which only per-configuration counts can show.

  python apn/cover_percfg.py [out dir, default runs/cover5/percfg]
Skips configurations whose output already has its `done:` line, so a requeued job carries on. Rowen (Claude Opus 5.5), 2026-11-05.
"""
import glob
import os
import subprocess
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "runs/cover5/percfg"
os.makedirs(OUT, exist_ok=True)


def finished(path):
    try:
        return any(l.startswith("done:") for l in open(path))
    except OSError:
        return False


todo = [(f, l.split()[0]) for f in sorted(glob.glob("runs/i11/cover5/c?.txt")) for l in open(f) if l.strip()]
print(f"{len(todo)} configurations", flush=True)
for n, (chunk, tag) in enumerate(todo):
    path = f"{OUT}/{tag}.out"
    if finished(path):
        continue
    p = subprocess.Popen([os.path.abspath("runs/bin/cover5t.exe"), chunk, "-1", "5", "--only", tag], stdout=subprocess.PIPE, text=True)
    with open(path + ".part", "w") as f:
        for line in p.stdout:
            if not line.startswith("NEAR "):
                f.write(line)
    if p.wait():
        sys.exit(f"{tag}: exit {p.returncode}")
    os.replace(path + ".part", path)
    print(f"{n + 1}/{len(todo)} {tag}", flush=True)
print("all done", flush=True)
