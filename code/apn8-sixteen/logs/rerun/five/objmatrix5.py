"""objmatrix5.py (in runs/rerun5): which fresh build gives which run binary's object file, function by function.
`dumpbin /disasm` listings of objects compared with addresses and raw bytes dropped, MSVC's hashed lambda names and
symbol-relative offsets normalised (six's objcmp.py), order ignored. For each run object: the builds whose function bodies are
identical as a multiset, else the closest and the names of the functions that differ. Claude Opus 5.5, 2026-11-11."""
import collections, glob, os, re, sys

os.chdir(os.path.dirname(os.path.abspath(__file__)))


def norm(s):
    s = re.sub(r"<lambda_[0-9a-f]+>", "<lambda>", s)
    s = re.sub(r"\$[0-9a-f]{16,}", "$H", s)
    return re.sub(r"\s+", " ", s).strip()


def funcs(path):
    out, name, body = {}, None, []
    for l in open(path, errors="replace"):
        l = l.rstrip("\n")
        m = re.match(r"^(\S.*):$", l)
        if m and not l.startswith(("Dump of", "File Type")):
            if name is not None:
                out.setdefault(name, []).append(tuple(body))
            name, body = norm(m[1]), []
            continue
        m = re.match(r"^\s+[0-9A-F]{16}: (?:[0-9A-F]{2} )+\s*(.*)$", l)
        if m and name is not None:
            body.append(norm(m[1]))
    if name is not None:
        out.setdefault(name, []).append(tuple(body))
    return out


def bag(f):
    return collections.Counter(x for v in f.values() for x in v)


new = {os.path.basename(p)[:-4]: funcs(p) for p in sorted(glob.glob("dis/new-*.dis") + glob.glob("dis/hist-*.dis"))}
run = {os.path.basename(p)[:-4]: funcs(p) for p in sorted(glob.glob("dis/run-*.dis"))}
nb = {k: bag(v) for k, v in new.items()}
for r, fr in run.items():
    cr = bag(fr)
    hits = [n for n, cn in nb.items() if cn == cr]
    best = min(nb, key=lambda n: sum(((nb[n] - cr) + (cr - nb[n])).values()))
    line = f"{r}: {sum(cr.values())} functions; identical to {hits or 'none'}"
    if not hits:
        fb = new[best]
        diff = sorted(n for n in set(fr) | set(fb) if collections.Counter(fr.get(n, [])) != collections.Counter(fb.get(n, [])))
        line += (f"; closest {best} ({sum((nb[best] - cr).values())} bodies only in it, {sum((cr - nb[best]).values())} only in the "
                 f"run's); functions that differ: {', '.join(d[:60] for d in diff[:12])}{' ...' if len(diff) > 12 else ''}")
    print(line)
