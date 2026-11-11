"""Is five's obstruction local? For a sample of search A's core representatives (dim 5) from every configuration, run
five_critical (Fable 5.2's critical-problem solver): the whole problem (expect UNSAT, as both searches found) and the
subproblem on the planes inside each hyperplane of V (Fable found all 255 SAT on three cores). Rowen (Claude Opus 5.5),
2026-11-09.

python apn/five_hyp.py OUT_DIR [per_config] [threads]     resumable: finished cores are skipped
python apn/five_hyp.py OUT_DIR --tally
"""
import concurrent.futures as cf
import glob
import os
import random
import re
import subprocess
import sys

EXE = os.path.join('runs', 'bin', 'five_critical.exe')


def from_coef(c):
    B, k = 0, 0
    for i in range(8):
        for j in range(i + 1, 8):
            if c >> k & 1:
                B |= 1 << (8 * i + j)
                B |= 1 << (8 * j + i)
            k += 1
    return B


def sample(per):
    out = []
    for path in sorted(glob.glob(os.path.join('runs', 'i12', 'reps', 'p5-*.txt'))):
        tag = os.path.basename(path)[:-4]
        lines = [l.split() for l in open(path, encoding='utf-8') if l.startswith('5 ')]
        rng = random.Random('five-hyp ' + tag)
        for idx in sorted(rng.sample(range(len(lines)), min(per, len(lines)))):
            forms = ['%016x' % from_coef(int(x, 16)) for x in lines[idx][1:6]]
            out.append((tag, idx, forms))
    return out


def run_one(item, outdir):
    tag, idx, forms = item
    path = os.path.join(outdir, f'{tag}_{idx}.txt')
    if os.path.exists(path):
        return path
    full = subprocess.run([EXE] + forms + ['1', '0'], capture_output=True, text=True).stdout
    hyp = subprocess.run([EXE] + forms + ['1', '0', '-1', '1'], capture_output=True, text=True).stdout
    keep = [l for l in full.splitlines() if l.startswith(('G5', 'iter', 'SAT', 'cap', 'core rank', 'no 3-dim'))]
    keep += [l for l in hyp.splitlines() if l.startswith(('sat ', 'H ')) and ('UNSAT' in l or 'cap' in l or l.startswith('sat '))]
    tmp = path + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as f:
        f.write(f'core {tag} line {idx} forms {" ".join(forms)}\n' + '\n'.join(keep) + '\n')
    os.replace(tmp, path)
    return path


def tally(outdir):
    n = unsat = sat_full = cap_full = 0
    hyp_bad, cores, g5 = [], [], []
    for path in sorted(glob.glob(os.path.join(outdir, 'p5-*.txt'))):
        text = open(path, encoding='utf-8').read()
        n += 1
        if 'SAT!' in text:
            sat_full += 1
            print('SAT (an APN S8 over this core!):', path)
        elif re.search(r'^cap', text, re.M):
            cap_full += 1
        else:
            unsat += 1
            cs = [int(m) for m in re.findall(r'-> core (\d+)', text)]
            if cs:
                cores.append(cs[-1])
        g = re.search(r'^G5 (\d+)', text, re.M)
        if g:
            g5.append(int(g.group(1)))
        h = re.search(r'^sat (\d+) unsat (\d+) cap (\d+)', text, re.M)
        if not h or h.group(1) != '255':
            hyp_bad.append((path, h.group(0) if h else 'no hyperplane line'))
    print(f'cores {n}: whole problem UNSAT {unsat}, SAT {sat_full}, cap {cap_full}')
    if g5:
        print(f'|G5| {min(g5)}..{max(g5)}')
    if cores:
        print(f'conflict sets after the first passes (not minimised): {min(cores)}..{max(cores)} planes')
    print(f'cores with every hyperplane subproblem SAT: {n - len(hyp_bad)} of {n}')
    for p, h in hyp_bad:
        print('  not all hyperplanes SAT:', os.path.basename(p), h, '|', local_reason(p))


def apply(f, a):
    s = 0
    for i in range(8):
        if a >> i & 1:
            s ^= (f >> (8 * i)) & 255
    return s


def local_reason(path):
    """For a core with UNSAT hyperplanes: E = the subspace common to them; whether every form of S5 vanishes on E (the 4-space
    lemma in README "Toward a proof"), and how E meets the radicals of S5's rank-4 elements."""
    text = open(path, encoding='utf-8').read()
    hs = [int(x) for x in re.findall(r'^H (\d+) 0 .*UNSAT', text, re.M)]
    E = [v for v in range(1, 256) if all(bin(v & h).count('1') % 2 == 0 for h in hs)]
    span = [0]
    for f in (int(x, 16) for x in re.search(r'forms (.*)', text).group(1).split()):
        span += [g ^ f for g in span]
    iso = all(bin(apply(f, a) & b).count('1') % 2 == 0 for f in span for a in E for b in E)
    rads = [set(v for v in range(1, 256) if apply(f, v) == 0) for f in span[1:]]
    meets = sorted(len(set(E) & r) for r in rads if len(r) == 15)
    dimE = (len(E) + 1).bit_length() - 1
    return f'common subspace dim {dimE}, every form of S5 vanishes on it: {iso}, points shared with each rank-4 radical {meets}'


if __name__ == '__main__':
    outdir = sys.argv[1]
    if len(sys.argv) > 2 and sys.argv[2] == '--tally':
        tally(outdir)
        sys.exit(0)
    per = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    threads = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    os.makedirs(outdir, exist_ok=True)
    items = sample(per)
    print(f'{len(items)} cores sampled', flush=True)
    with cf.ThreadPoolExecutor(threads) as ex:
        for i, _ in enumerate(ex.map(lambda it: run_one(it, outdir), items)):
            if (i + 1) % 50 == 0:
                print(f'{i + 1} done', flush=True)
    tally(outdir)
