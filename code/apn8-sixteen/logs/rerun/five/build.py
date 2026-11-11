"""Five's re-run builds and the provenance check's disassembly (run from the gleaner root, as a share-runner job).

bin/<name>.exe: the bundle's src/<name>.cpp, compiled from src/ under its own name with the Build line (and src/r5/ into bin/r5/).
hist/<label>/: sources as committed at the commits the run binaries were built from, compiled under their own names.
dis/<label>.dis: dumpbin /disasm of every object, fresh and kept (runs/bin/<exe>_.obj, the reviewer's in runs/check5/review/).
Claude Opus 5.5, for the bundle, 2026-11-11.
"""
import os, shutil, subprocess, sys

ROOT = os.getcwd()
R = os.path.join(ROOT, 'runs', 'rerun5')
SRC = os.path.join(ROOT, 'findings', 'code', 'apn8-sixteen', 'src')
VC = r'C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat'
FLAGS = ['/nologo', '/O2', '/openmp', '/EHsc', '/std:c++17', '/arch:AVX2']

out = subprocess.run(f'cmd /c ""{VC}" >nul 2>nul && set"', capture_output=True, text=True, shell=True).stdout
env = dict(l.split('=', 1) for l in out.splitlines() if '=' in l)
log = open(os.path.join(R, 'build.log'), 'w')


PATH = next(v for k, v in env.items() if k.upper() == 'PATH')


def run(cmd, cwd, **kw):
    exe = shutil.which(cmd[0], path=PATH)   # CreateProcess searches the parent's PATH, not env's
    r = subprocess.run([exe] + cmd[1:], cwd=cwd, env=env, capture_output=True, text=True, errors='replace', **kw)
    log.write(f'$ {" ".join(cmd)}  (in {os.path.relpath(cwd, ROOT)})\n{r.stdout}{r.stderr}# exit {r.returncode}\n')
    log.flush()
    if r.returncode:
        sys.exit(f'failed: {cmd}')
    return r


def cl(cwd, src, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    run(['cl'] + FLAGS + [src, f'/Fe{outdir}\\{name}.exe', f'/Fo{outdir}\\{name}_'], cwd)


def git_show(commit, path, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    data = subprocess.run(['git', 'show', f'{commit}:{path}'], cwd=ROOT, capture_output=True, check=True).stdout
    open(dst, 'wb').write(data)


BIN = os.path.join(R, 'bin')
for n in ('cover', 'pspread3', 'five_orbits', 'five_vs_b', 'five_critical', 'fivecheck_configs', 'fivecheck_sieve',
          'fivecheck_classes', 'rv5_indep', 'rvb5_ref', 'rvb5_gamma', 'rvb5_gamma5', 'defcheck', 'defcheck2'):
    cl(SRC, f'{n}.cpp', BIN, n)
cl(os.path.join(SRC, 'r5'), 'fivecheck_sieve.cpp', os.path.join(BIN, 'r5'), 'fivecheck_sieve')

# label: (commit, [files], main file)
HIST = {
    'cover-10d6aaf': ('10d6aaf', ['apn/cover.cpp'], 'cover.cpp'),
    'cover-cfc3795': ('cfc3795', ['apn/cover.cpp'], 'cover.cpp'),
    'orbits-cfc3795': ('cfc3795', ['apn/five_orbits.cpp'], 'five_orbits.cpp'),
    'sieve-bff5505': ('bff5505', ['apn/fivecheck_sieve.cpp'], 'fivecheck_sieve.cpp'),
    'classes-bff5505': ('bff5505', ['apn/fivecheck_classes.cpp', 'apn/fivecheck_sieve.cpp'], 'fivecheck_classes.cpp'),
    'classes-a238b03': ('a238b03', ['apn/fivecheck_classes.cpp', 'apn/fivecheck_sieve.cpp'], 'fivecheck_classes.cpp'),
    'classes-3f8c590': ('3f8c590', ['apn/fivecheck_classes.cpp', 'apn/fivecheck_sieve.cpp'], 'fivecheck_classes.cpp'),
}
for label, (commit, files, main) in HIST.items():
    d = os.path.join(R, 'hist', label)
    for f in files:
        git_show(commit, f, os.path.join(d, os.path.basename(f)))
    cl(d, main, d, 'h')
# five_vs_b as it ran: apn/five_vs_b.cpp from apn/, including the reviewer's file by its working-tree path
cl(os.path.join(ROOT, 'apn'), 'five_vs_b.cpp', os.path.join(R, 'hist', 'vsb-apn'), 'h')

DIS = os.path.join(R, 'dis')
os.makedirs(DIS, exist_ok=True)


def dis(obj, label):
    r = subprocess.run([shutil.which('dumpbin', path=PATH), '/nologo', '/disasm', obj], env=env, capture_output=True, text=True,
                       errors='replace')
    open(os.path.join(DIS, label + '.dis'), 'w').write(r.stdout)


for f in sorted(os.listdir(BIN)):
    if f.endswith('_.obj'):
        dis(os.path.join(BIN, f), 'new-' + f[:-5])
dis(os.path.join(BIN, 'r5', 'fivecheck_sieve_.obj'), 'new-r5-fivecheck_sieve')
for label in list(HIST) + ['vsb-apn']:
    dis(os.path.join(R, 'hist', label, 'h_.obj'), 'hist-' + label)
RUNBIN = os.path.join(ROOT, 'runs', 'bin')
for n in ('cover5u', 'cover5v', 'pspread3h', 'five_orbits', 'five_vs_b', 'five_critical', 'fivecheck_configs', 'fivecheck_configs_b',
          'fivecheck_r5', 'fivecheck_ctl', 'fivecheck_sieve', 'fivecheck_sieve_b', 'fivecheck_cl2', 'fivecheck_cl3',
          'fivecheck_classes', 'fivecheck_classes_b', 'fivecheck_frozencl'):
    dis(os.path.join(RUNBIN, n + '_.obj'), 'run-' + n)
REV = os.path.join(ROOT, 'runs', 'check5', 'review')
for n in ('rvb5_ref', 'rvb5_gamma', 'rvb5_ref2', 'rvb5_gamma5'):
    dis(os.path.join(REV, f'bin_{n}_.obj'), 'run-' + n)
print('built')
