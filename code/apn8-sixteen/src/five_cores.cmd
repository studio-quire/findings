@echo off
rem five_cores.cmd: each of five's configurations' cores (phase 1 only, NEAR lines, runs\i12\near\), then their orbits under Aut(P)
rem (runs\i12\reps\). five_orbits exits non-zero on any failed check, and this stops there. Rowen (Claude Opus 5.5), 2026-11-06.
cd /d %~dp0..
if not exist runs\i12\near mkdir runs\i12\near
if not exist runs\i12\reps mkdir runs\i12\reps
for /f "tokens=1" %%t in (runs\i12\p5.txt) do (
  runs\bin\cover5v.exe runs\i12\p5.txt -1 5 --only %%t --p1all --frange 0 0 > runs\i12\near\%%t.txt 2>&1
  if errorlevel 1 exit 1
  runs\bin\five_orbits.exe runs\i12\near\%%t.txt runs\i12\aut.txt %%t runs\i12\reps\%%t.txt >> runs\i12\reps\summary.txt 2>&1
  if errorlevel 1 exit 1
)
exit /b 0
