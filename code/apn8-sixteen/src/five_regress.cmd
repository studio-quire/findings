@echo off
rem five_regress.cmd: --frange --p1all slices against six's whole per-configuration run, on p6-101 (dim-4, 5 and 6 cores).
rem Rowen (Claude Opus 5.5), 2026-11-06.
cd /d %~dp0..
runs\bin\cover5u.exe runs\i11\p6.txt -1 5 --only p6-101-a24-w11 --p1all --frange 0 10000 2>&1 | findstr /v /b NEAR > runs\i12\reg\s0.out
runs\bin\cover5u.exe runs\i11\p6.txt -1 5 --only p6-101-a24-w11 --p1all --frange 10000 20000 2>&1 | findstr /v /b NEAR > runs\i12\reg\s1.out
runs\bin\cover5u.exe runs\i11\p6.txt -1 5 --only p6-101-a24-w11 --p1all --frange 20000 30000 2>&1 | findstr /v /b NEAR > runs\i12\reg\s2.out
exit /b 0
