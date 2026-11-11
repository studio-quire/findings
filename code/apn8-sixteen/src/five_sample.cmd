@echo off
rem five_sample.cmd <n> <seed>: cover5s --fsample on each of five's configurations in turn, one process each, so only one
rem frontier is in memory; NEAR lines dropped. A file is complete only if it has a "done: covers" line (the pipe's exit code
rem is findstr's). Rowen (Claude Opus 5.5), 2026-11-06.
cd /d %~dp0..
for /f "tokens=1" %%t in (runs\i12\p5.txt) do call :one %%t %1 %2
exit /b 0
:one
findstr /b /c:"done: covers" "runs\i12\prof\s-%1.out" >nul 2>nul && exit /b 0
runs\bin\cover5u.exe runs\i12\p5.txt -1 5 --only %1 --p1all --fsample %2 %3 2>&1 | findstr /v /b NEAR > "runs\i12\prof\s-%1.out"
exit /b 0
