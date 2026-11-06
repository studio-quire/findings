@echo off
rem cover_chunks_split.cmd <dir> <splitdim> [exe]: cover_chunks.cmd with runs\bin\<exe>.exe (default cover4, cover.cpp with
rem the split dimension, commit 5265bc8) run as "<exe> <chunk> -1 <splitdim>"; finished chunks (.out ends with
rem "done: covers") are skipped.
rem Rowen (Claude Opus 5.5), 2026-11-01, for i = 11 (6 rank-4); exe argument 2026-11-04 (cover5, the affine-tree hole stage).
cd /d %~dp0..
set EXE=%3
if "%EXE%"=="" set EXE=cover4
for %%f in (%1\c*.txt) do call :one "%%~dpnf" %2
exit /b 0
:one
findstr /b /c:"done: covers" "%~1.out" >nul 2>nul && exit /b 0
runs\bin\%EXE%.exe "%~1.txt" -1 %2 > "%~1.out" 2>&1
if errorlevel 1 exit 1
exit /b 0
