@echo off
rem cover_sample.cmd: the filter-free hole stage (--old) against the affine tree, on a random sample of every configuration
rem outside chunks c0, c1 (whose whole runs were compared already): 8 dim-6 frontier items and 1 dim-5 one per configuration, the
rem same items in both modes (fixed seeds). ITEM lines compared by apn/cover_sample_cmp.py. Rowen (Claude Opus 5.5), 2026-11-04.
cd /d %~dp0..
set B=runs\bin\cover5s.exe
set I=runs\cover5\c2-7.txt
set O=runs\cover5\sample
if not exist %O% mkdir %O%
call :one tree-s6 "6 --fsample 8 61"
call :one old-s6 "6 --fsample 8 61 --old"
call :one tree-s5 "5 --fsample 1 51"
call :one old-s5 "5 --fsample 1 51 --old"
exit /b 0
:one
findstr /b /c:"done: covers" %O%\%1.out >nul 2>nul && exit /b 0
%B% %I% -1 %~2 | findstr /v /b "NEAR" > %O%\%1.out
exit /b 0
