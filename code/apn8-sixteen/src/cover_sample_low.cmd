@echo off
rem cover_sample_low.cmd: the same comparison as cover_sample.cmd on the seven dimension-4 frontier items, which the sample's split at 5 and 6
rem can't draw: all of them (they sit on p6-101, p6-114, p6-304, p6-310, runs\cover5\sample\low4.txt). Rowen (Claude Opus 5.5), 2026-11-04.
cd /d %~dp0..
set B=runs\bin\cover5s.exe
set I=runs\cover5\sample\low4.txt
set O=runs\cover5\sample
call :one tree-s4 "4 --fsample 100 41"
call :one old-s4 "4 --fsample 100 41 --old"
exit /b 0
:one
findstr /b /c:"done: covers" %O%\%1.out >nul 2>nul && exit /b 0
%B% %I% -1 %~2 | findstr /v /b "NEAR" > %O%\%1.out
exit /b 0
