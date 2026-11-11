@echo off
rem The sets S of dimension 4 (five forms with one linear relation), Gamma-orbit representatives, all classes. T has dimension 4 here.
set OMP_NUM_THREADS=6
if not exist runs\check5\k4 mkdir runs\check5\k4
for /L %%c in (10,1,31) do call :one %%c
echo ALLDONE
exit /b 0
:one
set /a C=%1-10
set N=%1
if %C% LSS 10 (set F=runs/check5/k4_c0%C%.txt) else (set F=runs/check5/k4_c%C%.txt)
runs\bin\fivecheck_ctl.exe list 5 %C% %F% runs/check5/k4 pivot || exit 1
exit /b 0
