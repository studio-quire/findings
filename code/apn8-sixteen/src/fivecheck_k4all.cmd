@echo off
rem Invariance check of the Gamma reduction used for the dimension-4 sets: sieve all of them, no assembly.
set OMP_NUM_THREADS=6
set FC_NOASM=1
if not exist runs\check5\k4all mkdir runs\check5\k4all
runs\bin\fivecheck_cl3.exe k4list || exit /b 1
for /L %%c in (10,1,31) do call :one %%c
echo ALLDONE
exit /b 0
:one
set /a C=%1-10
if %C% LSS 10 (set F=runs/check5/k4all_c0%C%.txt) else (set F=runs/check5/k4all_c%C%.txt)
if exist %F% runs\bin\fivecheck_ctl.exe list 5 %C% %F% runs/check5/k4all pivot || exit 1
exit /b 0
