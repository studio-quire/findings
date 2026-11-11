@echo off
rem Control: r = 4, all Gamma-orbit representatives of form tuples in one class of four solids (argument: class index, two digits).
set OMP_NUM_THREADS=8
if not exist runs\check5\c4c%1 mkdir runs\check5\c4c%1
runs\bin\fivecheck_ctl.exe list 4 %2 runs/check5/reps_r4_c%1.txt runs/check5/c4c%1 pivot || exit /b 1
echo ALLDONE
