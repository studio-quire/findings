@echo off
rem The representatives carrying the known functions not yet reached in list order: same pipeline, order chosen with hindsight.
set OMP_NUM_THREADS=4
if not exist runs\check5\prio01 mkdir runs\check5\prio01
if not exist runs\check5\prio02 mkdir runs\check5\prio02
start /b "" cmd /c "runs\bin\fivecheck_ctl.exe list 4 1 runs/check5/prio_c01.txt runs/check5/prio01 pivot"
runs\bin\fivecheck_ctl.exe list 4 2 runs/check5/prio_c02.txt runs/check5/prio02 pivot || exit /b 1
echo ALLDONE
