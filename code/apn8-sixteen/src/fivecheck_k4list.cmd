@echo off
set OMP_NUM_THREADS=6
runs\bin\fivecheck_cl2.exe k4list || exit /b 1
echo ALLDONE
