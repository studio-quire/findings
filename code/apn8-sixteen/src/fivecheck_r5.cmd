@echo off
rem Five-check main run: every set S of dimension 5, in all 22 classes, resumable by slice (a slice = one choice of the first form).
set OMP_NUM_THREADS=12
if not exist runs\check5\r5 mkdir runs\check5\r5
for /L %%c in (0,1,21) do (
  runs\bin\fivecheck_r5.exe run 5 %%c 0 28 runs/check5/r5 onlyk=5 || exit /b 1
)
echo ALLDONE
