@echo off
rem Search B's definition checks again, with logs kept, for five's code bundle (2026-11-09; the first runs' logs weren't saved).
rem Binaries from the committed source: apn\build.cmd fivecheck_sieve fivecheck_sieve_b, apn\build.cmd fivecheck_classes fivecheck_classes_b.
rem Run from the gleaner root.
set O=runs\check5\bundle
if not exist %O% mkdir %O%
runs\bin\fivecheck_classes_b.exe gamma 4 > %O%\gamma4.txt 2>&1 || exit /b 1
runs\bin\fivecheck_classes_b.exe slevel 5 2000 1109 > %O%\slevel5.txt 2>&1 || exit /b 1
runs\bin\fivecheck_classes_b.exe slevel 4 2000 1109 > %O%\slevel4.txt 2>&1 || exit /b 1
for %%c in (0 2 4 8 10 13 14 19 20 21) do (
  runs\bin\fivecheck_sieve_b.exe flat 5 %%c 4 %%c909 > %O%\flat5_c%%c.txt 2>&1 || exit /b 1
)
for %%c in (0 1 2) do (
  runs\bin\fivecheck_sieve_b.exe flat 4 %%c 4 %%c909 > %O%\flat4_c%%c.txt 2>&1 || exit /b 1
)
runs\bin\fivecheck_classes_b.exe gamma 5 > %O%\gamma5.txt 2>&1 || exit /b 1
echo ALLDONE
