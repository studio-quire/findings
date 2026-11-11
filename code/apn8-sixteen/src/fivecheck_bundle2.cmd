@echo off
rem Search B's classification, covers, block view and planted test again, with outputs kept, for five's bundle (2026-11-09).
rem configs: written to a new file and compared with the run's runs\check5\configs.txt. Run from the gleaner root.
set O=runs\check5\bundle
if not exist %O% mkdir %O%
if not exist runs\bin\fivecheck_configs_b.exe exit /b 2
runs\bin\fivecheck_configs_b.exe %O%\configs-rerun.txt > %O%\configs.txt 2>&1 || exit /b 1
fc /b %O%\configs-rerun.txt runs\check5\configs.txt > %O%\configs-fc.txt 2>&1 || exit /b 1
python apn\fivecheck_covers.py > %O%\covers.txt 2>&1 || exit /b 1
python apn\fivecheck_blockview.py runs\check5\known4.txt > %O%\blockview.txt 2>&1 || exit /b 1
runs\bin\fivecheck_classes_b.exe plant runs\check5\known4.txt 4 > %O%\plant.txt 2>&1 || exit /b 1
runs\bin\fivecheck_classes_b.exe plant runs\check5\known4.txt 4 pivot > %O%\plant-pivot.txt 2>&1 || exit /b 1
echo ALLDONE
