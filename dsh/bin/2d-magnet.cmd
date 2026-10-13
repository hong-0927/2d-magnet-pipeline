@echo off
setlocal
rem 2d-magnet.cmd -- 2D magnet pipeline entry (run_2d_magnet.py --quiet)
rem Usage: 2d-magnet.cmd --formula MnI2 --element Mn --nonmag I --S 2.5 --moment 5.0 --U 4.0 --J1 -0.3
rem Pipeline root = two levels up from this file (..\..\ = <repo>).
set "ROOT=%~dp0..\.."
set "PY=C:\ProgramData\anaconda3\python.exe"
if not exist "%PY%" set "PY=python"
"%PY%" "%ROOT%\run_2d_magnet.py" --quiet %*
exit /b %errorlevel%
