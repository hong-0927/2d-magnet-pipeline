@echo off
setlocal
rem 2d-magnet.cmd -- 2D magnet pipeline entry (run_2d_magnet.py --quiet)
rem Usage: 2d-magnet.cmd --formula MnI2 --element Mn --nonmag I --S 2.5 --moment 5.0 --U 4.0 --J1 -0.3
rem Pipeline root = two levels up from this file (..\..\ = <repo>).
rem Python discovery order: %MAGNET_PYTHON% -> python -> py -3 -> bundled fallback.
set "ROOT=%~dp0..\.."
set "PY="
if defined MAGNET_PYTHON set "PY=%MAGNET_PYTHON%"
if not defined PY for %%I in (python.exe) do if not defined PY set "PY=%%~$PATH:I"
if not defined PY for %%I in (py.exe) do if not defined PY set "PY=%%~$PATH:I"
if not defined PY if exist "C:\ProgramData\anaconda3\python.exe" set "PY=C:\ProgramData\anaconda3\python.exe"
if not defined PY (
  echo [2d-magnet] cannot find Python. Install Python 3.9+ or set MAGNET_PYTHON.
  exit /b 1
)
"%PY%" "%ROOT%\run_2d_magnet.py" --quiet %*
exit /b %errorlevel%
