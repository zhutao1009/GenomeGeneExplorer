@echo off
rem ===========================================================================
rem  Genome Gene Explorer - launcher
rem
rem  NOTE: this file must keep CRLF line endings and stay ASCII-only.
rem  cmd.exe mis-parses LF-only batch files (especially inside "if (...)"
rem  blocks), which produces errors such as 'cho' is not recognized.
rem ===========================================================================
setlocal EnableExtensions EnableDelayedExpansion
title Genome Gene Explorer
cd /d "%~dp0"
if errorlevel 1 goto :no_dir

call "%~dp0_find_python.bat"

if defined GGE_PY goto :launch
if not defined GGE_PY_BASE goto :need_python

rem -- A Python 3.10+ exists but lacks PySide6: build a private venv from it --
set "PYB=!GGE_PY_PIP!"
if not defined PYB set "PYB=!GGE_PY_BASE!"

echo.
echo [setup] No ready-to-use environment found. Building .venv from:
echo         !PYB!
echo.
if not exist ".venv\Scripts\python.exe" "!PYB!" -m venv .venv
if not exist ".venv\Scripts\python.exe" goto :venv_failed

set "GGE_PY=%~dp0.venv\Scripts\python.exe"
"!GGE_PY!" -m pip --version >nul 2>nul
if errorlevel 1 (
    echo [setup] pip is missing in the venv, running ensurepip ...
    "!GGE_PY!" -m ensurepip --upgrade
)

echo [setup] Installing requirements - needs internet, about 200 MB, a few minutes.
"!GGE_PY!" -m pip install -r requirements.txt
if errorlevel 1 goto :pip_failed

:launch
echo.
echo [run] !GGE_PY!
"!GGE_PY!" "%~dp0tools\doctor.py" --quiet
"!GGE_PY!" "%~dp0run_app.py" %*
set "RC=%ERRORLEVEL%"
if not "%RC%"=="0" goto :app_failed
exit /b 0


:no_dir
echo.
echo [ERROR] Cannot switch to the program directory:
echo         %~dp0
echo If the folder is on a disconnected drive, reconnect it and try again.
pause
exit /b 1


:need_python
echo.
echo [ERROR] No Python 3.10 or newer was found on this computer.
echo.
echo   Option 1 - install Python (recommended):
echo     https://www.python.org/downloads/   (3.10 - 3.13, 64 bit)
echo     Tick "Add python.exe to PATH" during installation, then run start.bat
echo     again.
echo.
echo   Option 2 - use the prebuilt program, no Python needed:
echo     dist\GenomeGeneExplorer\GenomeGeneExplorer.exe
echo.
echo   Option 3 - Python lives somewhere unusual:
echo     set GGE_PYTHON=C:\path\to\python.exe
echo     start.bat
echo.
pause
exit /b 1


:venv_failed
echo.
echo [ERROR] Could not create a virtual environment with:
echo         !PYB!
echo.
echo   Try the prebuilt program instead - no Python needed:
echo     dist\GenomeGeneExplorer\GenomeGeneExplorer.exe
echo.
pause
exit /b 1


:pip_failed
echo.
echo [ERROR] pip could not install the dependencies.
echo   Likely no internet access, or a proxy / mirror problem.
echo.
echo   Retry behind a mirror:
echo     "!GGE_PY!" -m pip install -r requirements.txt -i https://pypi.tuna.tsinghua.edu.cn/simple
echo.
echo   Or use the prebuilt program - no Python needed:
echo     dist\GenomeGeneExplorer\GenomeGeneExplorer.exe
echo.
pause
exit /b 1


:app_failed
echo.
echo [ERROR] Genome Gene Explorer exited with code %RC%.
echo   Full diagnostic report:  "!GGE_PY!" tools\doctor.py
echo.
pause
exit /b %RC%
