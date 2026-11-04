@echo off
rem ===========================================================================
rem  Genome Gene Explorer - build the Windows executable (directory mode)
rem
rem  NOTE: keep CRLF line endings and ASCII-only text. cmd.exe mis-parses
rem  LF-only batch files and UTF-8 blocks.
rem ===========================================================================
setlocal EnableExtensions EnableDelayedExpansion
title Genome Gene Explorer - build EXE
cd /d "%~dp0"
if errorlevel 1 goto :no_dir

call "%~dp0_find_python.bat"
if not defined GGE_PY goto :need_env

set "PY=!GGE_PY!"
echo ============================================
echo  Genome Gene Explorer - build EXE
echo ============================================
echo Python: !PY!
echo.

echo [1/4] Checking PyInstaller ...
"!PY!" -m PyInstaller --version >nul 2>nul
if errorlevel 1 (
    echo       not installed - installing requirements-dev.txt ...
    "!PY!" -m pip install -r requirements-dev.txt
    if errorlevel 1 goto :pip_failed
)

echo [2/4] Running the test suite ...
"!PY!" -m pytest -q
if errorlevel 1 (
    echo [WARN] Some tests failed - the build continues anyway.
)

echo.
echo [3/4] Building dist\GenomeGeneExplorer\ ...
"!PY!" -m PyInstaller --noconfirm --clean GenomeGeneExplorer.spec
if errorlevel 1 goto :build_failed

echo.
echo [4/4] Smoke test on the packed executable ...
"!PY!" "%~dp0tools\exe_smoke_test.py"
if errorlevel 1 (
    echo [ERROR] The packed executable failed its self-test.
    pause
    exit /b 1
)

echo.
echo Done: dist\GenomeGeneExplorer\GenomeGeneExplorer.exe
echo No Python needed on the target machine - copy the whole folder.
echo.
pause
exit /b 0


:no_dir
echo [ERROR] Cannot switch to the program directory: %~dp0
pause
exit /b 1

:need_env
echo [ERROR] No usable Python environment found.
echo   Run start.bat once to create .venv, or see tools\doctor.py for details.
pause
exit /b 1

:pip_failed
echo [ERROR] Could not install requirements-dev.txt.
echo   Try: "!PY!" -m pip install -r requirements-dev.txt -i https://pypi.tuna.tsinghua.edu.cn/simple
pause
exit /b 1

:build_failed
echo [ERROR] PyInstaller failed - see the output above.
pause
exit /b 1
