@echo off
rem ===========================================================================
rem  Locate a usable Python interpreter for Genome Gene Explorer.
rem
rem  This file deliberately has NO setlocal: the caller must see the results.
rem
rem  Output variables (path to python.exe, WITHOUT surrounding quotes -
rem  the caller is responsible for quoting it):
rem    GGE_PY        interpreter that already has PySide6 + Bio + pyfaidx
rem    GGE_PY_PIP    any Python >= 3.10 that has a working pip
rem    GGE_PY_BASE   any Python >= 3.10 (pip not required)
rem
rem  Override:  set GGE_PYTHON=C:\path\to\python.exe
rem ===========================================================================
set "GGE_PY="
set "GGE_PY_PIP="
set "GGE_PY_BASE="

if defined GGE_PYTHON call :probe "!GGE_PYTHON!"

if not defined GGE_PY if exist ".venv\Scripts\python.exe" call :probe ".venv\Scripts\python.exe"

rem -- the official py launcher knows every registered Python ----------------
if not defined GGE_PY (
    for %%V in (3.13 3.12 3.11 3.10) do (
        if not defined GGE_PY call :probe_launcher %%V
    )
)

rem -- whatever "python" resolves to on PATH -------------------------------
if not defined GGE_PY call :probe_where_python

rem -- common installation directories -------------------------------------
if not defined GGE_PY call :probe_dirs

rem -- give up: report what we have ----------------------------------------
if defined GGE_PY exit /b 0
if defined GGE_PY_PIP set "GGE_PY_BASE=%GGE_PY_PIP%"
exit /b 0


rem ===========================================================================
:probe
rem %1 = path to python.exe (quoted when the path contains spaces)
set "CAND=%~1"
if "%CAND%"=="" goto :eof
if not exist "%CAND%" goto :eof
rem Skip the Windows Store stub: launching it pops up the Store instead of Python.
if not "%CAND:WindowsApps=%"=="%CAND%" goto :eof

"%CAND%" -c "import sys;raise SystemExit(0 if sys.version_info[:2]>=(3,10) else 1)" >nul 2>nul
if errorlevel 1 goto :eof
if not defined GGE_PY_BASE set "GGE_PY_BASE=%CAND%"

"%CAND%" -m pip --version >nul 2>nul
if not errorlevel 1 if not defined GGE_PY_PIP set "GGE_PY_PIP=%CAND%"

"%CAND%" -c "import PySide6, Bio, pyfaidx" >nul 2>nul
if errorlevel 1 goto :eof
set "GGE_PY=%CAND%"
goto :eof


:probe_launcher
rem %1 = version such as 3.13
for /f "delims=" %%P in ('py -%1 -c "import sys;print(sys.executable)" 2^>nul') do call :probe "%%P"
goto :eof


:probe_where_python
for /f "delims=" %%P in ('where python 2^>nul') do (
    if not defined GGE_PY call :probe "%%P"
)
goto :eof


:probe_dirs
for /d %%D in ("%LOCALAPPDATA%\Programs\Python\Python3*") do (
    if not defined GGE_PY if exist "%%~D\python.exe" call :probe "%%~D\python.exe"
)
for /d %%D in ("%ProgramFiles%\Python3*") do (
    if not defined GGE_PY if exist "%%~D\python.exe" call :probe "%%~D\python.exe"
)
for /d %%D in ("C:\Python3*") do (
    if not defined GGE_PY if exist "%%~D\python.exe" call :probe "%%~D\python.exe"
)
rem Interpreters that ship alongside a tool (no PATH entry, no registry entry).
for /d %%D in ("%USERPROFILE%\.workbuddy\binaries\python\versions\*") do (
    if not defined GGE_PY if exist "%%~D\python.exe" call :probe "%%~D\python.exe"
)
goto :eof
