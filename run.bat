@echo off
setlocal EnableExtensions

cd /d "%~dp0"

set "CXX=g++"
set "CXXFLAGS=-std=c++17 -O3 -fopenmp -Wall -Wno-sign-compare"
set "PY="
set "CONDA_ENV=graphic"
set "CONDA_PY=D:\coding\anaconda\envs\graphic\python.exe"
set "TRI_SRC=gen.cpp"
set "TRI_EXE=gen.exe"
set "TRI_ARGS=--radius 1.0 --points 820 --min-distance 0.058 --seed 3 --output ball_mesh.txt"
set "IMPLICIT_VIEWER=implicit_euler.py"
set "EXPERIMENT_SCRIPT=experiment.py"
set "EXPLICIT_VIEWER=explicit_euler.py"
set "RC=0"

if /I "%~1"=="clean" goto clean
if /I "%~1"=="tri" goto run_tri
if /I "%~1"=="tet" goto run_tri
if /I "%~1"=="gen" goto run_tri
if /I "%~1"=="fem" goto run_fem
if /I "%~1"=="implicit" goto run_fem
if /I "%~1"=="experiment" goto run_experiment
if /I "%~1"=="femtest" goto run_experiment
if /I "%~1"=="explicit" goto run_explicit
if /I "%~1"=="viewer" goto run_fem
if /I "%~1"=="all" goto run_all
if /I "%~1"=="help" goto help
if /I "%~1"=="/?" goto help
if not "%~1"=="" goto unknown
goto run_all

:clean
call :clean_pycache
set "RC=%ERRORLEVEL%"
goto end

:run_all
call :clean_pycache
call :build_and_run_tri
if errorlevel 1 (
    set "RC=%ERRORLEVEL%"
    goto end
)
call :compile_and_run_python "%IMPLICIT_VIEWER%"
set "RC=%ERRORLEVEL%"
goto end

:run_tri
call :build_and_run_tri
set "RC=%ERRORLEVEL%"
goto end

:run_fem
call :compile_and_run_python "%IMPLICIT_VIEWER%"
set "RC=%ERRORLEVEL%"
goto end

:run_experiment
call :compile_and_run_python "%EXPERIMENT_SCRIPT%"
set "RC=%ERRORLEVEL%"
goto end

:run_explicit
call :compile_and_run_python "%EXPLICIT_VIEWER%"
set "RC=%ERRORLEVEL%"
goto end

:clean_pycache
echo.
echo [clean] Removing Python cache files...
if exist "__pycache__" rmdir /s /q "__pycache__"
for /d /r %%D in (__pycache__) do (
    if exist "%%D" rmdir /s /q "%%D"
)
for /r %%F in (*.pyc *.pyo) do (
    if exist "%%F" del /f /q "%%F"
)
echo [clean] Removing g++ build outputs...
if exist "%TRI_EXE%" del /f /q "%TRI_EXE%"
if exist "a.exe" del /f /q "a.exe"
for %%F in (*.o *.obj *.a *.lib *.dll *.pdb *.ilk *.exp) do (
    if exist "%%F" del /f /q "%%F"
)
echo [clean] Removing generated mesh and UI state...
if exist "ball_mesh.txt" del /f /q "ball_mesh.txt"
if exist "imgui.ini" del /f /q "imgui.ini"
echo [clean] Done.
exit /b 0

:build_and_run_tri
echo.
echo [tri] Building %TRI_SRC%...
where %CXX% >nul 2>nul
if errorlevel 1 (
    echo [tri] Cannot find %CXX%. Please install MinGW-w64 or add g++ to PATH.
    exit /b 1
)
if not exist "%TRI_SRC%" (
    echo [tri] Cannot find %TRI_SRC%.
    exit /b 1
)
%CXX% %CXXFLAGS% "%TRI_SRC%" -o "%TRI_EXE%"
if errorlevel 1 (
    echo [tri] Build failed. If %TRI_SRC% is only a library, add a main^(^) before running it as an executable.
    exit /b 1
)
echo [tri] Running %TRI_EXE%...
"%TRI_EXE%" %TRI_ARGS%
exit /b %ERRORLEVEL%

:compile_and_run_python
set "PY_SCRIPT=%~1"
echo.
echo [fem] Compile-checking %PY_SCRIPT%...
call :find_python
if errorlevel 1 exit /b 1
if not exist "%PY_SCRIPT%" (
    echo [fem] Cannot find %PY_SCRIPT%.
    exit /b 1
)
%PY% -m py_compile "%PY_SCRIPT%"
if errorlevel 1 (
    echo [fem] Python compile check failed.
    exit /b 1
)
echo [fem] Running %PY_SCRIPT%...
%PY% "%PY_SCRIPT%"
exit /b %ERRORLEVEL%

:find_python
if exist "%CONDA_PY%" (
    set "PY="%CONDA_PY%""
    echo [fem] Using conda env %CONDA_ENV%: %CONDA_PY%
    exit /b 0
)
where conda >nul 2>nul
if not errorlevel 1 (
    set "PY=conda run -n %CONDA_ENV% python"
    echo [fem] Using conda env %CONDA_ENV% through conda run.
    exit /b 0
)
where python >nul 2>nul
if not errorlevel 1 (
    set "PY=python"
    echo [fem] Conda env %CONDA_ENV% not found; using python from PATH.
    exit /b 0
)
where py >nul 2>nul
if not errorlevel 1 (
    set "PY=py -3"
    echo [fem] Conda env %CONDA_ENV% not found; using py -3.
    exit /b 0
)
echo [fem] Cannot find conda env %CONDA_ENV% or Python. Please check %CONDA_PY%.
exit /b 1

:unknown
echo Unknown command: %~1
goto help

:help
echo.
echo Usage:
echo   run.bat           Clean, then run triangulation and implicit_euler.py
echo   run.bat clean     Remove __pycache__ and .pyc/.pyo files
echo   run.bat tri       Build and run gen.cpp
echo   run.bat fem       Compile-check and run implicit_euler.py
echo   run.bat experiment Compile-check and run experiment.py
echo   run.bat femtest    Alias for experiment
echo   run.bat explicit  Compile-check and run explicit_euler.py
echo   run.bat all       Clean, then run triangulation and implicit_euler.py
echo.
goto end

:end
endlocal & exit /b %RC%
