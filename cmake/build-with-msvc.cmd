@echo off
setlocal enabledelayedexpansion

if "%~2"=="" (
  echo Usage: build-with-msvc.cmd ^<configure-preset^> ^<build-dir^> [target]
  echo Example: build-with-msvc.cmd vcpkg-x64 build kinDS-demo
  exit /b 1
)

set "PRESET=%~1"
set "BUILDDIR=%~2"
set "TARGET=%~3"
if "%TARGET%"=="" set "TARGET=kinDS-demo"

set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%VSWHERE%" (
  echo Error: vswhere not found. Install Visual Studio 2022 with C++ tools.
  exit /b 1
)

for /f "usebackq tokens=*" %%i in (`"%VSWHERE%" -latest -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VSPATH=%%i"
if not defined VSPATH (
  echo Error: Visual Studio with MSVC x64 tools not found.
  exit /b 1
)

call "%VSPATH%\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 (
  echo Error: Failed to initialize MSVC environment.
  exit /b 1
)

pushd "%~dp0.."
cmake --preset "%PRESET%"
if errorlevel 1 (
  popd
  exit /b 1
)
cmake --build "%BUILDDIR%" --target "%TARGET%"
set "EXITCODE=%ERRORLEVEL%"
popd
exit /b %EXITCODE%
