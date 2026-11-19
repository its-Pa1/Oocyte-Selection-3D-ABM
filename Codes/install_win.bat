@echo off
:: ============================================================
::  install_simuscale_env.bat
::  One-stop installer:
::    - Detects if Visual Studio 2022 is already installed
::    - If NOT: downloads and installs VS 2022 Community silently
::    - If YES: skips VS install entirely
::    - Sets up vcpkg (VS built-in or standalone fallback)
::    - Installs all C++ / OpenGL libraries via vcpkg
:: ============================================================

setlocal enabledelayedexpansion

echo ============================================================
echo  SimuScale Environment Installer
echo  Visual Studio 2022 + vcpkg + OpenGL dependencies
echo ============================================================
echo.

:: ------------------------------------------------------------------
:: 1. Check for admin rights
:: ------------------------------------------------------------------
net session >nul 2>&1
if %errorLevel% neq 0 (
    echo [ERROR] Please run this script as Administrator.
    pause
    exit /b 1
)

:: ------------------------------------------------------------------
:: 2. Detect if Visual Studio 2022 is already installed
:: ------------------------------------------------------------------
echo [INFO] Checking for existing Visual Studio 2022 installation...
set "VS_INSTALLED=0"
set "VS_EDITION="

for %%E in (Community Professional Enterprise) do (
    if exist "C:\Program Files\Microsoft Visual Studio\2022\%%E\Common7\IDE\devenv.exe" (
        set "VS_INSTALLED=1"
        set "VS_EDITION=%%E"
    )
)

:: Also check BuildTools
if exist "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\Tools\VsDevCmd.bat" (
    set "VS_INSTALLED=1"
    set "VS_EDITION=BuildTools"
)

if "%VS_INSTALLED%"=="1" (
    echo [OK] Visual Studio 2022 %VS_EDITION% already installed — skipping download.
    goto :SETUP_VCPKG
)

:: ------------------------------------------------------------------
:: 3. VS NOT found — download and install VS 2022 Community
:: ------------------------------------------------------------------
echo [INFO] Visual Studio 2022 not found. Starting installation...
echo.

set "DLDIR=C:\VS_Install_Temp"
if not exist "%DLDIR%" mkdir "%DLDIR%"

set "VS_BOOTSTRAPPER=%DLDIR%\vs_Community.exe"
if not exist "%VS_BOOTSTRAPPER%" (
    echo [INFO] Downloading Visual Studio 2022 Community bootstrapper...
    powershell -NoProfile -Command ^
        "Invoke-WebRequest -Uri 'https://aka.ms/vs/17/release/vs_Community.exe' -OutFile '%VS_BOOTSTRAPPER%' -UseBasicParsing"
    if not exist "%VS_BOOTSTRAPPER%" (
        echo [ERROR] Download failed. Check your internet connection.
        pause
        exit /b 1
    )
    echo [OK] Bootstrapper downloaded.
) else (
    echo [OK] Bootstrapper already cached, skipping download.
)

echo.
echo [INFO] Installing Visual Studio 2022 Community (15-30 min)...
"%VS_BOOTSTRAPPER%" --wait --quiet --norestart ^
    --add Microsoft.VisualStudio.Workload.NativeDesktop ^
    --add Microsoft.VisualStudio.Workload.NativeCrossPlat ^
    --add Microsoft.VisualStudio.Component.VC.CMake.Project ^
    --add Microsoft.VisualStudio.Component.Git ^
    --add Microsoft.VisualStudio.Component.Vcpkg ^
    --includeRecommended

if %errorLevel% neq 0 (
    if %errorLevel% neq 3010 (
        echo [ERROR] Visual Studio installation failed (code: %errorLevel%)
        pause
        exit /b %errorLevel%
    )
)
set "VS_EDITION=Community"
echo [OK] Visual Studio 2022 Community installed.

:: Clean up bootstrapper temp folder
set /p CLEANUP="[?] Delete installer temp files in %DLDIR%? (y/n): "
if /i "%CLEANUP%"=="y" rd /s /q "%DLDIR%"

:: ------------------------------------------------------------------
:: 4. SETUP_VCPKG — locate or install vcpkg
:: ------------------------------------------------------------------
:SETUP_VCPKG
echo.
echo [INFO] Locating vcpkg...

set "VCPKG_EXE="

:: Try VS built-in vcpkg first (all editions)
for %%E in (Community Professional Enterprise) do (
    set "_P=C:\Program Files\Microsoft Visual Studio\2022\%%E\VC\vcpkg\vcpkg.exe"
    if exist "!_P!" set "VCPKG_EXE=!_P!"
)

:: Try BuildTools vcpkg
set "_P=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\vcpkg\vcpkg.exe"
if exist "%_P%" set "VCPKG_EXE=%_P%"

:: Try standalone vcpkg
if "%VCPKG_EXE%"=="" (
    if exist "C:\vcpkg\vcpkg.exe" set "VCPKG_EXE=C:\vcpkg\vcpkg.exe"
)

:: None found — clone and bootstrap standalone vcpkg
if "%VCPKG_EXE%"=="" (
    echo [INFO] vcpkg not found. Cloning standalone vcpkg to C:\vcpkg...

    where git >nul 2>&1
    if %errorLevel% neq 0 (
        :: Add VS-bundled git to PATH
        for %%E in (Community Professional Enterprise BuildTools) do (
            set "_GP=C:\Program Files\Microsoft Visual Studio\2022\%%E\Common7\IDE\CommonExtensions\Microsoft\TeamFoundation\Team Explorer\Git\cmd"
            if exist "!_GP!\git.exe" set "PATH=!_GP!;%PATH%"
        )
        if exist "C:\Program Files\Git\cmd\git.exe" set "PATH=C:\Program Files\Git\cmd;%PATH%"
    )

    git clone https://github.com/microsoft/vcpkg.git C:\vcpkg
    call C:\vcpkg\bootstrap-vcpkg.bat -disableMetrics
    set "VCPKG_EXE=C:\vcpkg\vcpkg.exe"
)

:: Derive VCPKG_ROOT from exe path
for %%F in ("%VCPKG_EXE%") do set "VCPKG_ROOT=%%~dpF"
if "%VCPKG_ROOT:~-1%"=="\" set "VCPKG_ROOT=%VCPKG_ROOT:~0,-1%"

echo [OK] vcpkg found at: %VCPKG_ROOT%
setx VCPKG_ROOT "%VCPKG_ROOT%" /M
set "PATH=%VCPKG_ROOT%;%PATH%"

:: ------------------------------------------------------------------
:: 5. Integrate vcpkg with Visual Studio / MSBuild
:: ------------------------------------------------------------------
echo.
echo [INFO] Integrating vcpkg with Visual Studio...
"%VCPKG_EXE%" integrate install

:: ------------------------------------------------------------------
:: 6. Install C++ / OpenGL libraries
:: ------------------------------------------------------------------
echo.
echo [INFO] Installing OpenGL...
"%VCPKG_EXE%" install opengl:x64-windows

echo.
echo [INFO] Installing GLFW3 (windowing + input)...
"%VCPKG_EXE%" install glfw3:x64-windows

echo.
echo [INFO] Installing GLEW (OpenGL extension wrangler)...
"%VCPKG_EXE%" install glew:x64-windows

echo.
echo [INFO] Installing GLM (math / linear algebra)...
"%VCPKG_EXE%" install glm:x64-windows

echo.
echo [INFO] Installing freeglut...
"%VCPKG_EXE%" install freeglut:x64-windows

:: ------------------------------------------------------------------
:: 7. Ensure CMake is on PATH
:: ------------------------------------------------------------------
echo.
where cmake >nul 2>&1
if %errorLevel% neq 0 (
    for %%E in (Community Professional Enterprise BuildTools) do (
        set "_CP=C:\Program Files\Microsoft Visual Studio\2022\%%E\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin"
        if exist "!_CP!\cmake.exe" (
            setx PATH "%PATH%;!_CP!" /M
            set "PATH=!_CP!;%PATH%"
            echo [OK] CMake added to PATH from VS %%E
        )
    )
) else (
    echo [OK] CMake already on PATH:
    cmake --version
)

:: ------------------------------------------------------------------
:: 8. Summary
:: ------------------------------------------------------------------
echo.
echo ============================================================
echo  Installation complete!
echo ============================================================
echo  Visual Studio 2022 %VS_EDITION%   [ready]
echo  vcpkg at: %VCPKG_ROOT%            [integrated]
echo  OpenGL / GLFW3 / GLEW / GLM       [installed]
echo  freeglut                           [installed]
echo ============================================================
echo.
echo  Next step: run  build_simuscale_windows.bat
echo.
echo  Or build manually:
echo    mkdir build ^&^& cd build
echo    cmake .. -DCMAKE_TOOLCHAIN_FILE="%VCPKG_ROOT%\scripts\buildsystems\vcpkg.cmake"
echo    cmake --build . --config Release
echo ============================================================
echo.
if "%errorLevel%"=="3010" (
    echo [WARN] A system reboot is recommended before building.
)
pause
