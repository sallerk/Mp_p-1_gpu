@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0"

rem --- locate Visual Studio's x64 build environment --------------------------
set "VCVARS="
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
rem  note: !VSWHERE! (delayed) not %VSWHERE% -- the path contains "(x86)" and a
rem  literal ")" inside a parenthesised block would terminate it early.
if exist "!VSWHERE!" (
    "!VSWHERE!" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath > "%TEMP%\_vspath.txt" 2>nul
    set /p VSPATH=<"%TEMP%\_vspath.txt"
    del /q "%TEMP%\_vspath.txt" 2>nul
    if exist "!VSPATH!\VC\Auxiliary\Build\vcvars64.bat" set "VCVARS=!VSPATH!\VC\Auxiliary\Build\vcvars64.bat"
)
if not defined VCVARS (
    for %%v in (2022 2019) do (
        for %%e in (Enterprise Professional Community BuildTools) do (
            if exist "%ProgramFiles%\Microsoft Visual Studio\%%v\%%e\VC\Auxiliary\Build\vcvars64.bat" (
                set "VCVARS=%ProgramFiles%\Microsoft Visual Studio\%%v\%%e\VC\Auxiliary\Build\vcvars64.bat"
            )
        )
    )
)
if not defined VCVARS (
    echo ERROR: could not find vcvars64.bat - install Visual Studio 2019/2022
    echo        with the "Desktop development with C++" workload.
    exit /b 1
)

rem  Captured before vcvars64.bat runs: it sets its own VCPKG_ROOT
rem  (pointing at the vcpkg bundled with Visual Studio itself), unconditionally
rem  overwriting anything already in the environment. An explicit override set
rem  before calling build.bat is kept here so it still wins.
set "USER_VCPKG_ROOT=%VCPKG_ROOT%"

echo Using %VCVARS%
call "%VCVARS%" >nul 2>nul
if errorlevel 1 ( echo ERROR: vcvars64 failed & exit /b 1 )

rem --- regenerate the kernel bundle from src/cl/*.cl --------------------------
rem  MSVC caps a string literal at 64 KB, so genbundle.py emits chunks.
where python >nul 2>nul
if not errorlevel 1 (
    python tools\genbundle.py || ( echo ERROR: genbundle.py failed & exit /b 1 )
) else (
    echo Skipping bundle regeneration ^(python not found^); using existing src\bundle.cpp
)

rem --- locate GMP (vcpkg, x64-windows-static triplet) -------------------------
rem  Added in 1.8: Gcd.cpp's gcdGmp delegates the production gcd(x-1,2^p-1)
rem  call to GMP's mpz_gcd -- ~11x faster than this project's own hand-rolled
rem  gcdHalf at production scale (see gcdGmp's own comment in Gcd.cpp for the
rem  numbers). Statically linked (x64-windows-static triplet), so there is no
rem  DLL to ship -- Mp_p-1_gpu.exe stays the single file every prior release
rem  was. Detected the same way vcvars64.bat is above: use what is already
rem  installed, fail with setup instructions if not -- this is not
rem  auto-installed, the same as Visual Studio itself is not. Three places
rem  are checked, in order: an explicit override (USER_VCPKG_ROOT, captured
rem  above), the vcpkg bundled with Visual Studio itself (VCVARS just set
rem  VCPKG_ROOT to it -- every machine that can already run this script has
rem  one, nothing extra to install to get vcpkg itself), and a vcpkg folder
rem  placed next to this project.
set "VCPKG_DIR="
if defined USER_VCPKG_ROOT if exist "%USER_VCPKG_ROOT%\installed\x64-windows-static\include\gmp.h" set "VCPKG_DIR=%USER_VCPKG_ROOT%"
if not defined VCPKG_DIR if defined VCPKG_ROOT if exist "%VCPKG_ROOT%\installed\x64-windows-static\include\gmp.h" set "VCPKG_DIR=%VCPKG_ROOT%"
if not defined VCPKG_DIR if exist "%~dp0vcpkg\installed\x64-windows-static\include\gmp.h" set "VCPKG_DIR=%~dp0vcpkg"
if not defined VCPKG_DIR (
    echo ERROR: GMP not found. Mp_p-1_gpu links GMP ^(x64-windows-static^) for its gcd step.
    echo        Install it once using the vcpkg that ships with Visual Studio:
    echo            "%VCPKG_ROOT%\vcpkg.exe" install gmp:x64-windows-static
    echo        ^(if that path doesn't exist, clone https://github.com/microsoft/vcpkg,
    echo        run bootstrap-vcpkg.bat, then vcpkg install gmp:x64-windows-static,
    echo        and either set VCPKG_ROOT to that folder or place it at
    echo        %~dp0vcpkg^)
    exit /b 1
)
set "GMP_INC=%VCPKG_DIR%\installed\x64-windows-static\include"
set "GMP_LIB=%VCPKG_DIR%\installed\x64-windows-static\lib\gmp.lib"

if not exist build mkdir build

rem --- compile ---------------------------------------------------------------
rem  No OpenCL SDK needed: src\clshim.cpp resolves the entry points from the
rem  driver's OpenCL.dll at run time.
set SRCS=^
 src\main.cpp ^
 src\clshim.cpp ^
 src\bundle.cpp ^
 src\common.cpp ^
 src\log.cpp ^
 src\File.cpp ^
 src\fs.cpp ^
 src\timeutil.cpp ^
 src\clwrap.cpp ^
 src\Queue.cpp ^
 src\Event.cpp ^
 src\TimeInfo.cpp ^
 src\Profile.cpp ^
 src\Kernel.cpp ^
 src\KernelCompiler.cpp ^
 src\AllocTrac.cpp ^
 src\Trig.cpp ^
 src\TrigBufCache.cpp ^
 src\FFTConfig.cpp ^
 src\state.cpp ^
 src\sha3.cpp ^
 src\gpuid.cpp ^
 src\TuneEntry.cpp ^
 src\CycleFile.cpp ^
 src\Args.cpp ^
 src\Gpu.cpp ^
 src\Proof.cpp ^
 src\md5.cpp ^
 src\Saver.cpp ^
 src\Signal.cpp ^
 src\tune.cpp ^
 src\Primes.cpp ^
 src\BigInt.cpp ^
 src\Gcd.cpp ^
 src\Parallel.cpp ^
 src\Pool.cpp ^
 src\testBigInt.cpp ^
 src\Bounds.cpp ^
 src\Config.cpp ^
 src\Worktodo.cpp ^
 src\PM1.cpp ^
 src\Selftest.cpp ^
 src\Stage2Plan.cpp ^
 src\Stage2Save.cpp ^
 src\Pp1Stage2Save.cpp ^
 src\Save.cpp ^
 src\Results.cpp

rem  No /Isrc: Windows include lookup is case-insensitive, so putting src\ on the
rem  include path makes the CRT's <csignal> -> <signal.h> resolve to our own
rem  src\Signal.h. Quoted includes already search the including file's directory.
cl /nologo /O2 /EHsc /std:c++20 /MT /W3 /D_CRT_SECURE_NO_WARNINGS /DNOMINMAX ^
   /I"%GMP_INC%" /Fe:Mp_p-1_gpu.exe /Fo:build\ %SRCS% "%GMP_LIB%"
if errorlevel 1 ( echo. & echo BUILD FAILED & exit /b 1 )

echo.
echo Built Mp_p-1_gpu.exe
