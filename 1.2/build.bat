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
 src\testBigInt.cpp ^
 src\Bounds.cpp ^
 src\Config.cpp ^
 src\PM1.cpp ^
 src\Selftest.cpp ^
 src\Stage2Plan.cpp ^
 src\Stage2Save.cpp ^
 src\Pp1Stage2Save.cpp ^
 src\Save.cpp

rem  No /Isrc: Windows include lookup is case-insensitive, so putting src\ on the
rem  include path makes the CRT's <csignal> -> <signal.h> resolve to our own
rem  src\Signal.h. Quoted includes already search the including file's directory.
cl /nologo /O2 /EHsc /std:c++20 /MT /W3 /D_CRT_SECURE_NO_WARNINGS /DNOMINMAX ^
   /Fe:Mp_p-1_gpu.exe /Fo:build\ %SRCS%
if errorlevel 1 ( echo. & echo BUILD FAILED & exit /b 1 )

echo.
echo Built Mp_p-1_gpu.exe
