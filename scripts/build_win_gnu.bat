@echo off

echo.

IF %1.==. GOTO No1

setlocal enabledelayedexpansion enableextensions

set LIST1=
for %%f in (lib\*.o) do set LIST1=!LIST1! %%f
set LIST1=%LIST1:~1%

set LIST2=
for %%f in (lapack\lib\*.o) do set LIST2=!LIST2! %%f
set LIST2=%LIST2:~1%

gfortran -Ilib -Ilapack/lib %1 !LIST2! !LIST1! -o %~n1.exe

echo.
echo.
echo Finished building program %~n1.exe.
echo.

GOTO End1

:No1
  ECHO No filename. Correct usage: build_program.bat filename.f90
GOTO End1

:End1