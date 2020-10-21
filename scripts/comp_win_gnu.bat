@echo off
echo.
:: Create lib folder
mkdir lib
mkdir lapack\lib
:: Build LAPACK routines
for %%f in (lapack/lapack_routine/*.f) do gfortran -c lapack/lapack_routine/%%f -o lapack/lib/%%~nf.o
:: Build FEANBEAM modules
for %%f in (src/*.f90) do gfortran -c src/%%f -Jlib -o lib/%%~nf.o

echo.
echo.
echo Finished building libraries.
echo.
