@echo off
:: Create lib folder
mkdir lib
mkdir lapack\lib
:: Build LAPACK routines
for %%f in (lapack/lapack_routine/*.f) do gfortran -c lapack/lapack_routine/%%f -o lapack/lib/%%~nf.o
:: Build FEANBEAM modules
for %%f in (source/*.f90) do gfortran -c source/%%f -Jlib -o lib/%%~nf.o