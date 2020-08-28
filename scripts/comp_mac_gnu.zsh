#!/bin/zsh

# Create lib folder
mkdir lib
mkdir lapack/lib

# Build LAPACK routines
for file in lapack/lapack_routine/*.f(:t:r); do
  gfortran -c lapack/lapack_routine/$file.f -o lapack/lib/$file.o
done

# Build FEANBEAM modules
for file in source/*.f90(:t:r); do
  gfortran -c source/$file.f90 -Jlib -o lib/$file.o
done

echo "Finished compiling libraries"