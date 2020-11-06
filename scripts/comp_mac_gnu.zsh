#!/bin/zsh

# Create lib folder
echo "\nCreating directories"
mkdir lib
mkdir lapack/lib

# Build LAPACK routines
echo "\n============================\n"
echo "Building LAPACK libraries\n"
for file in lapack/lapack_routine/*.f(:t:r); do
  echo "Building" $file.f
  gfortran -c lapack/lapack_routine/$file.f -o lapack/lib/$file.o
done

# Build BEAM App modules
echo "\n============================\n"
echo "Building BEAM App libraries\n"
for file in source/*.f90(:t:r); do
  echo "Building" $file.f90
  gfortran -c source/$file.f90 -Jlib -o lib/$file.o
done
echo "\n============================\n"
echo "Finished compiling libraries\n"
