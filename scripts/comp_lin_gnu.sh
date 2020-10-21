#!/bin/sh

# Create lib folder
mkdir lib
mkdir lapack/lib

# Build LAPACK routines
for file in lapack/lapack_routine/*.f
do
  fname="${file##*/}"
  name="${fname%.f}"
  gfortran -c lapack/lapack_routine/$name.f -o lapack/lib/$name.o
done

# Build FEANBEAM modules
for file in src/*.f90
do
  fname="${file##*/}"
  name="${fname%.f90}"
  gfortran -c src/$name.f90 -Jlib -o lib/$name.o
done

echo "Finished compiling libraries"
