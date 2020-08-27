mkdir lib
gfortran -c source/vector_algebra.f90 -Jlib -o lib/vector_algebra.o
gfortran -c source/legendre_gauss.f90 -Jlib -o lib/legendre_gauss.o
gfortran -c source/shape_functions.f90 -Jlib -o lib/shape_functions.o
gfortran -c source/beam.f90 -Jlib -o lib/beam.o
gfortran -c source/solver.f90 -Jlib -o lib/solver.o
gfortran -c source/mesher.f90 -Jlib -o lib/mesher.o