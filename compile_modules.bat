gfortran -c source/vector_algebra.f90 -Jlib -o lib/vector_algebra.o
gfortran -c source/legendre_gauss.f90 -Jlib -o lib/vector_algebra.o
gfortran -c source/shape_functions.f90 -Jlib -o lib/vector_algebra.o
gfortran -c source/beam.f90 -Jlib -o lib/vector_algebra.o
gfortran -c source/solver.f90 -Jlib -o lib/vector_algebra.o
gfortran -c source/mesher.f90 -Jlib -o lib/vector_algebra.o