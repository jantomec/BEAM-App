gfortran -ffree-fortran -c modules/vector_algebra core/vector_algebra.f90
gfortran -ffree-fortran -c modules/legendre_gauss core/legendre_gauss.f90
gfortran -ffree-fortran -c modules/shape_functions core/shape_functions.f90
gfortran -ffree-fortran -c modules/beam core/beam.f90
::gfortran -ffree-fortran -c modules/solver /I:eispack core/solver.f90
::gfortran -ffree-fortran -c modules/mesher core/mesher.f90