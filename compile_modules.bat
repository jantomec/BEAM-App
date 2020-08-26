ifort /Qmkl:sequential -c -o modules/vector_algebra /module:modules core/vector_algebra.f90
ifort /Qmkl:sequential -c -o modules/legendre_gauss /module:modules core/legendre_gauss.f90
ifort /Qmkl:sequential -c -o modules/shape_functions /module:modules core/shape_functions.f90
ifort /Qmkl:sequential -c -o modules/beam /module:modules core/beam.f90
ifort /Qmkl:sequential -c -o modules/solver /I:eispack /module:modules core/solver.f90
ifort /Qmkl:sequential -c -o modules/mesher /module:modules core/mesher.f90