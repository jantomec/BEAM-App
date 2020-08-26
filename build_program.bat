@echo off
IF %1.==. GOTO No1

ifort /Qmkl:sequential /F500000000 -module:modules %1 modules/mesher.obj modules/solver.obj modules/beam.obj modules/shape_functions.obj modules/legendre_gauss.obj modules/vector_algebra.obj
GOTO End1

:No1
  ECHO No filename. Correct usage: build_program.bat filename.f90
GOTO End1

:End1