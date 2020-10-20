! =============================================================================
!
! This is example 7.3 from Simo & Vu-Quoc [1985] - 
! 	A three-dimensional finite-strain rod model. 
! 	Part II: Computational aspects
!
! =============================================================================
! 
!	 			   | F						
!	    		   |						
!	        _ _ _ _v_ _ _ _					EA		= 170000
!	      /				    \				GAy		= 170
!	     /					 \				GAz		= 17000
!	    /					  \				Jx		= 170
!	   /					   \			EIyy	= 170
!	  /						    \			EIzz	= 170
!	 |							 |			N 		= 40
!	 |							 |			L 		= 1.5
!	 |							 |
!	 0                           -
!
! =============================================================================

program test03
	
	use solver
	use mesher
	use mesh_object
	use legendre_gauss
	
	implicit none
	
	integer, parameter :: Nele = 40
	integer, parameter :: order = 1
	integer, parameter :: Ngauss = 1
	double precision, parameter :: L = 1.5D0
	double precision, parameter :: PI = 4 * atan (1.0D0), DEG = PI / 180
	double precision, parameter :: Q0 = -1400.0D0
	double precision, parameter :: dS = 0.2D0
	integer, parameter :: MAXITER = 30
	double precision, parameter :: TOLER = 1D-8
	integer, parameter :: MAXSTEPS = 1000
	double precision, dimension (6), parameter :: material = (/ 170000.0D0, 170.0D0, 17000.0D0, 170.0D0, 170.0D0, 170.0D0 /)
	character (len = 12), parameter :: folder = 'arch-results'
	character (len = 25), parameter :: fname_format = '("arch_", I0.3, ".dat")'
	
	type (ElementMesh) :: mesh
	double precision :: lambda
	double precision, dimension (Ngauss) :: gpts, gwgts
	integer, parameter :: Nno = Nele * order + 1
	double precision, dimension (3, Nno) :: X, U
	double precision, dimension (3 * Nno) :: Uinc
	double precision, dimension (6, 6) :: C
	logical, dimension (6, Nno) :: DOF
	double precision, dimension (6, Nno) :: Q, QC
	double precision, dimension (Nele, 6, Ngauss) :: p, f
	double precision, dimension (Nele, 3, Ngauss) :: om
	double precision, dimension (Nele, Ngauss, 3, 3) :: rot
	double precision, dimension (6, Nno) :: R
	integer :: i, j, st, Niter, info
	character (len = 80) :: fname
	character (len = 80) :: fullname
	character (len = 80) :: wdir
	character (len = 160) :: command
	
	! =================================================
	! MESH
	mesh%Nno = Nno
	mesh%Nele = Nele
	mesh%order = order
	call legauss (Ngauss, gpts, gwgts)
	call arcmsh (L, -15.0D0 * DEG, 195.0D0 * DEG, mesh, rot, gpts)
	
	! =================================================
	! ELASTIC MODULI MATRIX
	C = 0.0D0  ! add material
	do i = 1, 6
		C (i, i) = material (i)
	end do
	
	! =================================================
	! DATA INITIALIZATION
	lambda = 0.0D0
	U = 0.0D0
	Uinc = 0.0D0
	p = 0.0D0
	f = 0.0D0
	om = 0.0D0
	
	st = getcwd(wdir)
	! Create directory for results
	write (command, '("MKDIR", X, A, "\" A)') trim(wdir), trim(folder)
	call execute_command_line (command, .TRUE.)
	
	! =================================================
	! BOUNDARY CONDITIONS
	DOF = .TRUE.
	DOF (:, 1) = .FALSE.
	DOF (:, Nno) = (/ .FALSE., .FALSE., .FALSE., .TRUE., .TRUE., .TRUE.  /)
	QC = 0.0D0
	Q = 0.0D0
	Q (1, Nno / 2 + 1) = Q0
	
	! =================================================
	! ARC-LENGTH ROUTINE
	do j = 0, MAXSTEPS
	
		if (j > 0) then
			write (6, '(/, "Step", X, I3)') j
			call arc_length_iter (mesh%ele, mesh%X0, Uinc, U, C, DOF, Q, QC, p, rot, om, f, R, lambda, dS, TOLER, MAXITER, Niter, info)
		end if
		
		X (1, :) = mesh%X0 (1, :) + U (1, :)
		X (2, :) = mesh%X0 (2, :) + U (2, :)
		X (3, :) = mesh%X0 (3, :) + U (3, :)
		
		write (fname, fname_format) j
		fullname = trim (folder)//'/'//trim (fname)
		open (unit = 12, file = trim (fullname), status = 'unknown', action = 'write')
		write (12, '("node", 19X, "X1", 19X, "X2", 19X, "X3", 19X, "U1", 19X, "U2", 19X, "U3", 19X, "R1", 19X, &
			"R2", 19X, "R3", 19X, "M1", 19X, "M2", 19X, "M3")')
		do i = 1, Nno
			write (12, '(I0.4, X, ES20.13, X, ES20.13, X, ES20.13, X, ES20.13, X, ES20.13, X, &
			ES20.13, X, ES20.13, X, ES20.13, X, ES20.13, X, ES20.13, X, ES20.13, X, &
			ES20.13)') i, X (1, i), X (2, i), X (3, i), U (1, i), U (2, i), U (3, i), R (1, i), &
			R (2, i), R (3, i), R (4, i), R (5, i), R (6, i)
		end do
		close(12)
		
		if (abs (lambda) > 1) exit
		
	end do
	
	if (j .eq. MAXSTEPS + 1) write (6, '(/, "Solution not converging fast enough")')
	
end program test03