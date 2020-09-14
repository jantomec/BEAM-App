program cantilever
	
	use solver
	use mesher
	use mesh_object
	
	implicit none
	
	integer, parameter :: Nele = 5
	integer, parameter :: order = 1
	integer, parameter :: Ngauss = 1
	integer, dimension (1), parameter :: Nsteps = (/ 1 /)
	double precision, parameter :: L = 1.0D1
	double precision, parameter :: PI = 4 * atan (1.0D1), Q0 = 8 * PI
	integer, parameter :: MAXITER = 50
	double precision, parameter :: TOLER = 1D-10
	double precision, dimension (6), parameter :: material = (/ 1.0D1, 1.0D1, 1.0D1, 1.0D1, 2.0D1, 1.0D1 /)
	character (len = *), parameter :: folder = 'cantilever'
	character (len = *), parameter :: fname_format = '("step", I0.3, ".dat")'
	
	type (ElementMesh) :: mesh
	integer, parameter :: Nno = Nele * order + 1
	double precision, dimension (3, Nno) :: X, U
	double precision, dimension (6, 6) :: C
	logical, dimension (6, Nno) :: DOF
	double precision, dimension (6, Nno) :: dU, Q
	double precision, dimension (Nele, 6, Ngauss) :: p, f
	double precision, dimension (Nele, 3, Ngauss) :: om
	double precision, dimension (Nele, Ngauss, 3, 3) :: rot
	double precision, dimension (6, Nno) :: R
	integer :: i, j, st, Niter, info
	character (len = 255) :: fname_long
	character (len = :), allocatable :: fname
	character (len = :), allocatable :: relpath
	character (len = 255) :: wdir_long
	character (len = :), allocatable :: wdir
	character (len = 255) :: command
	logical :: did_dir
	
	! =================================================
	! MESH
	mesh%Nno = Nno
	mesh%Nele = Nele
	mesh%order = order
	call lmsh (L, mesh)
	do i = 1, Nele
		do j = 1, Ngauss
			rot (i, j, 1, :)= (/ 1.0D1, 0.0D1, 0.0D1 /)
			rot (i, j, 2, :)= (/ 0.0D1, 1.0D1, 0.0D1 /)
			rot (i, j, 3, :)= (/ 0.0D1, 0.0D1, 1.0D1 /)
		end do
	end do
	
	! =================================================
	! ELASTIC MODULI MATRIX
	C = 0.0D1  ! add material
	do i = 1, 6
		C (i, i) = material (i)
	end do
	
	! =================================================
	! DATA INITIALIZATION
	U = 0.0D1
	p = 0.0D1
	f = 0.0D1
	om = 0.0D1
	
	st = getcwd (wdir_long)
	wdir = trim (wdir_long)
	
	write (command, '("if exist", X, A, "\", A, X, "rmdir", X, "/S", X, "/Q", X, A, "\" A)') wdir, folder, wdir, folder
	call execute_command_line (command, .TRUE.)
	
	write (command, '("mkdir", X, A, "\" A)') wdir, folder
	call execute_command_line (command, .TRUE.)
	
	! =================================================
	! BOUNDARY CONDITIONS
	DOF = .TRUE.
	DOF (:, 1) = .FALSE.
	dU = 0.0D1
	Q = 0.0D1
	
	! =================================================
	! FORCE CONTROL ROUTINE
	if (Nsteps (1) .ne. 0) then
		do j = 0, Nsteps (1)
			if (j > 0) then
				write (6, '(/, "Step", X, I3)') j
				Q (5, Nno) = Q (5, Nno) + Q0 / nsteps (1)
				call newton_iter (mesh%ele, mesh%X0, U, C, DOF, dU, Q, p, rot, om, f, R, TOLER, MAXITER, 'RSD', Niter, info, .TRUE.)
			end if
			
			X (1, :) = mesh%X0 (1, :) + U (1, :)
			X (2, :) = mesh%X0 (2, :) + U (2, :)
			X (3, :) = mesh%X0 (3, :) + U (3, :)
			
			write (fname_long, fname_format) j
			fname = trim (fname_long)
			relpath = folder//'/'//fname
			open (unit = 12, file = relpath, status = 'unknown', action = 'write')
			write (12, '("node", 19X, "X1", 19X, "X2", 19X, "X3", 19X, "U1", 19X, "U2", 19X, "U3", 19X, "R1", 19X, "R2", 19X, &
				"R3", 19X, "M1", 19X, "M2", 19X, "M3")')
			do i = 1, Nno
				write (12, '(I0.4, X, ES20.11E3, X, ES20.11E3, X, ES20.11E3, X, ES20.11E3, X, ES20.11E3, X, ES20.11E3, X, ES20.11E3, X, &
					ES20.11E3, X, ES20.11E3, X, ES20.11E3, X, ES20.11E3, X, ES20.11E3)') i, X (1, i), X (2, i), X (3, i), U (1, i), &
					U (2, i), U (3, i), R (1, i), R (2, i), R (3, i), R (4, i), R (5, i), R (6, i)
			end do
			close(12)
			
		end do
	end if
	
end program cantilever