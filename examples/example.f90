program example

	use solver
	
	implicit none
	
	integer, dimension (1, 2) :: ele
	double precision, dimension (3, 2) :: X0, U
	double precision, dimension (6, 6) :: C
	logical, dimension (6, 2) :: DOF
	double precision, dimension (6, 2) :: dU, Q, res
	double precision, dimension (1, 6, 1) :: p, f
	double precision, dimension (1, 3, 1) :: om
	double precision, dimension (1, 1, 3, 3) :: rot
	integer, parameter :: MAXITER = 20
	double precision, parameter :: TOLER = 1D-5
	integer :: i, Niter, errck
		
	ele (1, :) = (/ 1, 2 /)
	X0 (:, 1) = (/ 0.0D1, 0.0D1, 0.0D1 /)
	X0 (:, 2) = (/ 1.0D1, 2.0D1, 3.0D1 /)
	U = 0.0D1
	C = 0.0D1
	do i = 1, 6
		C (i, i) = 1.0D1
	end do
	DOF (:, 1) = .FALSE.
	DOF (:, 2) = .TRUE.
	dU = 0.0D1
	Q = 0.0D1
	Q (2, 2) = 0.001
	p = 0.0D1
	f = 0.0D1
	om = 0.0D1
	rot = 0.0D1
	do i = 1, 3
		rot (1, 1, i, i) = 1
	end do
	
	call newton_iter (ele, X0, U, C, DOF, dU, Q, p, rot, om, f, res, TOLER, MAXITER, 'RSD', Niter, errck, .TRUE.)
	
	print *, U (:, 2)
	print *, rot (1, 1, :, :)
	
end program example
