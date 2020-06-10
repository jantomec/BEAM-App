program example

	use f95_precision
	use solver
	
	implicit none
	
	integer, dimension (1, 2) :: ele
	real (DP), dimension (3, 2) :: X0, U
	real (DP), dimension (6, 6) :: C
	logical, dimension (6, 2) :: DOF
	real (DP), dimension (6, 2) :: dU, Q
	real (DP), dimension (1, 6, 1) :: p, f
	real (DP), dimension (1, 3, 1) :: om
	real (DP), dimension (1, 1, 3, 3) :: rot
	integer, parameter :: MAXITER = 20
	real (DP), parameter :: TOLER = 1E-5
	integer :: i
		
	ele (1, :) = (/ 1, 2 /)
	X0 (:, 1) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
	X0 (:, 2) = (/ 1.0_DP, 2.0_DP, 3.0_DP /)
	U = 0.0_DP
	C = 0.0_DP
	do i = 1, 6
		C (i, i) = 1.0_DP
	end do
	DOF (:, 1) = .FALSE.
	DOF (:, 2) = .TRUE.
	dU = 0.0_DP
	Q = 0.0_DP
	Q (2, 2) = 0.001
	p = 0.0_DP
	f = 0.0_DP
	om = 0.0_DP
	rot = 0.0_DP
	do i = 1, 3
		rot (1, 1, i, i) = 1
	end do
	
	call newton_iter (ele, X0, U, C, DOF, dU, Q, p, rot, om, f, TOLER, MAXITER)
	
	print *, U (:, 2)
	print *, rot (1, 1, :, :)
	
end program example
