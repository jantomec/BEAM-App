! Finite element solver
!
! author ....... Jan Tomec
! copyright .... Copyright 2020, Project THREAD - University of Rijeka, Faculty of Civil Engineering
! credits ...... Jan Tomec, Gordan Jelenić
! license ...... GPL
! version ...... 1.0.0
! maintainer ... Jan Tomec
! email ........ jan.tomec@gradri.uniri.hr
! status ....... Development
! date ......... 06/09/2020
!
! ------------------------------------------------------------------------------

module solver

	use lapack95
	use f95_precision
	use beam
	
	implicit none
	
	private
	
	public :: newton_iter
	
	contains
	
	! perform newton-raphson iteration
	!
	! ele .......... elements array (no ele, no nodes on ele)
	! X0 ........... initial coordinates (3, no all nodes)
	! U ............ displacement (3, no all nodes)
	! C ............ tangent elastic moduli (6, 6)
	! DOF .......... degrees of freedom (6, no all nodes)
	! dU ........... prescribed displacement/rotation (6, no all nodes)
	! Q ............ nodal load (6, no all nodes)
	! p ............ distributed load (no ele, 6, no gauss)
	! rot .......... rotations in gauss points (no ele, no gauss, 3, 3)
	! om ........... curvature vector (no ele, 3, no gauss)
	! f ............ internal force vector (no ele, 3, no gauss)
	!
	! modify U, rot, om, f
	subroutine newton_iter (ele, X0, U, C, DOF, dU, Q, p, rot, om, f, TOLER, MAXITER)
	
		implicit none
		
		integer, dimension (:, :), intent (in) :: ele  ! (no ele, no nodes on ele)
		real (DP), dimension (:, :), intent (in) :: X0  ! (3, no all nodes)
		real (DP), dimension (:, :), intent (inout) :: U  ! (3, no all nodes)
		real (DP), dimension (6, 6), intent (in) :: C
		logical, dimension (:, :), intent (in) :: DOF  ! (6, no all nodes)
		real (DP), dimension (:, :), intent (in) :: dU, Q  ! (6, no all nodes)
		real (DP), dimension (:, :, :), intent (in) :: p  ! (no ele, 6, no nodes on ele)
		real (DP), dimension (:, :, :, :), intent (inout) :: rot  ! (no ele, no gauss, 3, 3)
		real (DP), dimension (:, :, :), intent (inout) :: om  ! (no ele, 3, no gauss)
		real (DP), dimension (:, :, :), intent (inout) :: f  ! (no ele, 6, no gauss)
		integer :: nno, nele, ndof, i, j, k, info
		real (DP) :: TOLER, Rnorm
		integer :: MAXITER
		logical, dimension (6 * size (X0 (1, :))) :: DOFsel
		real (DP), dimension (6 * size (X0 (1, :)), 6 * size (X0 (1, :))) :: tangent
		real (DP), allocatable :: K1 (:, :), K2 (:, :)
		real (DP), dimension (6, size (X0 (1, :))) :: res
		real (DP), dimension (6 * size (X0 (1, :))) :: res2
		real (DP), allocatable :: res1 (:), R1 (:), R2 (:, :)
		real (DP), dimension (3, size (X0 (1, :))) :: X, dx, th
		real (DP), dimension (6 * size (X0 (1, :))) :: Fint, Fext, R
		integer, allocatable :: ipiv (:)
		
		nno = size (X0 (1, :))
		nele = size (ele (:, 1))
		
		res = dU
		
		DOFsel = pack (DOF, .TRUE.)
		ndof = count (DOFsel)
		allocate (K1 (6 * nno, ndof), K2 (ndof, ndof), res1 (ndof))
		allocate (R1 (ndof), R2 (ndof, 1), ipiv (ndof))
		
		write (6, '(3X, "n", X, "|", X, "residual vector norm")')
		write (6, '(A)'), repeat("-", 28)
		
		do i = 0, MAXITER-1
			if (i > 0) then
				res = 0.0_DP
				tangent = Kg (ele, X0, X, rot, C, f)  ! tangent
				
				do concurrent (j = 1:6 * nno)  ! keep only dof 1st step
					K1 (j, :) = pack (tangent (j, :), DOFsel)
				end do
				do concurrent (j = 1:ndof)  ! keep only dof 2nd step
					K2 (:, j) = pack (K1 (:, j), DOFsel)
				end do
				
				call dgetrf (ndof, ndof, K2, ndof, ipiv, info)  ! LU factorization
				if (info .ne. 0) then
					write (6, '("Singular matrix")')
					stop
				end if
				
				call dgetrs ('N', ndof, 1, K2, ndof, ipiv, R2, ndof, info)  ! solve system
				if (info .ne. 0) then
					write (6, '("Singular matrix")')
					stop
				end if
				
				res2 = 0.0_DP
				k = 1
				do j = 1, 6 * nno  ! overwrite the result
					if (DOFsel (j) .eq. .TRUE.) then
						res2 (j) = R2 (k, 1)
						k = k + 1
					end if
				end do
				
				res = reshape (res2, (/ 6, nno /))
			end if
			
			do concurrent (j = 1:3)
				dx (j, :) = res (j, :) 
				th (j, :) = res (j + 3, :) 
			end do
			
			U = U + dx
			X = X0 + U
			
			call curv (ele, X0, X, th, C, rot, om, f)
			
			Fint = Fintg (ele, X0, X, f)
			Fext = Fextg (ele, X0, Q, p)
			R = Fint - Fext
			R1 = pack (R, DOFsel)
			R2 (:, 1) = -R1  ! invert and verticalize residual
			Rnorm = norm2 (R1)
			write (6, '(I4, X, "|", X, ES20.13)'), i, Rnorm
			
			if (Rnorm < TOLER) then
				exit
			end if
			
		end do
		
		if (i .eq. MAXITER) then
			write (6, '("Not converging")')
			stop
		end if
		
	end subroutine newton_iter
	
end module