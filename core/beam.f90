! Beam element
!
! author ....... Jan Tomec
! copyright .... Copyright 2020, Project THREAD - University of Rijeka, Faculty of Civil Engineering
! credits ...... Jan Tomec, Gordan Jelenić
! license ...... GPL
! version ...... 1.0.0
! maintainer ... Jan Tomec
! email ........ jan.tomec@gradri.uniri.hr
! status ....... Development
! date ......... 06/05/2020
!
! ------------------------------------------------------------------------------

module beam
	
	use shape_functions
	use legendre_gauss
	use vector_algebra
	
	implicit none
	
	private
	
	public :: Kg, Fintg, Fextg, curv
	
	contains
	
	! compute length of an element
	function element_length (X0, g_ord)
		
		implicit none
		
		integer, intent (in) :: g_ord
		double precision, dimension (:, :), intent (in) :: X0
		double precision :: element_length
		integer :: e_ord, i
		double precision, dimension (g_ord) :: pts, wgts
		double precision, dimension (size (X0 (1, :)), g_ord) :: dN
		double precision, dimension (3, g_ord) :: dX0
		double precision, dimension (3) :: intg
		
		e_ord = size (X0 (1, :)) - 1
		
		call legauss (g_ord, pts, wgts)
		dN = shdfun (e_ord, pts)
		
		dX0 = matmul (X0, dN)
		
		do i = 1, 3
			intg (i) = dot_product (dX0 (i, :), wgts)
		end do
		element_length = norm2 (intg)
	
	end function element_length
	
	! compute single Pi matrix
	function Pi (rot)
	
		implicit none
		
		double precision, dimension (3, 3), intent (in) :: rot
		double precision, dimension (6, 6) :: Pi
		
		Pi = 0.0_DP
		Pi (1:3, 1:3) = rot
		Pi (4:6, 4:6) = rot
		
	end function Pi
	
	! compute single Xi matrix
	function Xi (dX, N, dN)
	
		implicit none
		
		double precision, dimension (:), intent (in) :: dX
		double precision, intent (in) :: N, dN
		double precision, dimension (6, 6) :: Xi
		double precision, dimension (3, 3) :: S
		integer :: i
		
		S = skew (dX)
		Xi = 0.0_DP
		do i = 1, 6
			Xi (i, i) = dN
		end do
		Xi (4, 2) = - N * S (1, 2)
		Xi (4, 3) = - N * S (1, 3)
		Xi (5, 1) = - N * S (2, 1)
		Xi (5, 3) = - N * S (2, 3)
		Xi (6, 1) = - N * S (3, 1)
		Xi (6, 2) = - N * S (3, 2)
		
	end function Xi
	
	! compute elemental internal forces in nodes
	function Finte (X0, X, f)
		
		implicit none
		
		double precision, dimension (:, :), intent (in) :: X0, X
		double precision, dimension (:, :), intent (in) :: f
		double precision, dimension (6, size (X0(1,:))) :: Finte
		integer :: g_ord, nno, e_ord, g, i
		double precision :: L
		double precision, dimension (size (f (1, :))) :: pts, wgts
		double precision, dimension (size (X (1, :)), size (f (1, :))) :: N, dN
		double precision, dimension (3, size (f (1, :))) :: dX
		double precision, dimension (3, 3) :: S
		double precision, dimension (6, 6) :: Xi_i
		
		g_ord = size (f (1, :))
		nno = size (X (1, :))
		e_ord = nno - 1
		
		L = element_length (X0, g_ord) 
		call legauss (g_ord, pts, wgts)
		N = shfun (e_ord, pts)
		dN = 2.0_DP / L * shdfun (e_ord, pts)
		
		dX = matmul (X, dN)

		Finte = 0.0_DP
		
		do g = 1, g_ord
			do i = 1, nno
				Xi_i = Xi (dX (:, g), N (i, g), dN (i, g))
				Finte (:, i) = Finte (:, i) + matmul (Xi_i, f (:, g)) * wgts (g)
			end do
		end do
		
		Finte = L / 2.0_DP * Finte
	
	end function Finte
	
	! compute elemental external forces in nodes
	function Fexte (X0, p)
		
		implicit none
		
		double precision, dimension (:, :), intent (in) :: X0
		double precision, dimension (:, :), intent (in) :: p
		double precision, dimension (6, size (X0 (1,:))) :: Fexte
		integer :: g_ord, nno, e_ord, g, i
		double precision :: L
		double precision, dimension (size (p (1, :))) :: pts, wgts
		double precision, dimension (size (X0 (1, :)), size (p (1, :))) :: N
		double precision, dimension (6, 6) :: H
		
		g_ord = size (p (1, :))
		nno = size (X0 (1, :))
		e_ord = nno - 1
		
		L = element_length (X0, g_ord) 
		call legauss (g_ord, pts, wgts)
		N = shfun (e_ord, pts)
		
		Fexte = 0.0_DP
		
		do g = 1, g_ord
			H = 0.0_DP
			do i = 1, nno
				H (1, 1) = N (i, g)
				H (2, 2) = N (i, g)
				H (3, 3) = N (i, g)
				H (4, 4) = N (i, g)
				H (5, 5) = N (i, g)
				H (6, 6) = N (i, g)
				Fexte (:, i) = Fexte (:, i) + matmul (H, p (:, g)) * wgts (g)
			end do
		end do
		
		Fexte = L / 2.0_DP * Fexte
				
	end function Fexte
	
	! compute elemental curvature, internal stresses and new rotation
	subroutine curve (X0, X, th, C, rot, om, f)
		
		implicit none
		
		double precision, dimension (:, :), intent (in) :: X0, X, th
		double precision, dimension (6, 6), intent (in) :: C
		double precision, dimension (:, :, :), intent (inout) :: rot
		double precision, dimension (:, :), intent (inout) :: om
		double precision, dimension (:, :), intent (out) :: f
		integer :: g_ord, nno, e_ord, g
		double precision :: L
		double precision, dimension (size (om (1, :))) :: pts, wgts
		double precision, dimension (size (X (1, :)), size (om (1, :))) :: N, dN
		double precision, dimension (3, size (om (1, :))) :: dX, t, dt
		double precision, dimension (size (om (1, :))) :: tn
		double precision, dimension (3, 3) :: R, rotinv
		double precision, dimension (3) ::  b1, b2, b3, a1, a2, a3, Gamma, kappa, fn, fm
		double precision, dimension (3), parameter :: E3 = (/ 0.0_DP, 0.0_DP, 1.0_DP /)
		
		g_ord = size (om (1, :))
		nno = size (X0 (1, :))
		e_ord = nno - 1
		
		L = element_length (X0, g_ord)
		call legauss (g_ord, pts, wgts)
		N = shfun (e_ord, pts)
		dN = 2.0_DP / L * shdfun (e_ord, pts)
		
		dX = matmul (X, dN)
		t = matmul (th, N)
		dt = matmul (th, dN)
		tn = norm2 (t, dim=1)
		
		do g = 1, g_ord
			R = rv2mat (t (:, g))
			rot (g, :, :) = matmul (R, rot (g, :, :))
			if (tn (g) == 0) then
				om (:, g) = om (:, g) + dt (:, g)
			else
				b1 = (1.0_DP - sin (tn (g)) / tn (g)) &
				   * dot_product (t (:, g), dt (:, g)) / tn (g) ** 2 * t (:, g)
				b2 = sin (tn (g)) / tn (g) * dt (:, g)
				b3 = (1.0_DP - cos (tn (g))) / tn (g) ** 2 &
				   * cross_product (t (:, g), dt (:, g))
				a1 = cos (tn (g)) * om (:, g)
				a2 = (1.0_DP - cos (tn (g))) / tn (g) ** 2 &
				   * dot_product (t (:, g), om (:, g)) * t (:, g)
				a3 = sin (tn (g)) / tn (g) * cross_product (t (:, g), om (:, g))
				om (:, g) = b1 + b2 + b3 + a1 + a2 + a3
			end if
			
			rotinv = transpose (rot (g, :, :))
			Gamma = matmul (rotinv, dX (:, g))
			Gamma = Gamma - E3
			kappa = matmul (rotinv, om (:, g))
			fn = matmul (C (1:3, 1:3), Gamma)
			fm = matmul (C (4:6, 4:6), kappa)
			f (1:3, g) = matmul (rot (g, :, :), fn)
			f (4:6, g) = matmul (rot (g, :, :), fm)
			
		end do
		
	end subroutine curve
	
	! compute stiffness of an element
	function Ke (X0, X, rot, C)
		
		implicit none
		
		double precision, dimension (:, :), intent (in) :: X0, X
		double precision, dimension (:, :, :), intent (in) :: rot
		double precision, dimension (6, 6), intent (in) :: C
		double precision, dimension (6 * size (X (1, :)), 6 * size (X (1, :))) :: Ke
		integer :: g_ord, nno, e_ord, g, i, j, ii, ij
		double precision :: L
		double precision, dimension (size (rot (:, 1, 1))) :: pts, wgts
		double precision, dimension (size (X (1, :)), size (rot (:, 1, 1))) :: N, dN
		double precision, dimension (3, size (rot (:, 1, 1))) :: dX
		double precision, dimension (6, 6) :: Pi_g, Xi_i, Xi_j, c_g
		
		g_ord = size (rot (:, 1, 1))
		nno = size (X0 (1, :))
		e_ord = nno - 1
		
		L = element_length (X0, g_ord)
		call legauss (g_ord, pts, wgts)
		N = shfun (e_ord, pts)
		dN = 2.0_DP / L * shdfun (e_ord, pts)
		
		dX = matmul (X, dN)
		
		Ke = 0.0_DP
		
		do g = 1, g_ord
			Pi_g = Pi (rot (g, :, :))
			c_g = matmul (matmul (Pi_g, C), transpose (Pi_g))
			do i = 1, nno
				ii = 6 * (i-1) + 1
				Xi_i = Xi (dX (:, g), N (i, g), dN (i, g))
				do j = 1, nno
					Xi_j = Xi (dX (:, g), N (j, g), dN (j, g))
					ij = 6 * (j-1) + 1
					Ke (ii:ii + 5, ij:ij + 5) = Ke (ii:ii + 5, ij:ij + 5) + &
						wgts (g) * matmul (matmul (Xi_i, c_g), transpose (Xi_j))
				end do
			end do
		end do
		
		Ke = L / 2.0_DP * Ke
	
	end function Ke
	
	! compute stiffness of an element
	function Kgeoe (X0, X, rot, f)
		
		implicit none
		
		double precision, dimension (:, :), intent (in) :: X0, X
		double precision, dimension (:, :, :), intent (in) :: rot
		double precision, dimension (:, :), intent (in) :: f
		double precision, dimension (6 * size (X (1, :)), 6 * size (X (1, :))) :: Kgeoe
		integer :: g_ord, nno, e_ord, g, i, j, ii, ij
		double precision :: L
		double precision, dimension (size (rot (:, 1, 1))) :: pts, wgts
		double precision, dimension (size (X (1, :)), size (rot (:, 1, 1))) :: N, dN
		double precision, dimension (6, size (rot (:, 1, 1))) :: dX
		
		g_ord = size (rot (:, 1, 1))
		nno = size (X0 (1, :))
		e_ord = nno - 1
		
		L = element_length (X0, g_ord)
		call legauss (g_ord, pts, wgts)
		N = shfun (e_ord, pts)
		dN = 2.0_DP / L * shdfun (e_ord, pts)
		
		dX = matmul (X, dN)
		
		Kgeoe = 0.0_DP
		
		do g = 1, g_ord
			do i = 1, nno
				ii = 6 * (i-1) + 1
				do j = 1, nno
					ij = 6 * (j-1) + 1
					Kgeoe (ii:ii + 2, ij + 3:ij + 5) = &
						Kgeoe (ii:ii + 2, ij + 3:ij + 5) + &
						wgts (g) * -skew (f (1:3, g)) * dN (i, g) * N (j, g)
					Kgeoe (ii + 3:ii + 5, ij:ij + 2) = &
						Kgeoe (ii + 3:ii + 5, ij:ij + 2) + &
						wgts (g) * skew (f (1:3, g)) * N (i, g) * dN (j, g)
					Kgeoe (ii + 3:ii + 5, ij + 3:ij + 5) = &
						Kgeoe (ii + 3:ii + 5, ij + 3:ij + 5) + &
						wgts (g) * ( &
							-skew (f (4:6, g)) * dN (i, g) * N (j, g) &
							+ matmul (skew (dX (:, g)), skew (f (1:3, g))) &
							* N (i, g) * N (j, g) &
						)
				end do
			end do
		end do
		
		Kgeoe = L / 2.0_DP * Kgeoe
	
	end function Kgeoe
	
	! compute stiffness of an element
	function Kg (ele, X0, X, rot, C, f)
	
		implicit none
		
		integer, dimension (:, :), intent (in) :: ele
		double precision, dimension (:, :), intent (in) :: X0, X
		double precision, dimension (:, :, :, :), intent (in) :: rot
		double precision, dimension (6, 6), intent (in) :: C
		double precision, dimension (:, :, :), intent (in) :: f
		double precision, dimension (6 * size (X (1, :)), 6 * size (X (1, :))) :: Kg
		integer :: nno, nele, neno, e, i, j, ei, ej
		integer, dimension (size (ele (1, :))) :: ee
		double precision, dimension (6 * size (ele (1, :)), 6 * size (ele (1, :))) :: Kee
				
		nno = size (X0 (1, :))
		nele = size (ele (:, 1))
				
		Kg = 0.0_DP
		
		do e = 1, nele
			ee = ele (e, :)
			neno = size (ee)
			Kee = Ke (X0 (:, ee), X (:, ee), rot (e, :, :, :), C) &
				+ Kgeoe (X0 (:, ee), X (:, ee), rot (e, :, :, :), f (e, :, :))
			do i = 1, neno
				ei = ee (i)
				do j = 1, neno
					ej = ee (j)
					Kg (6 * (ei - 1) + 1:6 * ei, 6 * (ej - 1) + 1:6 * ej) = &
						Kg (6 * (ei - 1) + 1:6 * ei, 6 * (ej - 1) + 1:6 * ej) &
						+ Kee (6 * (i - 1) + 1:6 * i, 6 * (j - 1) + 1:6 * j)
				end do
			end do
		end do
		
	end function Kg
	
	! compute global internal force vector
	function Fintg (ele, X0, X, f)
	
		implicit none
		
		integer, dimension (:, :), intent (in) :: ele
		double precision, dimension (:, :), intent (in) :: X0, X
		double precision, dimension (:, :, :), intent (in) :: f
		double precision, dimension (6 * size (X (1, :))) :: Fintg
		integer :: nno, nele, neno, e, i, ei
		integer, dimension (size (ele (1, :))) :: ee
		double precision, dimension (6, size (ele (1, :))) :: Fie
		
		nno = size (X0 (1, :))
		nele = size (ele (:, 1))
		
		Fintg = 0.0_DP
		
		do e = 1, nele
			ee = ele (e, :)
			neno = size (ee)
			Fie = Finte (X0 (:, ee), X (:, ee), f (e, :, :))
			do i = 1, neno
				ei = ee (i)
				Fintg (6 * (ei - 1) + 1:6 * ei) = &
					Fintg (6 * (ei - 1) + 1:6 * ei) + Fie (:, i)
			end do
		end do
	
	end function Fintg
	
	! compute global external force vector
	function Fextg (ele, X0, Q, p)
	
		implicit none
		
		integer, dimension (:, :), intent (in) :: ele
		double precision, dimension (:, :), intent (in) :: X0
		double precision, dimension (:, :), intent (in) :: Q
		double precision, dimension (:, :, :), intent (in) :: p
		double precision, dimension (6 * size (X0 (1, :))) :: Fextg
		integer :: nno, nele, neno, e, i, ei
		integer, dimension (size (ele (1, :))) :: ee
		double precision, dimension (6, size (ele (1, :))) :: Fee
		
		nno = size (X0 (1, :))
		nele = size (ele (:, 1))
		
		Fextg = 0.0_DP
		
		do e = 1, nele
			ee = ele (e, :)
			neno = size (ee)
			Fee = Fexte (X0 (:, ee), p (e, :, :))
			do i = 1, neno
				ei = ee (i)
				Fextg (6 * (ei - 1) + 1:6 * ei) = &
					Fextg (6 * (ei - 1) + 1:6 * ei) + Fee (:, i)			
			end do
		end do
		
		Fextg = Fextg + reshape (Q, (/ 6 * nno /))
	
	end function Fextg
	
	! compute global external force vector
	subroutine curv (ele, X0, X, th, C, rot, om, f)
	
		implicit none
		
		integer, dimension (:, :), intent (in) :: ele  ! (no ele, no nodes on ele)
		double precision, dimension (:, :), intent (in) :: X0, X, th  ! (3, no all nodes)
		double precision, dimension (6, 6), intent (in) :: C
		double precision, dimension (:, :, :, :), intent (inout) :: rot  ! (no ele, no gauss, 3, 3)
		double precision, dimension (:, :, :), intent (inout) :: om  ! (no ele, 3, no gauss)
		double precision, dimension (:, :, :), intent (out) :: f  ! (no ele, 6, no gauss)
		integer :: nno, nele, e
		integer, dimension (size (ele (1, :))) :: ee
		
		nno = size (X0 (1, :))
		nele = size (ele (:, 1))
		
		do e = 1, nele
			ee = ele (e, :)
			call curve ( &
				X0 (:, ee), X (:, ee), th (:, ee), C, &
				rot (e, :, :, :), om (e, :, :), f (e, :, :) &
			)
		end do
		
	end subroutine curv
	
	
end module
