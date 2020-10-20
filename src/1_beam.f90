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
	
	public :: assemble_external_force, assemble_internal_force, update_stress_strain
	
	contains
	
	! compute length of an element
	function element_length (X0, g_ord)
		
		implicit none
		
		integer, intent (in) :: g_ord
		double precision, dimension (:, :), intent (in) :: X0  ! (3, no nodes on ele)
		double precision :: element_length
		integer :: e_ord, i
		double precision, dimension (g_ord) :: pts, wgts
		double precision, dimension (size (X0 (1, :)), g_ord) :: dN  ! (no nodes on ele, no gauss)
		double precision, dimension (3, g_ord) :: dX0  ! (3, no gauss)
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
	function Pi_matrix (rot)
	
		implicit none
		
		double precision, dimension (3, 3), intent (in) :: rot
		double precision, dimension (6, 6) :: Pi_matrix
		
		Pi_matrix = 0.0D0
		Pi_matrix (1:3, 1:3) = rot
		Pi_matrix (4:6, 4:6) = rot
		
	end function Pi_matrix
	
	! compute single Xi matrix
	function Xi_matrix (dX, N, dN)
	
		implicit none
		
		double precision, dimension (3), intent (in) :: dX
		double precision, intent (in) :: N, dN
		double precision, dimension (6, 6) :: Xi_matrix
		double precision, dimension (3, 3) :: S
		integer :: i
		
		S = skew (dX)
		Xi_matrix = 0.0D0
		do i = 1, 6
			Xi_matrix (i, i) = dN
		end do
		Xi_matrix (4, 2) = - N * S (1, 2)
		Xi_matrix (4, 3) = - N * S (1, 3)
		Xi_matrix (5, 1) = - N * S (2, 1)
		Xi_matrix (5, 3) = - N * S (2, 3)
		Xi_matrix (6, 1) = - N * S (3, 1)
		Xi_matrix (6, 2) = - N * S (3, 2)
		
	end function Xi_matrix
	
	! compute elemental internal forces in nodes
	function internal_forces (X0, X, stress)
		
		implicit none
		
		double precision, dimension (:, :), intent (in) :: X0, X  ! (3, no nodes on ele)
		double precision, dimension (:, :), intent (in) :: stress  ! (6, no gauss)
		double precision, dimension (6, size (X0(1,:))) :: internal_forces  ! (6, no nodes on ele)
		integer :: g_ord, nno, e_ord, g, i
		double precision :: L
		double precision, dimension (size (stress (1, :))) :: pts, wgts  ! (no gauss)
		double precision, dimension (size (X (1, :)), size (stress (1, :))) :: N, dN  ! (no nodes on ele, no gauss)
		double precision, dimension (3, size (stress (1, :))) :: dX  ! (3, no gauss)
		double precision, dimension (3, 3) :: S
		double precision, dimension (6, 6) :: Xi_i
		
		g_ord = size (stress (1, :))
		nno = size (X (1, :))
		e_ord = nno - 1
		
		L = element_length (X0, g_ord) 
		call legauss (g_ord, pts, wgts)
		N = shfun (e_ord, pts)
		dN = 2.0D0 / L * shdfun (e_ord, pts)
		
		dX = matmul (X, dN)

		internal_forces = 0.0D0
		
		do g = 1, g_ord
			do i = 1, nno
				Xi_i = Xi_matrix (dX (:, g), N (i, g), dN (i, g))
				internal_forces (:, i) = internal_forces (:, i) + matmul (Xi_i, stress (:, g)) * wgts (g)
			end do
		end do
		
		internal_forces = L / 2.0D0 * internal_forces
	
	end function internal_forces
	
	! compute elemental external forces in nodes
	function external_forces (X0, pressure)
		
		implicit none
		
		double precision, dimension (:, :), intent (in) :: X0  ! (3, no nodes on ele)
		double precision, dimension (:, :), intent (in) :: pressure  ! (6, no gauss)
		double precision, dimension (6, size (X0 (1,:))) :: external_forces  ! (6, no nodes on ele)
		integer :: g_ord, nno, e_ord, g, i
		double precision :: L
		double precision, dimension (size (pressure (1, :))) :: pts, wgts  ! (no gauss)
		double precision, dimension (size (X0 (1, :)), size (pressure (1, :))) :: N  ! (no nodes onele, no gauss)
		double precision, dimension (6, 6) :: H
		
		g_ord = size (pressure (1, :))
		nno = size (X0 (1, :))
		e_ord = nno - 1
		
		L = element_length (X0, g_ord) 
		call legauss (g_ord, pts, wgts)
		N = shfun (e_ord, pts)
		
		external_forces = 0.0D0
		
		do g = 1, g_ord
			H = 0.0D0
			do i = 1, nno
				H (1, 1) = N (i, g)
				H (2, 2) = N (i, g)
				H (3, 3) = N (i, g)
				H (4, 4) = N (i, g)
				H (5, 5) = N (i, g)
				H (6, 6) = N (i, g)
				external_forces (:, i) = external_forces (:, i) + matmul (H, pressure (:, g)) * wgts (g)
			end do
		end do
		
		external_forces = L / 2.0D0 * external_forces
				
	end function external_forces
	
	! compute elemental curvature, internal stresses and new rotation
	subroutine curvature (X0, X, th, C, rot, om, stress)
		
		implicit none
		
		double precision, dimension (:, :), intent (in) :: X0, X, th  ! (3, no nodes on ele)
		double precision, dimension (6, 6), intent (in) :: C
		double precision, dimension (:, :, :), intent (inout) :: rot  ! (no gauss, 3, 3)
		double precision, dimension (:, :), intent (inout) :: om  ! (3, no gauss)
		double precision, dimension (:, :), intent (out) :: stress  ! (6, no gauss)
		integer :: g_ord, nno, e_ord, g
		double precision :: L
		double precision, dimension (size (om (1, :))) :: pts, wgts  ! (no gauss)
		double precision, dimension (size (X (1, :)), size (om (1, :))) :: N, dN  ! (no nodes ele, no gauss)
		double precision, dimension (3, size (om (1, :))) :: dX, t, dt  ! (3, no gauss)
		double precision, dimension (size (om (1, :))) :: tn  ! (no gauss)
		double precision, dimension (3, 3) :: R, rotinv
		double precision, dimension (3) ::  b1, b2, b3, a1, a2, a3, Gamma, kappa, fn, fm
		double precision, dimension (3), parameter :: E3 = (/ 0.0D0, 0.0D0, 1.0D0 /)
		
		g_ord = size (om (1, :))
		nno = size (X0 (1, :))
		e_ord = nno - 1
		
		L = element_length (X0, g_ord)
		call legauss (g_ord, pts, wgts)
		N = shfun (e_ord, pts)
		dN = 2.0D0 / L * shdfun (e_ord, pts)
		
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
				b1 = (1.0D0 - sin (tn (g)) / tn (g)) &
				   * dot_product (t (:, g), dt (:, g)) / tn (g) ** 2 * t (:, g)
				b2 = sin (tn (g)) / tn (g) * dt (:, g)
				b3 = (1.0D0 - cos (tn (g))) / tn (g) ** 2 &
				   * cross_product (t (:, g), dt (:, g))
				a1 = cos (tn (g)) * om (:, g)
				a2 = (1.0D0 - cos (tn (g))) / tn (g) ** 2 &
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
			stress (1:3, g) = matmul (rot (g, :, :), fn)
			stress (4:6, g) = matmul (rot (g, :, :), fm)
			
		end do
		
	end subroutine curvature
	
	! compute stiffness of an element
	function material_stiffness (X0, X, rot, C)
		
		implicit none
		
		double precision, dimension (:, :), intent (in) :: X0, X  ! (3, no nodes on ele)
		double precision, dimension (:, :, :), intent (in) :: rot  ! (no gauss, 3, 3)
		double precision, dimension (6, 6), intent (in) :: C
		double precision, dimension (6 * size (X (1, :)), 6 * size (X (1, :))) :: material_stiffness  ! (6*no nodes on ele, 
                                                                                                      !  6*no nodes on ele)
		integer :: g_ord, nno, e_ord, g, i, j, ii, ij
		double precision :: L
		double precision, dimension (size (rot (:, 1, 1))) :: pts, wgts  ! (no gauss)
		double precision, dimension (size (X (1, :)), size (rot (:, 1, 1))) :: N, dN  ! (no nodes on ele, no gauss)
		double precision, dimension (3, size (rot (:, 1, 1))) :: dX  ! (3, no gauss)
		double precision, dimension (6, 6) :: Pi_g, Xi_i, Xi_j, c_g
		
		g_ord = size (rot (:, 1, 1))
		nno = size (X0 (1, :))
		e_ord = nno - 1
		
		L = element_length (X0, g_ord)
		call legauss (g_ord, pts, wgts)
		N = shfun (e_ord, pts)
		dN = 2.0D0 / L * shdfun (e_ord, pts)
		
		dX = matmul (X, dN)
		
		material_stiffness = 0.0D0
		
		do g = 1, g_ord
			Pi_g = Pi_matrix (rot (g, :, :))
			c_g = matmul (matmul (Pi_g, C), transpose (Pi_g))
			do i = 1, nno
				ii = 6 * (i-1) + 1
				Xi_i = Xi_matrix (dX (:, g), N (i, g), dN (i, g))
				do j = 1, nno
					Xi_j = Xi_matrix (dX (:, g), N (j, g), dN (j, g))
					ij = 6 * (j-1) + 1
					material_stiffness (ii:ii + 5, ij:ij + 5) = material_stiffness (ii:ii + 5, ij:ij + 5) + &
						wgts (g) * matmul (matmul (Xi_i, c_g), transpose (Xi_j))
				end do
			end do
		end do
		
		material_stiffness = L / 2.0D0 * material_stiffness
	
	end function material_stiffness
	
	! compute stiffness of an element
	function geometrical_stiffness (X0, X, rot, stress)
		
		implicit none
		
		double precision, dimension (:, :), intent (in) :: X0, X  ! (3, no nodes on ele)
		double precision, dimension (:, :, :), intent (in) :: rot  ! (no gauss, 3, 3)
		double precision, dimension (:, :), intent (in) :: stress  ! (6, no gauss)
		double precision, dimension (6 * size (X (1, :)), 6 * size (X (1, :))) :: geometrical_stiffness  ! (6*no nodes on ele, 
                                                                                                         !  6*no nodes on ele)
		integer :: g_ord, nno, e_ord, g, i, j, ii, ij
		double precision :: L
		double precision, dimension (size (rot (:, 1, 1))) :: pts, wgts  ! (no gauss)
		double precision, dimension (size (X (1, :)), size (rot (:, 1, 1))) :: N, dN  ! (no nodes on ele, no gauss)
		double precision, dimension (6, size (rot (:, 1, 1))) :: dX  ! (6, no gauss)
		
		g_ord = size (rot (:, 1, 1))
		nno = size (X0 (1, :))
		e_ord = nno - 1
		
		L = element_length (X0, g_ord)
		call legauss (g_ord, pts, wgts)
		N = shfun (e_ord, pts)
		dN = 2.0D0 / L * shdfun (e_ord, pts)
		
		dX = matmul (X, dN)
		
		geometrical_stiffness = 0.0D0
		
		do g = 1, g_ord
			do i = 1, nno
				ii = 6 * (i-1) + 1
				do j = 1, nno
					ij = 6 * (j-1) + 1
					geometrical_stiffness (ii:ii + 2, ij + 3:ij + 5) = &
						geometrical_stiffness (ii:ii + 2, ij + 3:ij + 5) - &
						wgts (g) * skew (stress (1:3, g)) * dN (i, g) * N (j, g)
					geometrical_stiffness (ii + 3:ii + 5, ij:ij + 2) = &
						geometrical_stiffness (ii + 3:ii + 5, ij:ij + 2) + &
						wgts (g) * skew (stress (1:3, g)) * N (i, g) * dN (j, g)
					geometrical_stiffness (ii + 3:ii + 5, ij + 3:ij + 5) = &
						geometrical_stiffness (ii + 3:ii + 5, ij + 3:ij + 5) + &
						wgts (g) * ( &
							-skew (stress (4:6, g)) * dN (i, g) * N (j, g) &
							+ matmul (skew (dX (:, g)), skew (stress (1:3, g))) &
							* N (i, g) * N (j, g) &
						)
				end do
			end do
		end do
		
		geometrical_stiffness = L / 2.0D0 * geometrical_stiffness
	
	end function geometrical_stiffness
	
	! compute stiffness of an element
	function assemble_tangent (ele, X0, X, rot, C, stress) result (Kg)
	
		implicit none
		
		integer, dimension (:, :), intent (in) :: ele  ! (no ele, no nodes on ele)
		double precision, dimension (:, :), intent (in) :: X0, X  ! (3, no nodes)
		double precision, dimension (:, :, :, :), intent (in) :: rot  ! (no ele, no gauss, 3, 3)
		double precision, dimension (6, 6), intent (in) :: C
		double precision, dimension (:, :, :), intent (in) :: stress  ! (no ele, 6, no gauss)
		double precision, dimension (6 * size (X (1, :)), 6 * size (X (1, :))) :: Kg  ! (6*no nodes, 6*no nodes)
		integer :: nno, nele, neno, e, i, j, ei, ej
		integer, dimension (size (ele (1, :))) :: ee  ! (no nodes on ele)
		double precision, dimension (6 * size (ele (1, :)), 6 * size (ele (1, :))) :: Ke  ! (6*no nodes on ele, 
                                                                                           !  6*no nodes on ele)
        
		nno = size (X0 (1, :))
		nele = size (ele (:, 1))
				
		Kg = 0.0D0
		
		do e = 1, nele
			ee = ele (e, :)
			neno = size (ee)
			Ke = material_stiffness (X0 (:, ee), X (:, ee), rot (e, :, :, :), C) &
				+ geometrical_stiffness (X0 (:, ee), X (:, ee), rot (e, :, :, :), stress (e, :, :))
			do i = 1, neno
				ei = ee (i)
				do j = 1, neno
					ej = ee (j)
					Kg (6 * (ei - 1) + 1:6 * ei, 6 * (ej - 1) + 1:6 * ej) = &
						Kg (6 * (ei - 1) + 1:6 * ei, 6 * (ej - 1) + 1:6 * ej) &
						+ Ke (6 * (i - 1) + 1:6 * i, 6 * (j - 1) + 1:6 * j)
				end do
			end do
		end do
		
	end function assemble_tangent
	
	! compute global internal force vector
	function assemble_internal_force (ele, X0, X, stress) return (Fint)
	
		implicit none
		
		integer, dimension (:, :), intent (in) :: ele  ! (no ele, no nodes on ele)
		double precision, dimension (:, :), intent (in) :: X0, X  ! (3, no nodes)
		double precision, dimension (:, :, :), intent (in) :: stress  ! (no ele, 6, no gauss)
		double precision, dimension (6 * size (X (1, :))) :: Fint  ! (6*no nodes)
		integer :: nno, nele, neno, e, i, ei
		integer, dimension (size (ele (1, :))) :: ee
		double precision, dimension (6, size (ele (1, :))) :: Fie
		
		nno = size (X0 (1, :))
		nele = size (ele (:, 1))
		
		Fint = 0.0D0
		
		do e = 1, nele
			ee = ele (e, :)
			neno = size (ee)
			Fie = internal_forces (X0 (:, ee), X (:, ee), stress (e, :, :))
			do i = 1, neno
				ei = ee (i)
				Fint (6 * (ei - 1) + 1:6 * ei) = &
					Fint (6 * (ei - 1) + 1:6 * ei) + Fie (:, i)
			end do
		end do
	
	end function Fint
	
	! compute global external force vector
	function assemble_external_force (ele, X0, Q, pressure) return (Fext)
	
		implicit none
		
		integer, dimension (:, :), intent (in) :: ele
		double precision, dimension (:, :), intent (in) :: X0
		double precision, dimension (:, :), intent (in) :: Q
		double precision, dimension (:, :, :), intent (in) :: pressure
		double precision, dimension (6 * size (X0 (1, :))) :: Fext
		integer :: nno, nele, neno, e, i, ei
		integer, dimension (size (ele (1, :))) :: ee
		double precision, dimension (6, size (ele (1, :))) :: Fee
		
		nno = size (X0 (1, :))
		nele = size (ele (:, 1))
		
		Fext = 0.0D0
		
		do e = 1, nele
			ee = ele (e, :)
			neno = size (ee)
			Fee = external_forces (X0 (:, ee), pressure (e, :, :))
			do i = 1, neno
				ei = ee (i)
				Fext (6 * (ei - 1) + 1:6 * ei) = &
					Fext (6 * (ei - 1) + 1:6 * ei) + Fee (:, i)			
			end do
		end do
		
		Fext = Fext + reshape (Q, (/ 6 * nno /))
	
	end function Fext
	
	! compute global external force vector
	subroutine update_stress_strain (ele, X0, X, th, C, rot, om, stress)
	
		implicit none
		
		integer, dimension (:, :), intent (in) :: ele  ! (no ele, no nodes on ele)
		double precision, dimension (:, :), intent (in) :: X0, X, th  ! (3, no all nodes)
		double precision, dimension (6, 6), intent (in) :: C
		double precision, dimension (:, :, :, :), intent (inout) :: rot  ! (no ele, no gauss, 3, 3)
		double precision, dimension (:, :, :), intent (inout) :: om  ! (no ele, 3, no gauss)
		double precision, dimension (:, :, :), intent (out) :: stress  ! (no ele, 6, no gauss)
		integer :: nno, nele, e
		integer, dimension (size (ele (1, :))) :: ee
		
		nno = size (X0 (1, :))
		nele = size (ele (:, 1))
		
		do e = 1, nele
			ee = ele (e, :)
			call curvature ( &
				X0 (:, ee), X (:, ee), th (:, ee), C, &
				rot (e, :, :, :), om (e, :, :), stress (e, :, :) &
			)
		end do
		
	end subroutine curv
	
	
end module
