! Algebraic routines
!
! author ....... Jan Tomec
! copyright .... Copyright 2020, Project THREAD - University of Rijeka, Faculty of Civil Engineering
! credits ...... Jan Tomec, Gordan Jelenić
! license ...... GPL
! version ...... 1.0.0
! maintainer ... Jan Tomec
! email ........ jan.tomec@gradri.uniri.hr
! status ....... Development
! date ......... 06/04/2020
!
! ------------------------------------------------------------------------------

module vector_algebra
	
	use f95_precision

	implicit none
	
	private
	
	public :: cross_product, tensor_product, skew, rv2mat
	
	contains
	
	! compute cross product
	pure function cross_product(a, b)
	
		implicit none
		
		real (DP), dimension (3), intent (in) :: a, b
		real (DP), dimension (3) :: cross_product
		
		cross_product(1) = a(2) * b(3) - a(3) * b(2)
		cross_product(2) = a(3) * b(1) - a(1) * b(3)
		cross_product(3) = a(1) * b(2) - a(2) * b(1)
		
	end function cross_product
	
	! compute tensor product
	pure function tensor_product(a, b)
	
		implicit none
		
		real (DP), dimension (3), intent (in) :: a, b
		real (DP), dimension (3, 3) :: tensor_product
		integer :: i, j
		
		do i = 1, 3
			do j = 1, 3
				tensor_product (i, j) = a (i) * b (j)
			end do
		end do
		
	end function tensor_product
	
	! create skew-symmetric matrix from rotation vector
	! from Argyris_J.H.--An_excursion_into_large_rotations, p.88
	pure function skew (r)
		
		implicit none
		
		real (DP), dimension (3), intent (in) :: r
		real (DP), dimension (3, 3) :: skew
		
		skew (1, 1) = 0.0_DP
		skew (1, 2) = - r (3)
		skew (1, 3) = r (2)
		skew (2, 1) = r (3)
		skew (2, 2) = 0.0_DP
		skew (2, 3) = - r (1)
		skew (3, 1) = - r (2)
		skew (3, 2) = r (1)
		skew (3, 3) = 0.0_DP
		
	end function skew
	
	! create skew-symmetric matrix squared from rotation vector
	! from Argyris_J.H.--An_excursion_into_large_rotations, p.88
	pure function skew2 (r)
	
		implicit none
		
		real (DP), dimension (3), intent (in) :: r
		real (DP), dimension (3, 3) :: skew2
		
		skew2 (1, 1) = - (r (2) ** 2 + r (3) ** 2)
		skew2 (1, 2) = r (1) * r (2)
		skew2 (1, 3) = r (1) * r (3)
		skew2 (2, 1) = r (1) * r (2)
		skew2 (2, 2) = - (r (1) ** 2 + r (3) ** 2)
		skew2 (2, 3) = r (2) * r (3)
		skew2 (3, 1) = r (1) * r (3)
		skew2 (3, 2) = r (2) * r (3)
		skew2 (3, 3) = - (r (1) ** 2 + r (2) ** 2)
		
	end function skew2
	
	! create rotation matrix from rotation vector
	! from Argyris_J.H.--An_excursion_into_large_rotations, p.88
	function rv2mat (r) result (T)
		
		implicit none
		
		real (DP), dimension (3), intent (in) :: r
		real (DP), dimension (3, 3) :: T, T1, T2, I, S, S2
		real (DP) :: sin_r, norm_r, sin_hr, norm_hr
		integer :: j
		
		I = 0.0_DP
		do j = 1, 3
			I (j, j) = 1.0_DP
		end do
		
		norm_r = norm2 (r)
		if (norm_r .eq. 0.0_DP) then
			T = I
			return
		end if
		
		sin_r = sin (norm_r)
		norm_hr = norm2 (r / 2.0_DP)
		sin_hr = sin (norm_hr)
		
		S = skew (r)
		S2 = skew2 (r)
		
		T1 = sin_r / norm_r * S
		T2 = 0.5_DP * (sin_hr / norm_hr) ** 2 * S2
		T = I + T1 + T2
		
	end function rv2mat
	
end module