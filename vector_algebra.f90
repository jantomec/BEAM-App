module vector_algebra

	implicit none	
	
	contains
	
	! create skew-symmetric matrix from rotation vector
	! from Argyris_J.H.--An_excursion_into_large_rotations, p.88
	function skew (r)
		
		implicit none
		
		integer, parameter :: DP = selected_real_kind (p=15, r=307)
		
		real (DP), dimension (3) :: r
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
	function skew2 (r)
	
		implicit none
		
		integer, parameter :: DP = selected_real_kind (p=15, r=307)
		
		real (DP), dimension (3) :: r
		real (DP), dimension (3, 3) :: skew2
		
		skew2 (1, 1) = - (r (2) ** 2.0_DP + r (3) ** 2.0_DP)
		skew2 (1, 2) = r (1) * r (2)
		skew2 (1, 3) = r (1) * r (3)
		skew2 (2, 1) = r (1) * r (2)
		skew2 (2, 2) = - (r (1) ** 2.0_DP + r (3) ** 2.0_DP)
		skew2 (2, 3) = r (2) * r (3)
		skew2 (3, 1) = r (1) * r (3)
		skew2 (3, 2) = r (2) * r (3)
		skew2 (3, 3) = - (r (1) ** 2.0_DP + r (2) ** 2.0_DP)
		
	end function skew2
	
	! create rotation matrix from rotation vector
	! from Argyris_J.H.--An_excursion_into_large_rotations, p.88
	function rv2mat (r) result (T)
		
		implicit none
		
		integer, parameter :: DP = selected_real_kind (p=15, r=307)
		
		integer :: j
		real (DP) :: sin_r, norm_r, sin_hr, norm_hr
		real (DP), dimension (3) :: r
		real (DP), dimension (3, 3) :: T, T1, T2, I, S, S2
		
		norm_r = norm2 (r)
		sin_r = sin (norm_r)
		norm_hr = norm2 (r / 2.0_DP)
		sin_hr = sin (norm_hr)
		
		I = 0.0_DP
		do j = 1, 3
			I (j, j) = 1.0_DP
		end do
		
		S = skew (r)
		S2 = skew2 (r)
		
		T1 = sin_r / norm_r * S
		T2 = 0.5_DP * (sin_hr / norm_hr) ** 2.0_DP * S2
		T = I + T1 + T2
		
	end function rv2mat
	
end module