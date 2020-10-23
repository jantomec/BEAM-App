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
	
	implicit none
	
	private
	
	public :: cross_product, tensor_product, skew, exponentialMap, diagonalMatrix, identityMatrix
	
	contains
	
	! compute cross product
	pure function cross_product(a, b)
	
		implicit none
		
		double precision, dimension (3), intent (in) :: a, b
		double precision, dimension (3) :: cross_product
		
		cross_product(1) = a(2) * b(3) - a(3) * b(2)
		cross_product(2) = a(3) * b(1) - a(1) * b(3)
		cross_product(3) = a(1) * b(2) - a(2) * b(1)
		
	end function cross_product
	
	! compute tensor product
	pure function tensor_product(a, b)
	
		implicit none
		
		double precision, dimension (3), intent (in) :: a, b
		double precision, dimension (3, 3) :: tensor_product
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
		
		double precision, dimension (3), intent (in) :: r
		double precision, dimension (3, 3) :: skew
		
		skew (1, 1) = 0.0D0
		skew (1, 2) = - r (3)
		skew (1, 3) = r (2)
		skew (2, 1) = r (3)
		skew (2, 2) = 0.0D0
		skew (2, 3) = - r (1)
		skew (3, 1) = - r (2)
		skew (3, 2) = r (1)
		skew (3, 3) = 0.0D0
		
	end function skew
    
    pure function antiskew (R)
		
		implicit none
		
		double precision, dimension (3, 3), intent (in) :: R
		double precision, dimension (3) :: antiskew
		
		antiskew (1) = R (3, 2)
		antiskew (2) = R (1, 3)
		antiskew (3) = R (2, 3)
		
	end function antiskew
	
	! create skew-symmetric matrix squared from rotation vector
	! from Argyris_J.H.--An_excursion_into_large_rotations, p.88
	! pure function skew2 (r)
	
		! implicit none
		
		! double precision, dimension (3), intent (in) :: r
		! double precision, dimension (3, 3) :: skew2
		
		! skew2 (1, 1) = - (r (2) ** 2 + r (3) ** 2)
		! skew2 (1, 2) = r (1) * r (2)
		! skew2 (1, 3) = r (1) * r (3)
		! skew2 (2, 1) = r (1) * r (2)
		! skew2 (2, 2) = - (r (1) ** 2 + r (3) ** 2)
		! skew2 (2, 3) = r (2) * r (3)
		! skew2 (3, 1) = r (1) * r (3)
		! skew2 (3, 2) = r (2) * r (3)
		! skew2 (3, 3) = - (r (1) ** 2 + r (2) ** 2)
		
	! end function skew2
	
	! create rotation matrix from rotation vector
	! from Argyris_J.H.--An_excursion_into_large_rotations, p.88
	function exponentialMap (r) result (T)
		
		implicit none
		
		double precision, dimension (3), intent (in) :: r
		double precision, dimension (3, 3) :: T, T1, T2, I, S, S2
		double precision :: sin_r, norm_r, sin_hr, norm_hr
		integer :: j
		
		I = 0.0D0
		do j = 1, 3
			I (j, j) = 1.0D0
		end do
		
		norm_r = norm2 (r)
		if (norm_r .eq. 0.0D0) then
			T = I
			return
		end if
		
		sin_r = sin (norm_r)
		norm_hr = norm2 (r / 2.0D0)
		sin_hr = sin (norm_hr)
		
		S = skew (r)
		S2 = matmul (S, S)
		
		T1 = sin_r / norm_r * S
		T2 = 0.5D0 * (sin_hr / norm_hr) ** 2 * S2
		T = I + T1 + T2
		
	end function exponentialMap
    
    function logarithmicMap (R) result (t)
		
		implicit none
		
		double precision, dimension (3, 3), intent (in) :: R
		double precision, dimension (3, 3)              :: A, A2
		double precision, dimension (3, 3)              :: Tskew
		double precision, dimension (3)                 :: t
		double precision :: Anorm
		!integer :: j
		
        A = 0.5D0 * (R - transpose (R))
        A2 = matmul (A, A)
		Anorm = sqrt (-0.5D0 * trace (A2))
        Tskew = asin (Anorm) / Anorm * A
        t = antiskew (Tskew)
		
	end function logarithmicMap
    
    pure function diagonalMatrix (a)
	
		implicit none
		
		double precision, dimension (:), intent (in) :: a
		double precision, dimension (size (a), size (a)) :: diagonalMatrix
		integer :: i
		
		diagonalMatrix = 0.0D0
        do i = 1, 6
            diagonalMatrix (i, i) = a (i)
        end do
		
	end function diagonalMatrix
    
    pure function identityMatrix (n)
	
		implicit none
		
		integer, intent (in) :: n
		double precision, dimension (n, n) :: identityMatrix
		integer :: i
		
		identityMatrix = 0.0D0
        do i = 1, n
            identityMatrix (i, i) = 1.0D0
        end do
		
	end function identityMatrix
    
    pure function trace (A)
	
		implicit none
		
		double precision, dimension (:, :), intent (in) :: A
		double precision                                :: trace
		integer :: i
		
		trace = 0.0D0
        do i = 1, size (A (1, :))
            trace = trace + A (i, i)
        end do
		
	end function trace
	
end module