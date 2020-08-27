! Legendre - Gauss quadrature ordinates and weights
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

module legendre_gauss
		
	implicit none
	
	private
	
	public :: legauss
	
	contains
	
	! evaluate 1st shape function in given points (gauss points)
	pure subroutine legauss (n, x, w)
		
		implicit none
		
		integer, intent (in) :: n
		double precision, dimension (n), intent (out) :: x, w
		integer :: m, i, j
		double precision :: z, z1, p1, p2, p3, pp
		double precision, parameter :: TOL = 1.0E-30_DP, PI = 4.0_DP * atan (1.0_DP)
				
		m = (n + 1) / 2
		do i = 1, m
			z = cos (PI * (i-0.25_DP) / (n+0.5_DP))
			z1 = 0.0_DP
			do while (abs (z - z1) .gt. TOL)
				p1 = 1.0_DP
				p2 = 0.0_DP
				do j = 1, n
					p3 = p2
					p2 = p1
					p1 = ((2.0_DP * j - 1.0_DP) * z * p2 - (j - 1.0_DP) * p3) / j
				end do
				pp = n * (z * p1 - p2) / (z * z - 1.0_DP)
				z1 = z
				z = z1 - p1 / pp
			end do
			x (i) = -z
			x (n+1-i) = z
			w (i) = 2.0_DP / ((1.0_DP - z ** 2) * pp ** 2)
			w (n+1-i) = w(i)
		end do
		
	end subroutine legauss
	
end module
