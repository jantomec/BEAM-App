! Shape functions
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

module shape_functions

	use f95_precision
	
	implicit none
	
	private
	
	public :: shfun, shdfun
	
	contains
	
	! evaluate 1st shape function in given points (gauss points)
	pure function N1 (e_ord, pts)
		
		implicit none
		
		integer, intent (in) :: e_ord
		real (DP), dimension (:), intent (in) :: pts
		real (DP), dimension (size (pts)) :: N1
				
		if (e_ord == 1) then
			N1 = -0.5_DP * pts + 0.5_DP
		else if (e_ord == 2) then
			N1 = 0.5_DP * pts ** 2 - 0.5_DP * pts
		end if
		
	end function N1
	
	! evaluate 2nd shape function in given points (gauss points)
	pure function N2 (e_ord, pts)
		
		implicit none
		
		integer, intent (in) :: e_ord
		real (DP), dimension (:), intent (in) :: pts
		real (DP), dimension (size (pts)) :: N2
				
		if (e_ord == 1) then
			N2 = 0.5_DP * pts + 0.5_DP
		else if (e_ord == 2) then
			N2 = -1.0_DP * pts ** 2 + 1.0_DP
		end if
		
	end function N2
	
	! evaluate 3rd shape function in given points (gauss points)
	pure function N3 (e_ord, pts)
		
		implicit none
		
		integer, intent (in) :: e_ord
		real (DP), dimension (:), intent (in) :: pts
		real (DP), dimension (size (pts)) :: N3
				
		if (e_ord == 2) then
			N3 = 0.5_DP * pts ** 2 + 0.5_DP * pts
		end if
		
	end function N3
	
	! evaluate 1st shape function in given points (gauss points)
	pure function dN1 (e_ord, pts)
		
		implicit none
		
		integer, intent (in) :: e_ord
		real (DP), dimension (:), intent (in) :: pts
		real (DP), dimension (size (pts)) :: dN1
				
		if (e_ord == 1) then
			dN1 = -0.5_DP
		else if (e_ord == 2) then
			dN1 = pts - 0.5_DP
		end if
		
	end function dN1
	
	! evaluate 2nd shape function in given points (gauss points)
	pure function dN2 (e_ord, pts)
		
		implicit none
		
		integer, intent (in) :: e_ord
		real (DP), dimension (:), intent (in) :: pts
		real (DP), dimension (size (pts)) :: dN2
				
		if (e_ord == 1) then
			dN2 = 0.5_DP
		else if (e_ord == 2) then
			dN2 = -2.0_DP * pts
		end if
		
	end function dN2
	
	! evaluate 3rd shape function in given points (gauss points)
	pure function dN3 (e_ord, pts)
		
		implicit none
		
		integer, intent (in) :: e_ord
		real (DP), dimension (:), intent (in) :: pts
		real (DP), dimension (size (pts)) :: dN3
				
		if (e_ord == 2) then
			dN3 = pts + 0.5_DP
		end if
		
	end function dN3
	
	! evaluate shape function in given points (gauss points)
	function shfun (e_ord, pts)
		
		implicit none
		
		integer, intent (in) :: e_ord
		real (DP), dimension (:), intent (in) :: pts
		real (DP), dimension (e_ord + 1, size (pts)) :: shfun
		
		if (e_ord == 1) then
			shfun (1, :) = N1 (1, pts)
			shfun (2, :) = N2 (1, pts)
		else if (e_ord == 2) then
			shfun (1, :) = N1 (2, pts)
			shfun (2, :) = N2 (2, pts)
			shfun (3, :) = N3 (2, pts)
		else
			print *, 'Module: shape_functions, Function: shfun'
			print *, 'Message: Unexpected element_order.'
			stop
		end if		
	
	end function shfun
	
	! evaluate shape function derivatives in given points (gauss points)
	function shdfun (e_ord, pts)
		
		implicit none
		
		integer, intent (in) :: e_ord
		real (DP), dimension (:), intent (in) :: pts
		real (DP), dimension (e_ord + 1, size (pts)) :: shdfun
				
		if (e_ord == 1) then
			shdfun (1, :) = dN1 (1, pts)
			shdfun (2, :) = dN2 (1, pts)
		else if (e_ord == 2) then
			shdfun (1, :) = dN1 (2, pts)
			shdfun (2, :) = dN2 (2, pts)
			shdfun (3, :) = dN3 (2, pts)
		else
			print *, 'Module: shape_functions, Function: shdfun'
			print *, 'Message: Unexpected element_order.'
			stop
		end if		
	
	end function shdfun
	
end module
