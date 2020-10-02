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
	
	implicit none
	
	private
	
	public :: shfun, shdfun
	
	contains
	
	! evaluate shape function in given points (gauss points)
	function shfun (e_ord, pts)
		
		implicit none
		
		integer, intent (in) :: e_ord
		double precision, dimension (:), intent (in) :: pts
		double precision, dimension (e_ord + 1, size (pts)) :: shfun
		
		if (e_ord == 1) then
			shfun (1, :) = -0.5D0 * pts + 0.5D0
			shfun (2, :) = 0.5D0 * pts + 0.5D0
		else if (e_ord == 2) then
			shfun (1, :) = 0.5D0 * pts ** 2 - 0.5D0 * pts
			shfun (2, :) = -1.0D0 * pts ** 2 + 1.0D0
			shfun (3, :) = 0.5D0 * pts ** 2 + 0.5D0 * pts
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
		double precision, dimension (:), intent (in) :: pts
		double precision, dimension (e_ord + 1, size (pts)) :: shdfun
				
		if (e_ord == 1) then
			shdfun (1, :) = -0.5D0
			shdfun (2, :) = 0.5D0
		else if (e_ord == 2) then
			shdfun (1, :) = pts - 0.5D0
			shdfun (2, :) = -2.0D0 * pts
			shdfun (3, :) = pts + 0.5D0
		else
			print *, 'Module: shape_functions, Function: shdfun'
			print *, 'Message: Unexpected element_order.'
			stop
		end if		
	
	end function shdfun
	
end module
