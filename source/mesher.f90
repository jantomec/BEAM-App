! Mesh creator
!
! author ....... Jan Tomec
! copyright .... Copyright 2020, Project THREAD - University of Rijeka, Faculty of Civil Engineering
! credits ...... Jan Tomec, Gordan Jelenić
! license ...... GPL
! version ...... 1.0.0
! maintainer ... Jan Tomec
! email ........ jan.tomec@gradri.uniri.hr
! status ....... Development
! date ......... 06/11/2020
!
! ------------------------------------------------------------------------------
module mesh_object
	
	implicit none
	
	type ElementMesh
		integer :: Nno, Nele, order
		double precision, dimension (:, :), allocatable :: X0
		integer, dimension (:, :), allocatable :: ele
	end type ElementMesh
end module

module mesher
	
	use mesh_object
	use vector_algebra
	
	implicit none
	
	private
	
	double precision, parameter :: PI = 4 * atan (1.0)
	
	public :: lmsh, ramsh, arcmsh
	
	contains
	
	! create line mesh
	!
	! L ............ length of the beam
	! ele .......... elements array (no ele, no nodes on ele)
	! X0 ........... initial coordinates (3, no all nodes)
	!
	! modify ele, X0
	subroutine lmsh (L, mesh)
		
		double precision, intent (in) :: L
		type (ElementMesh), intent (inout) :: mesh
		integer :: i, j
		double precision :: dx
		
		allocate (mesh%X0 (3, mesh%Nno), mesh%ele (mesh%Nele, mesh%order + 1))
		
		mesh%X0 = 0.0
		dx = L / (mesh%Nno - 1)
		mesh%X0 (3, :) = (/ (j * dx, j = 0, mesh%Nno - 1) /)
		
		do i = 1, mesh%order + 1
			mesh%ele (:, i) = (/ (j, j = i, mesh%order * mesh%Nele + i - mesh%order, mesh%order) /)
		end do
		
	end subroutine lmsh
	
	! create right-angle mesh
	!
	! L ............ length of one beam
	! ele .......... elements array (no ele, no nodes on ele)
	! X0 ........... initial coordinates (3, no all nodes)
	!
	! modify ele, X0
	subroutine ramsh (L, mesh)
		
		double precision, intent (in) :: L
		type (ElementMesh), intent (inout) :: mesh
		integer :: i, j
		double precision :: dx
		
		allocate (mesh%X0 (3, mesh%Nno), mesh%ele (mesh%Nele, mesh%order + 1))
		
		mesh%X0 = 0.0
		dx = 2 * L / (mesh%Nno - 1)
		do j = 0, mesh%Nno - 1
			mesh%X0 (3, j + 1) = j * dx / dsqrt (2.0D1)
			if (j > mesh%Nno / 2) then
				mesh%X0 (1, j + 1) = L / dsqrt (2.0D1) - (j - mesh%Nno / 2) * dx / dsqrt (2.0D1)
			else
				mesh%X0 (1, j + 1) = j * dx / dsqrt (2.0D1)
			end if
		end do
		do i = 1, mesh%order + 1
			mesh%ele (:, i) = (/ (j, j = i, mesh%order * mesh%Nele + i - mesh%order, mesh%order) /)
		end do
		
	end subroutine ramsh
	
	! create arc mesh
	!
	! L ............ length of the beam
	! ele .......... elements array (no ele, no nodes on ele)
	! X0 ........... initial coordinates (3, no all nodes)
	!
	! modify ele, X0
	subroutine arcmsh (R, phi_i, phi_f, mesh, rot, gpts)
		
		double precision, intent (in) :: R, phi_i, phi_f
		type (ElementMesh), intent (inout) :: mesh
		double precision, dimension (:, :, :, :), intent (out) :: rot
		double precision, dimension (:), intent (in) :: gpts
		double precision, dimension (mesh%order + 1) :: e
		integer :: i, j
		double precision :: dphi, phi1, phi2, phig
		double precision, dimension (3) :: rotvec
		
		allocate (mesh%X0 (3, mesh%Nno), mesh%ele (mesh%Nele, mesh%order + 1))
		
		mesh%X0 = 0.0
		dphi = (phi_f - phi_i) / (mesh%Nno - 1)
		do i = 1, mesh%Nno
			mesh%X0 (3, i) = R * cos (phi_i + (i - 1) * dphi)
			mesh%X0 (1, i) = R * sin (phi_i + (i - 1) * dphi)
		end do
		
		do i = 1, mesh%order + 1
			mesh%ele (:, i) = (/ (j, j = i, mesh%order * mesh%Nele + i - mesh%order, mesh%order) /)
		end do
		
		do i = 1, mesh%Nele
			e = mesh%ele (i, :)
			phi1 = phi_i + (e (1) - 1) * dphi
			phi2 = phi_i + (e (mesh%order + 1) - 1) * dphi
			do j = 1, size (gpts)
				phig = ((gpts (j) + 1.0) / 2.0 * (phi2 - phi1) + phi1) + PI / 2
				rotvec = (/ 0.0D1, phig, 0.0D1 /)
				rot (i, j, :, :) = rv2mat (rotvec)
			end do
		end do
		
	end subroutine arcmsh
	
	! read mesh from .mesh file  -  MISSING
	!subroutine loadmesh (filename, mesh)
	!	
	!	character (len = *), intent (in) :: filename
	!	type (ElementMesh), intent (out) :: mesh
	!	integer :: i
	!	
	!end subroutine loadmesh
	
	
end module