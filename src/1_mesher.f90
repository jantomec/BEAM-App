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
! date ......... 11/06/2020
!
! ------------------------------------------------------------------------------
module mesher
	
	use mesh_objects
	use vector_algebra
    use legendre_gauss
	
	implicit none
	
	private
	
	double precision, parameter :: PI = 4 * atan (1.0D0)
	
	public :: lineMesh, arcMesh!, ramsh, loadmesh
	
	contains
	
	function lineMesh (L, noElements, elementOrder, gaussOrder, properties, pressure)
		
		type (ElementMesh)                                 :: lineMesh
        double precision,                   intent (in)    :: L
        integer,                            intent (in)    :: noElements
        integer,                            intent (in)    :: elementOrder
        integer,                            intent (in)    :: gaussOrder
        type (ElementProperties),           intent (in)    :: properties
        double precision, dimension (6, gaussOrder), optional :: pressure
		
		double precision,   dimension (3, elementOrder*noElements+1) :: coordinates
		type (LineElement), dimension (noElements)                   :: elements
		integer,            dimension (elementOrder+1)               :: nodes
        
		integer            :: i, j
		double precision   :: dx
		double precision, dimension (gaussOrder, 3,3) :: identityMatrix
        
        identityMatrix          = 0.0D0
        identityMatrix (:, 1,1) = 1.0D0
        identityMatrix (:, 2,2) = 1.0D0
        identityMatrix (:, 3,3) = 1.0D0
        
		coordinates = 0.0D0
		dx = L / (elementOrder*noElements)
		coordinates (1, :) = (/ (j * dx, j = 0, elementOrder*noElements) /)
		        
        do i = 1, noElements
            nodes = (/ (j, j = elementOrder*(i-1)+1, elementOrder*i+1) /)
            
            call elements (i)%init (nodes=nodes, properties=properties, rotationMatrix=identityMatrix, pressure=pressure)
        end do
        
        call lineMesh%init (coordinates, elements)
		
	end function lineMesh
    
    function arcMesh (radius, phi_i, phi_f, noElements, elementOrder, gaussOrder, properties, pressure)
		
		type (ElementMesh)                                 :: arcMesh
        double precision,                   intent (in)    :: radius, phi_i, phi_f
        integer,                            intent (in)    :: noElements
        integer,                            intent (in)    :: elementOrder
        integer,                            intent (in)    :: gaussOrder
        type (ElementProperties),           intent (in)    :: properties
        double precision, dimension (6, gaussOrder), optional :: pressure
		
		double precision,   dimension (3, elementOrder*noElements+1) :: coordinates
		type (LineElement), dimension (noElements)                   :: elements
		integer,            dimension (elementOrder+1)               :: nodes
        
		integer            :: i, j
		double precision   :: dphi, phi1, phi2, phig
		double precision, dimension (gaussOrder, 3,3) :: rotationMatrix
        
        double precision, parameter :: PI = 4 * atan (1.0D0)
        double precision, dimension (gaussOrder) :: pts, wgts
        double precision, dimension (3) :: rotvec
        
		coordinates = 0.0D0
		dphi = (phi_f - phi_i) / (elementOrder*noElements)
		do i = 1, elementOrder*noElements + 1
			coordinates (1, i) = radius * cos (phi_i + (i - 1) * dphi)
			coordinates (3, i) = radius * sin (phi_i + (i - 1) * dphi)
		end do
        
        call legauss (gaussOrder, pts, wgts)
        
        do i = 1, noElements
            nodes = (/ (j, j = elementOrder*(i-1)+1, elementOrder*i+1) /)
            phi1 = phi_i + (nodes (1) - 1) * dphi
			phi2 = phi_i + (nodes (elementOrder + 1) - 1) * dphi
            do j = 1, gaussOrder
				phig = ((pts (j) + 1.0D0) / 2.0D0 * (phi2 - phi1) + phi1) + PI / 2
                rotvec = (/ 0.0D0, -phig, 0.0D0 /)
				rotationMatrix (j, :, :) = rv2mat (rotvec)
                print *, phig
			end do
            call elements (i)%init (nodes=nodes, properties=properties, rotationMatrix=rotationMatrix, pressure=pressure)
        end do
        
        call arcMesh%init (coordinates, elements)
        		
	end function arcMesh
	
	! create right-angle mesh
	!
	! L ............ length of one beam
	! ele .......... elements array (no ele, no nodes on ele)
	! X0 ........... initial coordinates (3, no all nodes)
	!
	! modify ele, X0
	! subroutine ramsh (L, mesh)
		
		! double precision, intent (in) :: L
		! type (ElementMesh), intent (inout) :: mesh
		! integer :: i, j
		! double precision :: dx
		
		! allocate (mesh%X0 (3, mesh%Nno), mesh%ele (mesh%Nele, mesh%order + 1))
		
		! mesh%X0 = 0.0D0
		! dx = 2 * L / (mesh%Nno - 1)
		! do j = 0, mesh%Nno - 1
			! mesh%X0 (3, j + 1) = j * dx / dsqrt (2.0D0)
			! if (j > mesh%Nno / 2) then
				! mesh%X0 (1, j + 1) = L / dsqrt (2.0D0) - (j - mesh%Nno / 2) * dx / dsqrt (2.0D0)
			! else
				! mesh%X0 (1, j + 1) = j * dx / dsqrt (2.0D0)
			! end if
		! end do
		! do i = 1, mesh%order + 1
			! mesh%ele (:, i) = (/ (j, j = i, mesh%order * mesh%Nele + i - mesh%order, mesh%order) /)
		! end do
		
	! end subroutine ramsh
	
	! create arc mesh
	!
	! L ............ length of the beam
	! ele .......... elements array (no ele, no nodes on ele)
	! X0 ........... initial coordinates (3, no all nodes)
	!
	! modify ele, X0
	! subroutine arcmsh (radius, phi_i, phi_f, mesh, rot, gpts)
		
		! double precision, intent (in) :: radius, phi_i, phi_f
		! type (ElementMesh), intent (inout) :: mesh
		! double precision, dimension (:, :, :, :), intent (out) :: rot
		! double precision, dimension (:), intent (in) :: gpts
		! double precision, dimension (mesh%order + 1) :: e
		! integer :: i, j
		! double precision :: dphi, phi1, phi2, phig
		! double precision, dimension (3) :: rotvec
		
		! allocate (mesh%X0 (3, mesh%Nno), mesh%ele (mesh%Nele, mesh%order + 1))
		
		! mesh%X0 = 0.0D0
		! dphi = (phi_f - phi_i) / (mesh%Nno - 1)
		! do i = 1, mesh%Nno
			! mesh%X0 (3, i) = radius * cos (phi_i + (i - 1) * dphi)
			! mesh%X0 (1, i) = radius * sin (phi_i + (i - 1) * dphi)
		! end do
		
		! do i = 1, mesh%order + 1
			! mesh%ele (:, i) = (/ (j, j = i, mesh%order * mesh%Nele + i - mesh%order, mesh%order) /)
		! end do
		
		! do i = 1, mesh%Nele
			! e = mesh%ele (i, :)
			! phi1 = phi_i + (e (1) - 1) * dphi
			! phi2 = phi_i + (e (mesh%order + 1) - 1) * dphi
			! do j = 1, size (gpts)
				! phig = ((gpts (j) + 1.0D0) / 2.0D0 * (phi2 - phi1) + phi1) + PI / 2
				! rotvec = (/ 0.0D0, phig, 0.0D0 /)
				! rot (i, j, :, :) = rv2mat (rotvec)
			! end do
		! end do
		
	! end subroutine arcmsh
	
	! read mesh from .mesh file  -  MISSING
	!subroutine loadmesh (filename, mesh)
	!	
	!	character (len = *), intent (in) :: filename
	!	type (ElementMesh), intent (out) :: mesh
	!	integer :: i
	!	
	!end subroutine loadmesh
	
	
end module