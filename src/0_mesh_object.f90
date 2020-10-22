! Support objects
!
! author .............. Jan Tomec
! copyright ........... Copyright 2020, Project THREAD - University of Rijeka, Faculty of Civil Engineering
! credits ............. Jan Tomec, Gordan Jelenić
! license ............. GPL
! version ............. 1.0.0
! maintainer .......... Jan Tomec
! email ............... jan.tomec@gradri.uniri.hr
! status .............. Development
! date created ........ 04/06/2020
! date modified (1) ... 22/10/2020
!
! ------------------------------------------------------------------------------

module mesh_objects
	
	implicit none
    
    private
    
    public :: LineElement, ElementMesh, LineElement_init, ElementMesh_init
	
	type LineElement
        
        integer :: NoNodes, NoGauss
        integer,          dimension (:),       allocatable :: Nodes
		double precision, dimension (6, 6)                 :: C
		double precision, dimension (:, :, :), allocatable :: RotationMatrix  ! (NoGauss, 3, 3)
		double precision, dimension (:, :),    allocatable :: Strain, Stress, Pressure  ! (6, NoGauss)
                
	end type LineElement
    
    type ElementMesh
    
		integer :: NoNodes, NoElements
		double precision,   dimension (:, :), allocatable :: Coordinates
		double precision,   dimension (:, :), allocatable :: Displacements
		double precision,   dimension (:, :), allocatable :: Positions
		type (LineElement), dimension (:),    allocatable :: Elements    
        
	end type ElementMesh
    
    contains
    
    function LineElement_init (nodes, C, rotationMatrix, strain, stress, pressure) result (self)
        
        implicit none
        
        type (LineElement)                              :: self
        
        integer,          dimension (:)                 :: nodes
		double precision, dimension (6, 6)              :: C
		double precision, dimension (:, :, :)           :: rotationMatrix
		double precision, dimension (:, :), optional    :: strain, stress, pressure
        
        self%NoNodes = size (nodes)
        self%NoGauss = size (rotationMatrix(:, 1, 1))
        
        allocate (self%Nodes (self%NoNodes))
        allocate (self%RotationMatrix (self%NoGauss, 3, 3))
        allocate (self%Strain (6, self%NoGauss))
        allocate (self%Stress (6, self%NoGauss))
        allocate (self%Pressure (6, self%NoGauss))
        
        self%nodes = nodes
        self%C = C
        self%RotationMatrix = rotationMatrix
        
        self%Strain = 0.0D0
        self%Stress = 0.0D0
        self%Pressure = 0.0D0
        
        if (present (strain)) self%Strain = strain
        if (present (stress)) self%Stress = stress
        if (present (pressure)) self%Pressure = pressure
        
    end function
    
	function ElementMesh_init (coordinates, elements) result (self)
        
        implicit none
        
        type (ElementMesh)                              :: self
        
        double precision,  dimension (:, :)             :: coordinates
		type(LineElement), dimension (:)                :: elements
        
        self%NoNodes = size (coordinates(1, :))
        self%NoElements = size (elements)
        
        allocate (self%Coordinates (3, self%NoNodes))
        allocate (self%Elements (self%NoElements))
        
        self%Coordinates = coordinates
        self%Elements = elements
        self%Displacements = 0.0D0 * coordinates
        self%Positions = coordinates
        
    end function
    
end module