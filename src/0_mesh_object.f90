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
    
    public :: ElementProperties, LineElement, ElementMesh
	
    type ElementProperties
        
        double precision :: Area
        double precision :: Density
        double precision :: ElasticModulus
        double precision :: ShearModulus
        double precision :: InertiaPrimary
        double precision :: InertiaSecondary
        double precision :: InertiaTorsion
        double precision :: ShearCoefficient
        
	end type ElementProperties
    
	type LineElement
        
        integer :: NoNodes, NoGauss
        integer,          dimension (:),       allocatable :: Nodes
		double precision, dimension (6, 6)                 :: C
		double precision, dimension (6, 6)                 :: InertiaMatrix  ! (NoGauss, 3, 3)
		double precision, dimension (:, :, :), allocatable :: RotationMatrix  ! (NoGauss, 3, 3)
		double precision, dimension (:, :),    allocatable :: Strain, Stress, Pressure  ! (6, NoGauss)
        
        double precision :: A    ! Area
        double precision :: rho  ! Density
        double precision :: E    ! ElasticModulus
        double precision :: G    ! ShearModulus
        double precision :: I1   ! Inertia around the primary main axis
        double precision :: I2   ! Inertia around the secondary main axis
        double precision :: Ip   ! Inertia around the secondary main axis
        double precision :: It   ! Torsional inertia
        double precision :: k    ! ShearCoefficient
        
        contains
        
        procedure :: init => LineElement_init
        
	end type LineElement
    
    type ElementMesh
    
		integer :: NoNodes, NoElements
		double precision,   dimension (:, :), allocatable :: Coordinates
		double precision,   dimension (:, :), allocatable :: Displacements
		double precision,   dimension (:, :), allocatable :: Positions
		type (LineElement), dimension (:),    allocatable :: Elements
        
        contains
        
        procedure :: init => ElementMesh_init
        
	end type ElementMesh
    
    contains
    
    subroutine LineElement_init (self, nodes, properties, rotationMatrix, strain, stress, pressure)
        
        implicit none
        
        class (LineElement), intent (inout)             :: self
        
        integer,          dimension (:)                 :: nodes
		type (ElementProperties)                        :: properties
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
        self%RotationMatrix = rotationMatrix
        
        self%Strain = 0.0D0
        self%Stress = 0.0D0
        self%Pressure = 0.0D0
        
        if (present (strain)) self%Strain = strain
        if (present (stress)) self%Stress = stress
        if (present (pressure)) self%Pressure = pressure
        
        self%A   = properties%Area
        self%rho = properties%Density
        self%E   = properties%ElasticModulus
        self%G   = properties%ShearModulus
        self%I1  = properties%InertiaPrimary
        self%I2  = properties%InertiaSecondary
        self%Ip  = self%I1 + self%I2
        self%It  = properties%InertiaTorsion
        self%k   = properties%ShearCoefficient
        
        self%C        = 0.0D0
        self%C (1, 1) = self%G * self%A * self%k
        self%C (2, 2) = self%G * self%A * self%k
        self%C (3, 3) = self%E * self%A
        self%C (4, 4) = self%E * self%I1
        self%C (5, 5) = self%E * self%I2
        self%C (6, 6) = self%G * self%It
        
        !self%InertiaMatrix
        
    end subroutine
    
	subroutine ElementMesh_init (self, coordinates, elements)
        
        implicit none
        
        class (ElementMesh), intent (inout)             :: self
        
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
        
    end subroutine
    
end module