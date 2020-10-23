! =============================================================================
!
! This is example 7.3 from Simo & Vu-Quoc [1985] - 
!   A three-dimensional finite-strain rod model. 
!   Part II: Computational aspects
!
! =============================================================================
!
!                  | F
!                  |
!                  v
!              *********                 EA      = 17000
!          ***           ***             GAy     = 17000
!        **                 **           GAz     = 17000
!      **                     **         Jx      = 170
!     **                       **        EIyy    = 170
!    **                         **       EIzz    = 170
!    *                           *       N       = 40
!    0                           -       L       = 1.5
!
! =============================================================================

program arch
    
    use solver
    use mesher
    use mesh_objects
    use vector_algebra
    use files
    
    implicit none
    
    double precision, parameter :: PI = 4 * atan (1.0D0), DEG = PI / 180
    double precision, parameter :: dS = 0.2D0
    integer,          parameter :: MAXITER = 20
    double precision, parameter :: TOLER = 1D-8
    integer,          parameter :: MAXSTEPS = 1000
    character (len = 12), parameter :: folder = 'arch-results'
    
    type (ElementMesh) :: mesh
    type (ElementProperties) :: properties
    double precision :: lambda
    double precision, dimension (:),    allocatable :: Uinc
    logical,          dimension (:,:),  allocatable :: DOF
    double precision, dimension (:, :), allocatable :: Q, QC, R
    integer :: j, noIter, info
    
    ! =================================================
    ! MATERIAL AND GEOMETRIC PROPERTIES
    properties%Area             = 100.0D0
    ! properties%Density          = 0.0D0
    properties%ElasticModulus   = 1.7D2
    properties%ShearModulus     = 1.7D2
    properties%InertiaPrimary   = 1.0D0
    properties%InertiaSecondary = 1.0D0
    properties%InertiaTorsion   = 1.0D0
    properties%ShearCoefficient = 1.0D0
    
    ! =================================================
    ! MESH
    mesh = arcMesh (radius=1.5D0, phi_i=-1.5D1*DEG, phi_f=1.95D2*DEG, &
                    noElements=40, elementOrder=1, gaussOrder=1, properties=properties)
                    
    ! =================================================
    ! DATA INITIALIZATION
    allocate (Uinc (3 * mesh%NoNodes))
    allocate (DOF (6, mesh%NoNodes))
    allocate (Q (6, mesh%NoNodes), QC (6, mesh%NoNodes), R (6, mesh%NoNodes))
    
    lambda = 0.0D0
    Uinc = 0.0D0
    Q = 0.0D0
    QC = 0.0D0
    R = 0.0D0
    DOF = .TRUE.
    
    ! =================================================
    ! BOUNDARY CONDITIONS
    DOF (:, 1) = .FALSE.
    DOF (:, mesh%NoNodes) = (/ .FALSE., .FALSE., .FALSE., .TRUE., .TRUE., .TRUE.  /)
    Q (3, mesh%NoNodes / 2 + 1) = -1.4D3
    
    ! =================================================
    ! SET UP RESULT FILES   
    call removeFolder (folder)
    call createFolder (folder)
    
    ! =================================================
    ! ARC-LENGTH ROUTINE
    do j = 0, MAXSTEPS
    
        if (j > 0) then
            write (6, '(/, "Step", X, I3)') j
            call arcLength_iter (mesh, DOF, Uinc, Q, QC, R, noIter, lambda, dS, TOLER, MAXITER, info)
            if (info .EQ. 1) then
                print *, "dS too large."
                stop
            end if
        end if
        
        call writeResults (mesh, R, lambda, j, folder)
        
        if (abs (lambda) > 1) exit
        
    end do
    
    if (j .eq. MAXSTEPS + 1) write (6, '(/, "Solution not converging fast enough")')
    
end program arch