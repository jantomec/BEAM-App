program cantilever
    
    use solver
    use mesher
    use mesh_objects
    use vector_algebra
    use files
    
    implicit none
    
    integer,          parameter :: noSteps = 1
    double precision, parameter :: PI = 4 * atan (1.0D0), Q0 = 8 * PI
    character (len = *), parameter :: folder = 'cantilever-results'
    character (len = *), parameter :: fname_format = '("step", I0.3, ".dat")'
    
    type (ElementMesh) :: mesh
    type (ElementProperties) :: properties
    
    logical,          dimension (:, :), allocatable :: DOF
    double precision, dimension (:, :), allocatable :: Uload, Q, R
    integer :: j, noIter
    
    ! =================================================
    ! MATERIAL AND GEOMETRIC PROPERTIES
    properties%Area             = 1.0D0
    ! properties%Density          = 0.0D0
    properties%ElasticModulus   = 1.0D0
    properties%ShearModulus     = 1.0D0
    properties%InertiaPrimary   = 2.0D0
    properties%InertiaSecondary = 1.0D0
    properties%InertiaTorsion   = 1.0D0
    properties%ShearCoefficient = 1.0D0
        
    ! =================================================
    ! MESH  
    mesh = lineMesh (L=1.0D0, noElements=5, elementOrder=1, gaussOrder=1, properties=properties)
    
    ! =================================================
    ! DATA INITIALIZATION
    allocate (DOF (6, mesh%NoNodes))
    allocate (Uload (6, mesh%NoNodes), Q (6, mesh%NoNodes), R (6, mesh%NoNodes))
    
    Uload = 0.0D0
    Q = 0.0D0
    R = 0.0D0
    DOF = .TRUE.
    
    ! =================================================
    ! BOUNDARY CONDITIONS
    DOF (:, 1) = .FALSE.
    
    ! =================================================
    ! SET UP RESULT FILES   
    call removeFolder (folder)
    call createFolder (folder)
    
    ! =================================================
    ! FORCE CONTROL ROUTINE
    do j = 0, noSteps
        if (j > 0) then
        
            write (6, '(/, "Step", X, I3)') j
            
            Q (5, mesh%NoNodes) = Q (5, mesh%NoNodes) + Q0 / noSteps
            
            call newton_iter (mesh=mesh, DOF6=DOF, Uload=Uload, Q=Q, R=R, noIter=noIter)
            
        end if
        
        call writeResults (mesh, R, j, folder)
        
    end do
    
end program cantilever