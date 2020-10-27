program spaghetti
    
    use solver
    use mesher
    use mesh_objects
    use vector_algebra
    use files
    
    implicit none
    
    double precision, parameter :: PI = 4 * atan (1.0D0), h = 0.005D0, beta = 0.25D0, gamma = 0.5D0
    character (len = *), parameter :: folder = 'spaghetti-results'
    character (len = *), parameter :: fname_format = '("step", I0.3, ".dat")'
    
    type (ElementMesh) :: mesh
    type (ElementProperties) :: properties
    
    logical,          dimension (:, :), allocatable :: DOF
    double precision, dimension (:, :), allocatable :: Uload, Q, R
    double precision :: t
    integer          :: j, i, noIter
    
    ! =================================================
    ! MATERIAL AND GEOMETRIC PROPERTIES
    properties%Area             = 1.0D0
    properties%Density          = 1.0D0
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
    Q (5, mesh%NoNodes) = 8 * PI
    
    ! =================================================
    ! SET UP RESULT FILES   
    call removeFolder (folder)
    call createFolder (folder)
    
    ! =================================================
    ! FORCE CONTROL ROUTINE
    do j = 0, 2000 ! max noumber of steps
        t = t + h
        if (j > 0) then
        
            write (6, '(/, "Step", X, I3)') j
            
            call dynamic_iter (mesh=mesh, DOF6=DOF, Uload=Uload, Q=Q, R=R, h=h, beta=beta, gamma=gamma, noIter=noIter)
            
        end if
        
        do i = 1, mesh%NoElements
            mesh%Elements(i)%RotationMatrixLastConverged = mesh%Elements(i)%RotationMatrix
        end do
        
        call writeResults (mesh, R, j, folder)
        
        if (t >= 1.0D0) exit
        
    end do
    
end program spaghetti