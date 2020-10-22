program cantilever
    
    use solver
    use mesher
    use mesh_objects
    use vector_algebra
    use files
    
    implicit none
    
    integer,          parameter :: noElements = 5
    integer,          parameter :: elementOrder = 1
    integer,          parameter :: gaussOrder = 1
    integer,          parameter :: noSteps = 1
    integer,          parameter :: noNodes = noElements*elementOrder + 1
    double precision, parameter :: L = 1.0D0
    double precision, parameter :: PI = 4 * atan (1.0D0), Q0 = 8 * PI
    double precision, dimension (6), parameter :: material = (/ 1.0D0, 1.0D0, 1.0D0, 1.0D0, 2.0D0, 1.0D0 /)
    character (len = *), parameter :: folder = 'cantilever-results'
    character (len = *), parameter :: fname_format = '("step", I0.3, ".dat")'
    
    type (ElementMesh) :: mesh
    
    double precision, dimension (6, 6) :: C
    logical, dimension (6, noNodes) :: DOF
    double precision, dimension (6, noNodes) :: Uload, Q, R
    integer :: j, noIter
    
    ! =================================================
    ! ELASTIC MODULI MATRIX
    C = diagonalMatrix (material)
    
    ! =================================================
    ! MESH  
    mesh = lineMesh (L=L, noElements=noElements, elementOrder=elementOrder, gaussOrder=gaussOrder, C=C)
    
    ! =================================================
    ! SET UP RESULT FILES   
    call removeFolder (folder)
    call createFolder (folder)
    
    ! =================================================
    ! DATA INITIALIZATION
    Uload = 0.0D0
    Q = 0.0D0
    R = 0.0D0
    
    ! =================================================
    ! BOUNDARY CONDITIONS
    DOF = .TRUE.
    DOF (:, 1) = .FALSE.
    
    ! =================================================
    ! FORCE CONTROL ROUTINE
    do j = 0, noSteps
        if (j > 0) then
        
            write (6, '(/, "Step", X, I3)') j
            
            Q (5, noNodes) = Q (5, noNodes) + Q0 / noSteps
            
            call newton_iter (mesh=mesh, DOF6=DOF, Uload=Uload, Q=Q, R=R, noIter=noIter)
            
        end if
        
        call writeResults (mesh, R, j, folder)
        
    end do
    
end program cantilever