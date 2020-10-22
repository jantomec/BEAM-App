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
!           _ _ _ _v_ _ _ _                 EA      = 170000
!         /                 \               GAy     = 170
!        /                   \              GAz     = 17000
!       /                     \             Jx      = 170
!      /                       \            EIyy    = 170
!     /                         \           EIzz    = 170
!    |                           |          N       = 40
!    |                           |          L       = 1.5
!    |                           |
!    0                           -
!
! =============================================================================

program arch
    
    use solver
    use mesher
    use mesh_objects
    use vector_algebra
    use files
    
    implicit none
    
    integer, parameter :: noElements = 40
    integer, parameter :: elementOrder = 1
    integer, parameter :: gaussOrder = 1
    double precision, parameter :: L = 1.5D0
    double precision, parameter :: PI = 4 * atan (1.0D0), DEG = PI / 180
    double precision, parameter :: Q0 = -1400.0D0
    double precision, parameter :: dS = 0.2D0
    integer, parameter :: MAXITER = 20
    double precision, parameter :: TOLER = 1D-8
    integer, parameter :: MAXSTEPS = 1000
    double precision, dimension (6), parameter :: material = (/ 170000.0D0, 170.0D0, 17000.0D0, 170.0D0, 170.0D0, 170.0D0 /)
    character (len = 12), parameter :: folder = 'arch-results'
    
    type (ElementMesh) :: mesh
    double precision :: lambda
    integer, parameter :: noNodes = noElements * elementOrder + 1
    double precision, dimension (3 * noNodes) :: Uinc
    double precision, dimension (6, 6) :: C
    logical, dimension (6, noNodes) :: DOF
    double precision, dimension (6, noNodes) :: Q, QC, R
    integer :: i, j, noIter, info
    
    ! =================================================
    ! DATA INITIALIZATION
    lambda = 0.0D0
    Uinc = 0.0D0
    Q = 0.0D0
    QC = 0.0D0
    R = 0.0D0
    
    ! =================================================
    ! ELASTIC MODULI MATRIX
    C = diagonalMatrix (material)
    
    ! =================================================
    ! MESH
    mesh = arcMesh (L, -15.0D0 * DEG, 195.0D0 * DEG, noElements, elementOrder, gaussOrder, C)
    
    ! =================================================
    ! SET UP RESULT FILES   
    call removeFolder (folder)
    call createFolder (folder)
    
    ! =================================================
    ! BOUNDARY CONDITIONS
    DOF = .TRUE.
    DOF (:, 1) = .FALSE.
    DOF (:, noNodes) = (/ .FALSE., .FALSE., .FALSE., .TRUE., .TRUE., .TRUE.  /)
    Q (1, noNodes / 2 + 1) = Q0
    
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
        
        call writeResults (mesh, R, j, folder) ! missing lambda
        
        ! write (fname, fname_format) j
        ! fullname = trim (folder)//'/'//trim (fname)
        ! open (unit = 12, file = trim (fullname), status = 'unknown', action = 'write')
        ! write (12, '("node", 19X, "X1", 19X, "X2", 19X, "X3", 19X, "U1", 19X, "U2", 19X, "U3", 19X, "R1", 19X, &
            ! "R2", 19X, "R3", 19X, "M1", 19X, "M2", 19X, "M3")')
        ! do i = 1, noNodes
            ! write (12, '(I0.4, X, ES20.13, X, ES20.13, X, ES20.13, X, ES20.13, X, ES20.13, X, &
            ! ES20.13, X, ES20.13, X, ES20.13, X, ES20.13, X, ES20.13, X, ES20.13, X, &
            ! ES20.13)') i, X (1, i), X (2, i), X (3, i), U (1, i), U (2, i), U (3, i), R (1, i), &
            ! R (2, i), R (3, i), R (4, i), R (5, i), R (6, i)
        ! end do
        ! close(12)
        
        if (abs (lambda) > 1) exit
        
    end do
    
    if (j .eq. MAXSTEPS + 1) write (6, '(/, "Solution not converging fast enough")')
    
end program arch