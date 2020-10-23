! Finite element solver
!
! author .............. Jan Tomec
! Valueright ........... Valueright 2020, Project THREAD - University of Rijeka, Faculty of Civil Engineering
! credits ............. Jan Tomec, Gordan Jelenić
! license ............. GPL
! version ............. 1.0.0
! maintainer .......... Jan Tomec
! email ............... jan.tomec@gradri.uniri.hr
! status .............. Development
! date created ........ 06/09/2020
! date modified (1) ... 14/09/2020
! date modified (2) ... 02/10/2020
!
! ------------------------------------------------------------------------------

module solver

    use beam
    use mesh_objects
    
    implicit none
    
    private
    
    public :: newton_iter, arcLength_iter, dynamic_iter!, newton_iter_det, arc_length_iter
    
    contains
        
    subroutine begin_table (typ)
        
        implicit none
        
        character :: typ
        
        if (typ .eq. 'A') then      
            write (6, '(/, 3X, "n", X, "|", X, "residual vector norm", X, "|", X, "load parameter")')
            write (6, '(A)') repeat("-", 45)
        else
            write (6, '(/, 3X, "n", X, "|", X, "residual vector norm")')
            write (6, '(A)') repeat("-", 28)
        end if
        
    end subroutine begin_table
    
    subroutine solve (A, b, dof)
    
        implicit none
        
        double precision, dimension (:, :)                      :: A     ! (6*no nodes, 6*no nodes)
        double precision, dimension (:, :)    , intent (inout)  :: b     ! (6, no nodes)
        logical         , dimension (:)       , intent (in)     :: dof   ! (6*no nodes)
        
        integer :: ndof, i, errck
        integer, dimension (size (dof)) :: ipiv
        double precision, dimension (size (dof)) :: bflat
        
        external dgetrf, dgetrs
        
        ndof = size (dof)
        bflat = pack (b, .TRUE.)
        
        do i = 1, ndof
            if (.NOT. dof (i)) then
                A (:, i) = 0.0D0
                A (i, :) = 0.0D0
                A (i, i) = 1.0D0
            end if
        end do
        
        call dgetrf (ndof, ndof, A, ndof, ipiv, errck)  ! LU factorization
        if (errck .ne. 0) then
            print *, 'Module: solver, Function: dgetrf'
            print *, 'Message: LU factorization error.'
            stop
        end if
        
        call dgetrs ('N', ndof, 1, A, ndof, ipiv, bflat, ndof, errck)  ! solve system
        if (errck .ne. 0) then
            print *, 'Module: solver, Function: dgetrs'
            print *, 'Message: Invert matrix error.'
            stop
        end if
        
        b = reshape (bflat, (/ 6, size (b (1,:)) /))
        
    end subroutine solve
     
    subroutine solve2 (A, b, dof)
    
        implicit none
        
        double precision, dimension (:, :)                  :: A     ! (6*no nodes, 6*no nodes)
        double precision, dimension (:, :), intent (inout)  :: b     ! (6*no nodes, 2)
        logical,          dimension (:),    intent (in)     :: dof   ! (6*no nodes)
        
        integer :: ndof, i, errck
        integer, dimension (size (dof)) :: ipiv
        
        external dgetrf, dgetrs
        
        ndof = size (dof)
                
        do i = 1, ndof
            if (.NOT. dof (i)) then
                A (:, i) = 0.0D0
                A (i, :) = 0.0D0
                A (i, i) = 1.0D0
            end if
        end do
        
        call dgetrf (ndof, ndof, A, ndof, ipiv, errck)  ! LU factorization
        if (errck .ne. 0) then
            print *, 'Module: solver, Function: dgetrf'
            print *, 'Message: LU factorization error.'
            stop
        end if
                
        call dgetrs ('N', ndof, 2, A, ndof, ipiv, b, ndof, errck)  ! solve system
        if (errck .ne. 0) then
            print *, 'Module: solver, Function: dgetrs'
            print *, 'Message: Invert matrix error.'
            stop
        end if
        
    end subroutine solve2
    
    ! perform newton-raphson iteration
    !
    ! ele .......... elements array (no ele, no nodes on ele)
    ! X0 ........... initial coordinates (3, no all nodes)
    ! U ............ displacement (3, no all nodes)
    ! C ............ tangent elastic moduli (6, 6)
    ! DOF6 ......... degrees of freedom (6, no all nodes)
    ! Uload ........ prescribed displacement/rotation (6, no all nodes)
    ! Q ............ nodal force load (6, no all nodes)
    ! pressure ..... distributed load (no ele, 6, no gauss)
    ! rot .......... rotations in gauss points (no ele, no gauss, 3, 3)
    ! om ........... curvature vector (no ele, 3, no gauss)
    ! stress ....... internal force vector (no ele, 3, no gauss)
    ! TOLER ........ tolerance for convergence
    ! MAXITER ...... maximum no of iterations
    ! TEST ......... convergence test ('RSD' - resdidual, 'DSP' - displacement)
    ! noIter ........ integer, number of iterations used to converge
    ! prints ....... boolean, indicating if intermediate errck statements should be printed
    !
    ! modify U, rot, om, stress
    subroutine newton_iter (mesh, DOF6, Uload, Q, R, noIter, TOLER, MAXITER, TEST, PRINTS)
    
        implicit none
        
        type (ElementMesh),                            intent (inout) :: mesh
        logical,          dimension (6, mesh%NoNodes), intent (in)    :: DOF6
        double precision, dimension (6, mesh%NoNodes), intent (in)    :: Uload, Q
        double precision, dimension (6, mesh%NoNodes), intent (out)   :: R
        integer,                                       intent (out)   :: noIter
        
        double precision,    optional :: TOLER
        integer,             optional :: MAXITER
        character (len = 3), optional :: TEST
        logical,             optional :: PRINTS
        double precision              :: TOLERValue
        integer                       :: MAXITERValue
        character (len = 3)           :: TESTValue
        logical                       :: PRINTSValue
        
        double precision, dimension (6 * mesh%NoNodes, 6 * mesh%NoNodes) :: tangent
        double precision, dimension (3, mesh%NoNodes)                    :: X, dU, dth
        double precision, dimension (6, mesh%NoNodes)                    :: Fint, Fext
        logical,          dimension (6 * mesh%NoNodes)                   :: dof
        double precision                                                 :: convtest
        integer                                                          :: i, j, k, errck
        
        if (      present(TOLER))   TOLERValue   = TOLER
        if (.NOT. present(TOLER))   TOLERValue   = 1.0D-8
        if (      present(MAXITER)) MAXITERValue = MAXITER
        if (.NOT. present(MAXITER)) MAXITERValue = 30
        if (      present(TEST))    TESTValue    = TEST
        if (.NOT. present(TEST))    TESTValue    = 'RSD'
        if (      present(prints))  PRINTSValue  = PRINTS
        if (.NOT. present(prints))  PRINTSValue  = .TRUE.
        
        R = Uload  ! fill R with displacement load
        dof = pack (DOF6, .TRUE.)
        
        if (PRINTSValue) call begin_table ('N')
        
        do i = 0, MAXITERValue-1
            
            noIter = i + 1
            
            if (i > 0) then
                tangent = assemble_tangent (mesh)  ! tangent
                
                call solve (tangent, R, dof)  ! solve the system, fill R with results
                
            end if
            
            dU = R (1:3, :)  ! divide displacements and rotations
            dth = R (4:6, :)
            
            mesh%Displacements = mesh%Displacements + dU
            mesh%Positions = mesh%Coordinates + mesh%Displacements
                        
            if (TESTValue .eq. 'DSP') convtest = norm2 (R)  ! R are incremental updates
            
            call update_stress_strain (mesh, dth)
                        
            Fint = assemble_internal_force (mesh)
            Fext = assemble_external_force (mesh, Q)
            R = Fext - Fint  ! compute residual, fill R with residual forces
            
            do j = 1, 6
                do k = 1, mesh%NoNodes
                    if (.NOT. DOF6 (j, k)) R(j, k) = 0.0D0
                end do
            end do
            
            if (TESTValue .eq. 'RSD') convtest = norm2 (R)  ! R are residual forces
                        
            if (PRINTSValue) write (6, '(I4, X, "|", X, ES20.13)') i, convtest
            if (convtest < TOLERValue) exit
            
        end do
                
        if (i .eq. MAXITERValue) then
            if (PRINTSValue) write (6, '(/, "Not converging")')
            stop
        end if
        
    end subroutine newton_iter
    
    subroutine dynamic_iter (mesh, DOF6, Uload, Q, R, h, beta, gamma, noIter, TOLER, MAXITER, TEST, PRINTS)
    
        implicit none
        
        type (ElementMesh),                            intent (inout) :: mesh
        logical,          dimension (6, mesh%NoNodes), intent (in)    :: DOF6
        double precision, dimension (6, mesh%NoNodes), intent (in)    :: Uload, Q
        double precision, dimension (6, mesh%NoNodes), intent (out)   :: R
        double precision,                              intent (in)    :: h, beta, gamma
        integer,                                       intent (out)   :: noIter
        
        double precision,    optional :: TOLER
        integer,             optional :: MAXITER
        character (len = 3), optional :: TEST
        logical,             optional :: PRINTS
        double precision              :: TOLERValue
        integer                       :: MAXITERValue
        character (len = 3)           :: TESTValue
        logical                       :: PRINTSValue
        
        double precision, dimension (6 * mesh%NoNodes, 6 * mesh%NoNodes) :: tangent
        double precision, dimension (3, mesh%NoNodes)                    :: X, dU, dth
        double precision, dimension (6, mesh%NoNodes)                    :: Fint, Fext, Fine
        logical,          dimension (6 * mesh%NoNodes)                   :: dof
        double precision                                                 :: convtest
        integer                                                          :: i, j, k, errck
        
        if (      present(TOLER))   TOLERValue   = TOLER
        if (.NOT. present(TOLER))   TOLERValue   = 1.0D-8
        if (      present(MAXITER)) MAXITERValue = MAXITER
        if (.NOT. present(MAXITER)) MAXITERValue = 30
        if (      present(TEST))    TESTValue    = TEST
        if (.NOT. present(TEST))    TESTValue    = 'RSD'
        if (      present(prints))  PRINTSValue  = PRINTS
        if (.NOT. present(prints))  PRINTSValue  = .TRUE.
        
        R = Uload  ! fill R with displacement load
        dof = pack (DOF6, .TRUE.)
        
        if (PRINTSValue) call begin_table ('N')
        
        do i = 0, MAXITERValue-1
            
            noIter = i + 1
            
            if (i > 0) then
                tangent = assemble_tangent_dynamic (mesh, h, beta, gamma)  ! tangent
                
                call solve (tangent, R, dof)  ! solve the system, fill R with results
                
            end if
            
            dU = R (1:3, :)  ! divide displacements and rotations
            dth = R (4:6, :)
            
            mesh%Displacements = mesh%Displacements + dU
            mesh%Velocities = mesh%Velocities + gamma / (h*beta) * dU
            mesh%Accelerations = mesh%Accelerations + gamma / (h**2*beta) * dU
            mesh%Positions = mesh%Coordinates + mesh%Displacements
                        
            if (TESTValue .eq. 'DSP') convtest = norm2 (R)  ! R are incremental updates
            
            call update_stress_strain (mesh, dth)
            call updateDynamics (mesh, dth, h, beta, gamma)
                        
            Fint = assemble_internal_force (mesh)
            Fext = assemble_external_force (mesh, Q)
            Fine = assemble_inertial_force (mesh)
            R = Fine + Fext - Fint  ! compute residual, fill R with residual forces
            
            do j = 1, 6
                do k = 1, mesh%NoNodes
                    if (.NOT. DOF6 (j, k)) R (j, k) = 0.0D0
                end do
            end do
            
            if (TESTValue .eq. 'RSD') convtest = norm2 (R)  ! R are residual forces
                        
            if (PRINTSValue) write (6, '(I4, X, "|", X, ES20.13)') i, convtest
            if (convtest < TOLERValue) exit
            
        end do
                
        if (i .eq. MAXITERValue) then
            if (PRINTSValue) write (6, '(/, "Not converging")')
            stop
        end if
        
    end subroutine dynamic_iter
    
    recursive function det_rosetta ( mat, n ) result( accum )
        integer :: n
        double precision  :: mat(n, n)
        double precision    :: submat(n-1, n-1), accum
        integer :: i, sgn

        if ( n == 1 ) then
            accum = mat(1,1)
        else
            accum = 0.0D0
            sgn = 1
            do i = 1, n
                submat( 1:n-1, 1:i-1 ) = mat( 2:n, 1:i-1 )
                submat( 1:n-1, i:n-1 ) = mat( 2:n, i+1:n )

                accum = accum + sgn * mat(1, i) * det_rosetta( submat, n-1 )
                sgn = - sgn
            enddo
        endif
    end function
        
    ! perform arc-length iteration
    !
    ! ele .......... elements array (no ele, no nodes on ele)
    ! X0 ........... initial coordinates (3, no all nodes)
    ! U ............ displacement (3, no all nodes)
    ! C ............ tangent elastic moduli (6, 6)
    ! DOF6 .......... degrees of freedom (6, no all nodes)
    ! Uload ........ prescribed displacement/rotation (6, no all nodes)
    ! Q ............ nodal load (6, no all nodes)
    ! pressure ............ distributed load (no ele, 6, no gauss)
    ! rot .......... rotations in gauss points (no ele, no gauss, 3, 3)
    ! om ........... curvature vector (no ele, 3, no gauss)
    ! stress ............ internal force vector (no ele, 3, no gauss)
    ! resout ....... residual vector (6, no all nodes)
    ! lambda ....... load multiplier
    ! dS ........... convergence radius
    ! TOLER ........ tolerance for convergence
    ! MAXITER ...... maximum no of iterations
    ! TEST ......... convergence test ('RSD' - resdidual, 'DSP' - displacement)
    !
    ! modify U, rot, om, stress
    subroutine arcLength_iter (mesh, DOF6, Uinc, Q, QC, R, noIter, lambda, dS, TOLER, MAXITER, errck)
    
        implicit none
        
        type (ElementMesh),                            intent (inout) :: mesh
        logical,          dimension (6, mesh%NoNodes), intent (in)    :: DOF6
        double precision, dimension (6, mesh%NoNodes), intent (in)    :: Q, QC
        double precision, dimension (6, mesh%NoNodes), intent (out)   :: R
        integer,                                       intent (out)   :: noIter
        double precision,                              intent (inout) :: lambda
        double precision,                              intent (in)    :: dS
        
        double precision,    optional :: TOLER
        integer,             optional :: MAXITER
        double precision              :: TOLERValue
        integer                       :: MAXITERValue
        
        double precision, dimension (6 * mesh%NoNodes, 6 * mesh%NoNodes) :: tangent
        double precision, dimension (3, mesh%NoNodes)                    :: X, dU, dth, dUF, dUR, thF, thR
        double precision, dimension (6 * mesh%NoNodes)                   :: resFflat, resRflat
        double precision, dimension (6 * mesh%NoNodes, 2)                :: R2
        double precision, dimension (6, mesh%NoNodes)                    :: resF, resR
        double precision, dimension (3 * mesh%NoNodes)                   :: Uinc, dUFflat, dURflat
        double precision, dimension (6, mesh%NoNodes)                    :: Fint, Fext, FC
        logical,          dimension (6 * mesh%NoNodes)                   :: dof
        double precision :: convtest, dlambda, dlambda1, dlambda2, discriminant, UincdotdUF, a1, a2, a3
        integer          :: i, j, errck
        double precision, dimension (2) :: dlambda_test
        integer, dimension (1) :: dlambda_test_res
        
        if (      present(TOLER))   TOLERValue   = TOLER
        if (.NOT. present(TOLER))   TOLERValue   = 1.0D-8
        if (      present(MAXITER)) MAXITERValue = MAXITER
        if (.NOT. present(MAXITER)) MAXITERValue = 30
        
        errck = 0
                        
        dof = pack (DOF6, .TRUE.)
        
        Fint = assemble_internal_force (mesh)
        Fext = assemble_external_force (mesh, Q)
        FC = assemble_external_force (mesh, QC)
        
        R = FC + lambda * Fext - Fint
        R2 (:, 1) = pack (Fext, .TRUE.)
        R2 (:, 2) = pack (R, .TRUE.)
        
        do j = 1, 6*mesh%NoNodes
            if (.NOT. dof (j)) R2 (j, :) = 0.0D0
        end do
        
        call begin_table ('A')
        write (6, '(I4, X, "|", X, ES20.13, X, "|", X, F14.10)') 0, norm2 (R), lambda
        
        do i = 1, MAXITERValue
            
            noIter = i
                    
            tangent = assemble_tangent (mesh)
            
            call solve2 (tangent, R2, dof)
                        
            resFflat = R2 (:, 1)
            resRflat = R2 (:, 2)
            
            resF = reshape (resFflat, (/ 6, mesh%NoNodes /))
            resR = reshape (resRflat, (/ 6, mesh%NoNodes /))
            
            do concurrent (j = 1:3)
                dUF (j, :) = resF (j, :)
                thF (j, :) = resF (j + 3, :)
                dUR (j, :) = resR (j, :)
                thR (j, :) = resR (j + 3, :)
            end do
            
            dUFflat = pack (dUF, .TRUE.)
            dURflat = pack (dUR, .TRUE.)
            
            a1 = dot_product (dUFflat, dUFflat)
            
            if (i .eq. 1) then
                UincdotdUF = dot_product (Uinc, dUFflat)
                dlambda = sign (dS / sqrt (a1), UincdotdUF)
                Uinc = 0.0D0
            else
                a2 = 2 * dot_product(Uinc + dURflat, dUFflat)
                a3 = dot_product (Uinc + dURflat, Uinc + dURflat) - dS ** 2
                
                if (a2 ** 2 - 4 * a1 * a3 < 0) then
                    errck = 1
                    exit
                end if
                
                discriminant = sqrt (a2 ** 2 - 4 * a1 * a3)
                dlambda1 = (-a2 - discriminant) / (2 * a1)
                dlambda2 = (-a2 + discriminant) / (2 * a1)
                
                dlambda_test = (/ &
                    dot_product (Uinc + dURflat + dlambda1 * dUFflat, Uinc), &
                    dot_product (Uinc + dURflat + dlambda2 * dUFflat, Uinc) &
                /)
                dlambda_test_res = maxloc (dlambda_test)
                if (dlambda_test_res (1) .eq. 1) then
                    dlambda = dlambda1
                else
                    dlambda = dlambda2
                end if
            end if
            
            lambda = lambda + dlambda
            Uinc = Uinc + dURflat + dlambda * dUFflat
            mesh%Displacements = mesh%Displacements + dUR + dlambda * dUF
            mesh%Positions = mesh%Coordinates + mesh%Displacements
            dth = thR + dlambda * thF
            call update_stress_strain (mesh, dth)
            Fint = assemble_internal_force (mesh)
            
            R = FC + lambda * Fext - Fint
            R2 (:, 1) = pack (Fext, .TRUE.)
            R2 (:, 2) = pack (R, .TRUE.)

            do j = 1, 6*mesh%NoNodes
                if (.NOT. dof (j)) R2(j, :) = 0.0D0
            end do
                        
            convtest = norm2 (R2 (:, 2)) / norm2 (R2 (:, 1))
            write (6, '(I4, X, "|", X, ES20.13, X, "|", X, F15.11)') i, convtest, lambda
            if (convtest < TOLERValue) exit
            
        end do
                
        if (i .eq. MAXITERValue + 1) then
            write (6, '(/, "Not converging")')
            stop
        end if
        
    end subroutine arcLength_iter
    
    
end module