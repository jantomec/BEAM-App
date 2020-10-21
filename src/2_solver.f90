! Finite element solver
!
! author .............. Jan Tomec
! copyright ........... Copyright 2020, Project THREAD - University of Rijeka, Faculty of Civil Engineering
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
    
    implicit none
    
    private
    
    public :: newton_iter, newton_iter_det, arc_length_iter
    
    contains
    
    subroutine pack2 (A, DOF, B)
        
        implicit none
        
        double precision, dimension (:, :), intent (in) :: A
        logical, dimension (:), intent (in) :: DOF
        double precision, dimension (:, :), intent (out) :: B
        double precision, dimension (size (A (:, 1)), size (B (1, :))) :: AB
        integer :: j
        
        do j = 1, size (A (:, 1))
            AB (j, :) = pack (A (j, :), DOF)
        end do
        do j = 1, size (B (1, :))
            B (:, j) = pack (AB (:, j), DOF)
        end do
    
    end subroutine pack2
    
    subroutine logicwrite (A, DOF, B)
        
        implicit none
        
        double precision, dimension (:), intent (in) :: A
        logical, dimension (:), intent (in) :: DOF
        double precision, dimension (:), intent (out) :: B
        integer :: j, k
        
        k = 1
        do j = 1, size (B)  ! overwrite the result
            if (DOF (j) .eqv. .TRUE.) then
                B (j) = A (k)
                k = k + 1
            end if
        end do
        
    end subroutine logicwrite
    
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
    
    subroutine solve (A, b, c, dof, errck)
    
        implicit none
        
        double precision, dimension (:, :)    , intent (in)     :: A     ! (6*no nodes, 6*no nodes)
        double precision, dimension (:)       , intent (inout)  :: b     ! (6*no nodes)
        double precision, dimension (:)       , intent (in)     :: dof   ! (6*no nodes)
        integer                               , intent (out)    :: errck
        
        integer :: ndof, i, info
        integer, dimension (size (dof)) :: ipiv
        
        ndof = size (dof)
        
        do i = 1, ndof
            if (.NOT. dof (i)) then
                tangent (:, i) = 0.0D0
                tangent (i, :) = 0.0D0
                tangent (i, i) = 1.0D0
            end if
        end do
        
        call dgetrf (ndof, ndof, A, ndof, ipiv, info)  ! LU factorization
        if (info .ne. 0) then
            errck = 2
            exit
        end if
        
        call dgetrs ('N', ndof, 1, tangent, ndof, ipiv, b, ndof, info)  ! solve system
        if (info .ne. 0) then
            errck = 3
            exit
        end if
        
    end subroutine solve
    
    ! perform newton-raphson iteration
    !
    ! ele .......... elements array (no ele, no nodes on ele)
    ! X0 ........... initial coordinates (3, no all nodes)
    ! U ............ displacement (3, no all nodes)
    ! C ............ tangent elastic moduli (6, 6)
    ! DOF .......... degrees of freedom (6, no all nodes)
    ! Uload ........ prescribed displacement/rotation (6, no all nodes)
    ! Q ............ nodal load (6, no all nodes)
    ! pressure ............ distributed load (no ele, 6, no gauss)
    ! rot .......... rotations in gauss points (no ele, no gauss, 3, 3)
    ! om ........... curvature vector (no ele, 3, no gauss)
    ! stress ............ internal force vector (no ele, 3, no gauss)
    ! resout ....... residual vector (6, no all nodes)
    ! TOLER ........ tolerance for convergence
    ! MAXITER ...... maximum no of iterations
    ! TEST ......... convergence test ('RSD' - resdidual, 'DSP' - displacement)
    ! Niter ........ integer, number of iterations used to converge
    ! errck ........ integer, indication of error:
    !                    0 = no error
    !                    2 = LU Factorization error (dgetrf)
    !                    3 = Solver error (dgetrs)
    ! prints ....... boolean, indicating if intermediate info statements should be printed
    !
    ! modify U, rot, om, stress
    subroutine newton_iter (ele, X0, U, C, DOF, Uload, Q, pressure, rot, om, stress, resout, TOLER, MAXITER, TEST, Niter, prints)
    
        implicit none
        
        integer, dimension (:, :), intent (in) :: ele  ! (no ele, no nodes on ele)
        double precision, dimension (:, :), intent (in) :: X0  ! (3, no all nodes)
        double precision, dimension (:, :), intent (inout) :: U  ! (3, no all nodes)
        double precision, dimension (6, 6), intent (in) :: C
        logical, dimension (:, :), intent (in) :: DOF  ! (6, no all nodes)
        double precision, dimension (:, :), intent (in) :: Uload, Q  ! (6, no all nodes)
        double precision, dimension (:, :, :), intent (in) :: pressure  ! (no ele, 6, no nodes on ele)
        double precision, dimension (:, :, :, :), intent (inout) :: rot  ! (no ele, no gauss, 3, 3)
        double precision, dimension (:, :, :), intent (inout) :: om  ! (no ele, 3, no gauss)
        double precision, dimension (:, :, :), intent (inout) :: stress  ! (no ele, 6, no gauss)
        double precision, dimension (:, :), intent (out) :: resout  ! (6, no all nodes)
        double precision, intent (in) :: TOLER
        integer, intent (in) :: MAXITER
        character (len = 3), intent (in) :: TEST
        integer, intent (out) :: Niter
        logical, intent (in) :: prints
        
        integer :: nno, ndof, i, j, k, info
        double precision, dimension (6 * size (X0 (1, :)), 6 * size (X0 (1, :))) :: tangent
        double precision, dimension (3, size (X0 (1, :))) :: X, dU, dth
        double precision, dimension (6 * size (X0 (1, :))) :: Fint, Fext, R, resflat
        double precision, dimension (6, size (X0 (1, :))) :: res
        integer, dimension (6 * size (X0 (1, :))) :: ipiv
        logical, dimension (6 * size (X0 (1, :))) :: DOFflat
        double precision :: convtest
        
        external dgetrf, dgetrs
                
        nno = size (X0 (1, :))
        ndof = 6 * nno
        
        res = Uload
        DOFflat = pack (DOF, .TRUE.)
        
        if (prints) call begin_table ('N')
        
        do i = 0, MAXITER-1
            
            Niter = i + 1
            
            if (i > 0) then
                tangent = assemble_tangent (ele, X0, X, rot, C, stress)  ! tangent
                
                do j = 1, ndof
                    if (.NOT. DOFflat (j)) then
                        tangent (:, j) = 0.0D0
                        tangent (j, :) = 0.0D0
                        tangent (j, j) = 1.0D0
                    end if
                end do
                
                call solve ('N', ndof, 1, tangent, ndof, ipiv, R, ndof)
                                
                resflat = R
                res = reshape (resflat, (/ 6, nno /))
            end if
            
            dU = res (1:3, :)
            dth = res (4:6, :)
            
            U = U + dU
            X = X0 + U
                        
            call update_stress_strain (ele, X0, X, dth, C, rot, om, stress)
            
            Fint = assemble_internal_force (ele, X0, X, stress)
            Fext = assemble_external_force (ele, X0, Q, pressure)
            R = Fext - Fint
            
            do j = 1, ndof
                if (.NOT. DOFflat (j)) R(j) = 0.0D0
            end do
            
            if (TEST .eq. 'RSD') convtest = norm2 (R)
            if (TEST .eq. 'DSP') convtest = norm2 (resflat)
            
            if (prints) write (6, '(I4, X, "|", X, ES20.13)') i, convtest
            if (convtest < TOLER) exit
                        
        end do
        
        resout = reshape (R, (/ 6, nno /))
        
        if (i .eq. MAXITER) then
            if (prints) write (6, '(/, "Not converging")')
            stop
        end if
        
    end subroutine newton_iter
    
    ! perform newton-raphson iteration with determinant sign print
    !
    ! ele .......... elements array (no ele, no nodes on ele)
    ! X0 ........... initial coordinates (3, no all nodes)
    ! U ............ displacement (3, no all nodes)
    ! C ............ tangent elastic moduli (6, 6)
    ! DOF .......... degrees of freedom (6, no all nodes)
    ! Uload ........ prescribed displacement/rotation (6, no all nodes)
    ! Q ............ nodal load (6, no all nodes)
    ! pressure ............ distributed load (no ele, 6, no gauss)
    ! rot .......... rotations in gauss points (no ele, no gauss, 3, 3)
    ! om ........... curvature vector (no ele, 3, no gauss)
    ! stress ............ internal force vector (no ele, 3, no gauss)
    ! resout ....... residual vector (6, no all nodes)
    ! TOLER ........ tolerance for convergence
    ! MAXITER ...... maximum no of iterations
    ! TEST ......... convergence test ('RSD' - resdidual, 'DSP' - displacement)
    !
    ! modify U, rot, om, stress
    subroutine newton_iter_det (ele, X0, U, C, DOF, Uload, Q, pressure, rot, om, stress, resout, TOLER, MAXITER, TEST, Niter, errck)
    
        implicit none
        
        integer, dimension (:, :), intent (in) :: ele  ! (no ele, no nodes on ele)
        double precision, dimension (:, :), intent (in) :: X0  ! (3, no all nodes)
        double precision, dimension (:, :), intent (inout) :: U  ! (3, no all nodes)
        double precision, dimension (6, 6), intent (in) :: C
        logical, dimension (:, :), intent (in) :: DOF  ! (6, no all nodes)
        double precision, dimension (:, :), intent (in) :: Uload, Q  ! (6, no all nodes)
        double precision, dimension (:, :, :), intent (in) :: pressure  ! (no ele, 6, no nodes on ele)
        double precision, dimension (:, :, :, :), intent (inout) :: rot  ! (no ele, no gauss, 3, 3)
        double precision, dimension (:, :, :), intent (inout) :: om  ! (no ele, 3, no gauss)
        double precision, dimension (:, :, :), intent (inout) :: stress  ! (no ele, 6, no gauss)
        double precision, dimension (:, :), intent (out) :: resout  ! (6, no all nodes)
        double precision, intent (in) :: TOLER
        integer, intent (in) :: MAXITER
        character (len = 3), intent (in) :: TEST
        integer, intent (out) :: Niter
        integer, intent (out) :: errck
        
        integer :: nno, ndof, i, j, info
        logical, dimension (6 * size (X0 (1, :))) :: DOFsel
        double precision, dimension (6 * size (X0 (1, :)), 6 * size (X0 (1, :))) :: tangent
        double precision, allocatable :: K2 (:, :)
        double precision, dimension (6, size (X0 (1, :))) :: res
        double precision, dimension (6 * size (X0 (1, :))) :: res2
        double precision, allocatable :: R1 (:), R2 (:, :)
        double precision, dimension (3, size (X0 (1, :))) :: X, dU, dth
        double precision, dimension (6 * size (X0 (1, :))) :: Fint, Fext, R
        integer, allocatable :: ipiv (:)
        double precision :: convtest, tandet
        
        errck = 0.0D0
        
        nno = size (X0 (1, :))
        
        res = Uload
        
        DOFsel = pack (DOF, .TRUE.)
        ndof = count (DOFsel)
        allocate (K2 (ndof, ndof))
        allocate (R1 (ndof), R2 (ndof, 1), ipiv (ndof))
        
        call begin_table ('N')
        
        do i = 0, MAXITER-1
            
            Niter = i + 1
            
            if (i > 0) then
                res = 0.0D0
                tangent = assemble_tangent (ele, X0, X, rot, C, stress)  ! tangent
                
                call pack2 (tangent, DOFsel, K2)
                
                call dgetrf (ndof, ndof, K2, ndof, ipiv, info)  ! LU factorization
                if (info .ne. 0) then
                    errck = 2
                    exit
                end if  
                
                tandet = 1.0D0
                do j = 1, ndof  ! compute sign of determinant
                    tandet = K2 (j, j) / abs (K2 (j, j)) * tandet
                    if (j .ne. ipiv (j)) then
                        tandet = -1 * tandet
                    end if
                end do
                
                write (6, '("Sign of determinant:", X, F5.1)') tandet
                                
                call dgetrs ('N', ndof, 1, K2, ndof, ipiv, R2, ndof, info)  ! solve system
                if (info .ne. 0) then
                    errck = 3
                    exit
                end if  
                
                res2 = 0.0D0
                call logicwrite (R2 (:, 1), DOFsel, res2)
                res = reshape (res2, (/ 6, nno /))
                
            end if
            
            do concurrent (j = 1:3)
                dU (j, :) = res (j, :) 
                dth (j, :) = res (j + 3, :)
            end do
            
            U = U + dU
            X = X0 + U
            
            call update_stress_strain (ele, X0, X, dth, C, rot, om, stress)
            
            Fint = assemble_external_force (ele, X0, X, stress)
            Fext = assemble_internal_force (ele, X0, Q, pressure)
            R = Fint - Fext
            R1 = pack (R, DOFsel)
            
            if (TEST .eq. 'RSD') convtest = norm2 (R1)
            if (TEST .eq. 'DSP') convtest = norm2 (R2)
            
            write (6, '(I4, X, "|", X, ES20.13)') i, convtest
            if (convtest < TOLER) exit
            
            R2 (:, 1) = -R1  ! invert and verticalize residual
                        
        end do
        
        resout = reshape (R, (/ 6, nno /))
        
        if (i .eq. MAXITER) then
            write (6, '(/, "Not converging")')
            stop
        end if
        
    end subroutine newton_iter_det
    
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
    ! DOF .......... degrees of freedom (6, no all nodes)
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
    subroutine arc_length_iter (ele, X0, Uinc, U, C, DOF, Q, QC, pressure, rot, om, stress, resout, lambda, dS, TOLER, MAXITER, Niter, errck)
    
        implicit none
        
        integer, dimension (:, :), intent (in) :: ele  ! (no ele, no nodes on ele)
        double precision, dimension (:, :), intent (in) :: X0  ! (3, no all nodes)
        double precision, dimension (:), intent (inout) :: Uinc  ! (3 * no all nodes)
        double precision, dimension (:, :), intent (inout) :: U  ! (3, no all nodes)
        double precision, dimension (6, 6), intent (in) :: C
        logical, dimension (:, :), intent (in) :: DOF  ! (6, no all nodes)
        double precision, dimension (:, :), intent (in) :: Q, QC  ! (6, no all nodes)
        double precision, dimension (:, :, :), intent (in) :: pressure  ! (no ele, 6, no nodes on ele)
        double precision, dimension (:, :, :, :), intent (inout) :: rot  ! (no ele, no gauss, 3, 3)
        double precision, dimension (:, :, :), intent (inout) :: om  ! (no ele, 3, no gauss)
        double precision, dimension (:, :, :), intent (inout) :: stress  ! (no ele, 6, no gauss)
        double precision, dimension (:, :), intent (out) :: resout  ! (6, no all nodes)
        double precision, intent (inout) :: lambda
        double precision, intent (in) :: dS
        double precision, intent (in) :: TOLER
        integer, intent (in) :: MAXITER
        integer, intent (out) :: Niter
        integer, intent (out) :: errck
        
        integer :: nno, ndof, i, j, info
        logical, dimension (6 * size (X0 (1, :))) :: DOFsel
        double precision, dimension (6 * size (X0 (1, :)), 6 * size (X0 (1, :))) :: tangent
        double precision, allocatable :: K2 (:, :)
        double precision, dimension (6 * size (X0 (1, :))) :: res2F, res2R
        double precision, dimension (6, size (X0 (1, :))) :: resF, resR
        double precision, allocatable :: R2 (:, :)
        double precision, dimension (3, size (X0 (1, :))) :: X, dU, dth, U0, dUF, dUR, thF, thR
        double precision, dimension (3 * size (X0 (1, :))) :: dUFflat, dURflat
        double precision, dimension (6 * size (X0 (1, :))) :: Fint, Fext, FC, R
        integer, allocatable :: ipiv (:)
        double precision :: convtest, dlambda, dlambda1, dlambda2, discriminant, UincdotdUF, sig, a1, a2, a3
        double precision, dimension (2) :: dlambda_test
        integer, dimension (1) :: dlambda_test_res
        
        errck = 0
        
        nno = size (X0 (1, :))
                
        X = X0 + U
        DOFsel = pack (DOF, .TRUE.)
        ndof = count (DOFsel)
        allocate (K2 (ndof, ndof))
        allocate (R2 (ndof, 2), ipiv (ndof))
        
        Fint = assemble_external_force (ele, X0, X, stress)
        Fext = assemble_internal_force (ele, X0, Q, pressure)
        FC = assemble_internal_force (ele, X0, QC, pressure)
        R = Fint - lambda * Fext - FC
        R2 (:, 1) = pack (Fext, DOFsel)
        R2 (:, 2) = pack (-R, DOFsel)
        
        call begin_table ('A')
        write (6, '(I4, X, "|", X, ES20.13, X, "|", X, F14.10)') 0, norm2 (R), lambda
        
        do i = 1, MAXITER
            
            Niter = i
        
            res2F = 0.0D0
            res2R = 0.0D0
            
            tangent = assemble_tangent (ele, X0, X, rot, C, stress)  ! tangent
            call pack2 (tangent, DOFsel, K2)
            
            call dgetrf (ndof, ndof, K2, ndof, ipiv, info)  ! LU factorization
            if (info .ne. 0) then
                errck = 2
                exit
            end if          
            
            call dgetrs ('N', ndof, 2, K2, ndof, ipiv, R2, ndof, info)  ! solve system
            if (info .ne. 0) then
                errck = 3
                exit
            end if
            
            call logicwrite (R2 (:, 1), DOFsel, res2F)
            call logicwrite (R2 (:, 2), DOFsel, res2R)
            resF = reshape (res2F, (/ 6, nno /))
            resR = reshape (res2R, (/ 6, nno /))
            
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
            U = U + dUR + dlambda * dUF
            X = X0 + U
            dth = thR + dlambda * thF
            call update_stress_strain (ele, X0, X, dth, C, rot, om, stress)
            Fint = assemble_external_force (ele, X0, X, stress)
            R = Fint - lambda * Fext - FC
            R2 (:, 1) = pack (Fext, DOFsel)
            R2 (:, 2) = pack (-R, DOFsel)
            convtest = norm2 (R2 (:, 2)) / norm2 (R2 (:, 1))
            write (6, '(I4, X, "|", X, ES20.13, X, "|", X, F15.11)') i, convtest, lambda
            if (convtest < TOLER) exit
            
        end do
        
        resout = reshape (R, (/ 6, nno /))
        
        if (i .eq. MAXITER + 1) then
            write (6, '(/, "Not converging")')
            stop
        end if
        
    end subroutine arc_length_iter
    
    
end module