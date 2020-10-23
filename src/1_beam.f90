! Beam element
!
! author ....... Jan Tomec
! copyright .... Copyright 2020, Project THREAD - University of Rijeka, Faculty of Civil Engineering
! credits ...... Jan Tomec, Gordan Jelenić
! license ...... GPL
! version ...... 1.0.0
! maintainer ... Jan Tomec
! email ........ jan.tomec@gradri.uniri.hr
! status ....... Development
! date ......... 06/05/2020
!
! ------------------------------------------------------------------------------

module beam
    
    use shape_functions
    use legendre_gauss
    use vector_algebra
    use mesh_objects
    
    implicit none
    
    private
    
    public :: assemble_tangent, assemble_external_force, assemble_internal_force, update_stress_strain
    
    contains
    
    ! compute length of an element
    function element_length (coordinates, element)
        
        implicit none
        
        double precision, dimension (:, :)                                  :: coordinates
        type (LineElement)                                                  :: element
        double precision                                                    :: element_length
        
        double precision, dimension (element%NoGauss)                       :: pts, wgts
        double precision, dimension (element%NoNodes, element%NoGauss)      :: dN
        double precision, dimension (3, element%NoNodes)                    :: X0
        double precision, dimension (3, element%NoGauss)                    :: dX0
        double precision, dimension (3)                                     :: intg
        integer                                                             :: i
        
        call legauss (element%NoGauss, pts, wgts)
        dN = shdfun (element%NoNodes, pts)
        X0 = coordinates (:, element%Nodes)
        
        dX0 = matmul (X0, dN)
        
        do i = 1, 3
            intg (i) = dot_product (dX0 (i, :), wgts)
        end do
        element_length = norm2 (intg)
    
    end function element_length
    
    ! compute single Pi matrix
    function Pi_matrix (rot)
    
        implicit none
        
        double precision, dimension (3, 3), intent (in) :: rot
        double precision, dimension (6, 6) :: Pi_matrix
        
        Pi_matrix = 0.0D0
        Pi_matrix (1:3, 1:3) = rot
        Pi_matrix (4:6, 4:6) = rot
        
    end function Pi_matrix
    
    ! compute single Xi matrix
    function Xi_matrix (dX, N, dN)
    
        implicit none
        
        double precision, dimension (3), intent (in) :: dX
        double precision, intent (in) :: N, dN
        double precision, dimension (6, 6) :: Xi_matrix
        double precision, dimension (3, 3) :: S
        integer :: i
        
        S = skew (dX)
        Xi_matrix = 0.0D0
        do i = 1, 6
            Xi_matrix (i, i) = dN
        end do
        Xi_matrix (4, 2) = - N * S (1, 2)
        Xi_matrix (4, 3) = - N * S (1, 3)
        Xi_matrix (5, 1) = - N * S (2, 1)
        Xi_matrix (5, 3) = - N * S (2, 3)
        Xi_matrix (6, 1) = - N * S (3, 1)
        Xi_matrix (6, 2) = - N * S (3, 2)
        
    end function Xi_matrix
    
    ! compute elemental internal forces in nodes
    function internal_forces (coordinates, positions, element) result (F)
        
        implicit none
        
        double precision, dimension (:, :)                              :: coordinates, positions
        type (LineElement)                                              :: element
        double precision, dimension (6, element%NoNodes)                :: F
        
        double precision, dimension (3, element%NoNodes)                :: X
        double precision                                                :: L
        double precision, dimension (element%NoGauss)                   :: pts, wgts
        double precision, dimension (element%NoNodes, element%NoGauss)  :: N, dN
        double precision, dimension (3, element%NoGauss)                :: dX
        double precision, dimension (3, 3)                              :: S
        double precision, dimension (6, 6)                              :: Xi_i
        integer                                                         :: g, i
        
        X = positions (:, element%Nodes)
        
        L = element_length (coordinates, element)
        call legauss (element%NoGauss, pts, wgts)
        N = shfun (element%NoNodes, pts)
        dN = 2.0D0 / L * shdfun (element%NoNodes, pts)
        
        dX = matmul (X, dN)

        F = 0.0D0
        
        do g = 1, element%NoGauss
            do i = 1, element%NoNodes
                Xi_i = Xi_matrix (dX (:, g), N (i, g), dN (i, g))
                F (:, i) = F (:, i) + matmul (Xi_i, element%Stress (:, g)) * wgts (g)
            end do
        end do
        
        F = L / 2.0D0 * F
    
    end function internal_forces
    
    ! compute elemental external forces in nodes
    function external_forces (coordinates, element) result (F)
        
        implicit none
        
        double precision, dimension (:, :)                                      :: coordinates
        type (LineElement)                                                      :: element
        double precision, dimension (6, element%NoNodes)                        :: F
        
        double precision                                                        :: L
        double precision, dimension (element%NoGauss)                           :: pts, wgts
        double precision, dimension (element%NoNodes, element%NoGauss)          :: N
        double precision, dimension (6, 6)                                      :: H
        integer                                                                 :: g, i
                
        L = element_length (coordinates, element)
        call legauss (element%NoGauss, pts, wgts)
        N = shfun (element%NoNodes, pts)
        
        F = 0.0D0
        
        do g = 1, element%NoGauss
            H = 0.0D0
            do i = 1, element%NoNodes
                H (1, 1) = N (i, g)
                H (2, 2) = N (i, g)
                H (3, 3) = N (i, g)
                H (4, 4) = N (i, g)
                H (5, 5) = N (i, g)
                H (6, 6) = N (i, g)
                F (:, i) = F (:, i) + matmul (H, element%Pressure (:, g)) * wgts (g)
            end do
        end do
        
        F = L / 2.0D0 * F
                
    end function external_forces
    
    ! compute elemental curvature, internal stresses and new rotation
    subroutine curvature (coordinates, positions, element, theta)
        
        implicit none
        
        double precision, dimension (:, :)                                      :: coordinates
        double precision, dimension (:, :)                                      :: positions
        type (LineElement)                                                      :: element
        double precision, dimension (:, :)                                      :: theta
        
        double precision, dimension (3, element%NoNodes)                        :: X
        double precision                                                        :: L
        double precision, dimension (element%NoGauss)                           :: pts, wgts, tn
        double precision, dimension (element%NoNodes, element%NoGauss)          :: N, dN
        double precision, dimension (3, element%NoGauss)                        :: dX, t, dt
        integer                                                                 :: g
        
        double precision, dimension (3, 3)          :: R, rotinv
        double precision, dimension (3)             :: om, b1, b2, b3, a1, a2, a3, Gamma, kappa, fn, fm
        double precision, dimension (3), parameter  :: E3 = (/ 1.0D0, 0.0D0, 0.0D0 /)
        
        X = positions (:, element%Nodes)
        
        L = element_length (coordinates, element)
        call legauss (element%NoGauss, pts, wgts)
        N = shfun (element%NoNodes, pts)
        dN = 2.0D0 / L * shdfun (element%NoNodes, pts)
        
        dX = matmul (X, dN)
        t = matmul (theta, N)
        dt = matmul (theta, dN)
        tn = norm2 (t, dim=1)
                
        do g = 1, element%NoGauss
            om = element%strain (4:6, g)
            R = exponentialMap (t (:, g))
            element%RotationMatrix (g, :, :) = matmul (R, element%RotationMatrix (g, :, :))
            if (tn (g) == 0) then
                om = om + dt (:, g)
            else
                b1 = (1.0D0 - sin (tn (g)) / tn (g)) &
                   * dot_product (t (:, g), dt (:, g)) / tn (g) ** 2 * t (:, g)
                b2 = sin (tn (g)) / tn (g) * dt (:, g)
                b3 = (1.0D0 - cos (tn (g))) / tn (g) ** 2 &
                   * cross_product (t (:, g), dt (:, g))
                a1 = cos (tn (g)) * om
                a2 = (1.0D0 - cos (tn (g))) / tn (g) ** 2 &
                   * dot_product (t (:, g), om) * t (:, g)
                a3 = sin (tn (g)) / tn (g) * cross_product (t (:, g), om)
                om = b1 + b2 + b3 + a1 + a2 + a3
            end if
            
            rotinv = transpose (element%RotationMatrix (g, :, :))
            Gamma = matmul (rotinv, dX (:, g))
            Gamma = Gamma - E3
            kappa = matmul (rotinv, om)
            element%strain (1:3, g) = matmul (element%RotationMatrix (g, :, :), Gamma)
            element%strain (4:6, g) = om
                        
            fn = matmul (element%C (1:3, 1:3), Gamma)
            fm = matmul (element%C (4:6, 4:6), kappa)
            element%Stress (1:3, g) = matmul (element%RotationMatrix (g, :, :), fn)
            element%Stress (4:6, g) = matmul (element%RotationMatrix (g, :, :), fm)
            
        end do
        
    end subroutine curvature
    
    function material_stiffness (coordinates, positions, element) result (K)
        
        implicit none
        
        double precision, dimension (:, :)                                      :: coordinates
        double precision, dimension (:, :)                                      :: positions
        type (LineElement)                                                      :: element
        double precision, dimension (6 * element%NoNodes, 6 * element%NoNodes)  :: K
        
        double precision, dimension (3, element%NoNodes)                        :: X
        double precision                                                        :: L
        double precision, dimension (element%NoGauss)                           :: pts, wgts
        double precision, dimension (element%NoNodes, element%NoGauss)          :: N, dN
        double precision, dimension (3, element%NoGauss)                        :: dX
        integer                                                                 :: g, i, j, ii, ij
        double precision, dimension (6, 6) :: Pi_g, Xi_i, Xi_j, c_g
                
        X = positions (:, element%Nodes)
        
        L = element_length (coordinates, element)
        call legauss (element%NoGauss, pts, wgts)
        N = shfun (element%NoNodes, pts)
        dN = 2.0D0 / L * shdfun (element%NoNodes, pts)
        
        dX = matmul (X, dN)
        
        K = 0.0D0
        
        do g = 1, element%NoGauss
            Pi_g = Pi_matrix (element%RotationMatrix (g, :, :))
            c_g = matmul (matmul (Pi_g, element%C), transpose (Pi_g))
            do i = 1, element%NoNodes
                ii = 6 * (i-1) + 1
                Xi_i = Xi_matrix (dX (:, g), N (i, g), dN (i, g))
                do j = 1, element%NoNodes
                    Xi_j = Xi_matrix (dX (:, g), N (j, g), dN (j, g))
                    ij = 6 * (j-1) + 1
                    K (ii:ii + 5, ij:ij + 5) = K (ii:ii + 5, ij:ij + 5) + &
                        wgts (g) * matmul (matmul (Xi_i, c_g), transpose (Xi_j))
                end do
            end do
        end do
        
        K = L / 2.0D0 * K
    
    end function material_stiffness
    
    function geometrical_stiffness (coordinates, positions, element) result (K)
        
        implicit none
        
        double precision, dimension (:, :)                                      :: coordinates
        double precision, dimension (:, :)                                      :: positions
        type (LineElement)                                                      :: element
        double precision, dimension (6 * element%NoNodes, 6 * element%NoNodes)  :: K
        
        double precision, dimension (3, element%NoNodes)                        :: X
        double precision                                                        :: L
        double precision, dimension (element%NoGauss)                           :: pts, wgts
        double precision, dimension (element%NoNodes, element%NoGauss)          :: N, dN
        double precision, dimension (3, element%NoGauss)                        :: dX
        integer                                                                 :: g, i, j, ii, ij
                
        X = positions (:, element%Nodes)
        
        L = element_length (coordinates, element)
        call legauss (element%NoGauss, pts, wgts)
        N = shfun (element%NoNodes, pts)
        dN = 2.0D0 / L * shdfun (element%NoNodes, pts)
        
        dX = matmul (X, dN)
        
        K = 0.0D0
        
        do g = 1, element%NoGauss
            do i = 1, element%NoNodes
                ii = 6 * (i-1) + 1
                do j = 1, element%NoNodes
                    ij = 6 * (j-1) + 1
                    K (ii:ii + 2, ij + 3:ij + 5) = &
                        K (ii:ii + 2, ij + 3:ij + 5) - &
                        wgts (g) * skew (element%Stress (1:3, g)) * dN (i, g) * N (j, g)
                    K (ii + 3:ii + 5, ij:ij + 2) = &
                        K (ii + 3:ii + 5, ij:ij + 2) + &
                        wgts (g) * skew (element%Stress (1:3, g)) * N (i, g) * dN (j, g)
                    K (ii + 3:ii + 5, ij + 3:ij + 5) = &
                        K (ii + 3:ii + 5, ij + 3:ij + 5) + &
                        wgts (g) * ( &
                            -skew (element%Stress (4:6, g)) * dN (i, g) * N (j, g) &
                            + matmul (skew (dX (:, g)), skew (element%Stress (1:3, g))) &
                            * N (i, g) * N (j, g) &
                        )
                end do
            end do
        end do
        
        K = L / 2.0D0 * K
    
    end function geometrical_stiffness
    
    function mass_stiffness (h, beta, gamma, coordinates, element) result (K)
        
        implicit none
        
        type (LineElement)                 :: element
        double precision                   :: h, beta, gamma
        double precision, dimension (:, :) :: coordinates
        
        double precision, dimension (6 * element%NoNodes, 6 * element%NoNodes)  :: K
        
        double precision                                                        :: L
        double precision, dimension (element%NoGauss)                           :: pts, wgts
        double precision, dimension (element%NoNodes, element%NoGauss)          :: N, dN
        double precision, dimension (6, 6)                                      :: m
        double precision, dimension (3)                                         :: m11
        double precision, dimension (3, 3)                                      :: m1, m2, m3, m4, m3T, m3tp
        double precision                                                        :: t
        integer                                                                 :: g, i, j, ii, ij
        
        L = element_length (coordinates, element)
        call legauss (element%NoGauss, pts, wgts)
        N = shfun (element%NoNodes, pts)
        dN = 2.0D0 / L * shdfun (element%NoNodes, pts)
                
        K = 0.0D0
        
        do g = 1, element%NoGauss
            do i = 1, element%NoNodes
                ii = 6 * (i-1) + 1
                do j = 1, element%NoNodes
                    ij = 6 * (j-1) + 1
                    m = 0.0D0
                    m (1:3, 1:3) = wgts (g) * ( 1 /(h**2 * beta) * element%rho * element%A * N (i, g) * N (j, g) )
                    m11 = matmul ( &
                        element%RotationMatrix (g, :, :), &
                        matmul (element%InertiaMatrix, element%AngularAcceleration (:, g)) + &
                            cross_product (element%AngularVelocity (:, g), &
                                matmul (element%InertiaMatrix, element%AngularVelocity (:, g)) &
                            ) &
                    )
                    m1 = skew (m11)
                    m2 = matmul ( &
                        element%RotationMatrix (g, :, :), &
                        element%InertiaMatrix - &
                            h * gamma * skew (matmul (element%InertiaMatrix, element%AngularVelocity (:, g))) + &
                            h * gamma * matmul (skew (element%AngularVelocity (:, g)), element%InertiaMatrix) &
                    )
                    t = norm2 (element%Rotation (:, g))
                    m3tp = tensor_product(element%Rotation (:, g), element%Rotation (:, g)) / (t**2)
                    m3T = m3tp + t/2 / tan (t/2) * (identityMatrix (3) - m3tp) - 0.5D0 * skew (element%Rotation (:, g))
                    m3 = matmul (element%RotationMatrixLastConverged (g, :, :), m3T)
                    m4 = matmul(-m1 + 1 /(h**2 * beta) * m2, m3)
                    m (4:6, 4:6) = wgts (g) * ( m4 * N (i, g) * N (j, g) )
                    K (ii:ii + 5, ij:ij + 5) = K (ii:ii + 5, ij:ij + 5) + m
                end do
            end do
        end do
        
        K = L / 2.0D0 * K
    
    end function mass_stiffness
    
    ! compute stiffness of an element
    function assemble_tangent (mesh) result (Kg)
    
        implicit none
        
        type (ElementMesh)                                                      :: mesh
        double precision, dimension (6 * mesh%NoNodes, 6 * mesh%NoNodes)        :: Kg
        
        type (LineElement)  :: element
        integer             :: e, i, j, ei, ej
        double precision, dimension (:,:), allocatable :: Ke
        
        Kg = 0.0D0
        
        do e = 1, mesh%NoElements
            element = mesh%Elements (e)
            allocate (Ke (6*element%NoNodes, 6*element%NoNodes))
            Ke = material_stiffness (mesh%Coordinates, mesh%Positions, element) &
                + geometrical_stiffness (mesh%Coordinates, mesh%Positions, element)
            do i = 1, element%NoNodes
                ei = element%Nodes (i)
                do j = 1, element%NoNodes
                    ej = element%Nodes (j)
                    Kg (6 * (ei - 1) + 1:6 * ei, 6 * (ej - 1) + 1:6 * ej) = &
                        Kg (6 * (ei - 1) + 1:6 * ei, 6 * (ej - 1) + 1:6 * ej) &
                        + Ke (6 * (i - 1) + 1:6 * i, 6 * (j - 1) + 1:6 * j)
                end do
            end do
            deallocate (Ke)
        end do
        
    end function assemble_tangent
    
    function assemble_tangent_dynamic (mesh, h, beta, gamma) result (Kg)
    
        implicit none
        
        type (ElementMesh) :: mesh
        double precision   :: h, beta, gamma
        
        double precision, dimension (6 * mesh%NoNodes, 6 * mesh%NoNodes) :: Kg
        
        type (LineElement)  :: element
        integer             :: e, i, j, ei, ej
        double precision, dimension (:,:), allocatable :: Ke
        
        Kg = 0.0D0
        
        do e = 1, mesh%NoElements
            element = mesh%Elements (e)
            allocate (Ke (6*element%NoNodes, 6*element%NoNodes))
            Ke = material_stiffness (mesh%Coordinates, mesh%Positions, element) + &
                     geometrical_stiffness (mesh%Coordinates, mesh%Positions, element) + &
                     mass_stiffness (h, beta, gamma, mesh%Coordinates, element)
            do i = 1, element%NoNodes
                ei = element%Nodes (i)
                do j = 1, element%NoNodes
                    ej = element%Nodes (j)
                    Kg (6 * (ei - 1) + 1:6 * ei, 6 * (ej - 1) + 1:6 * ej) = &
                        Kg (6 * (ei - 1) + 1:6 * ei, 6 * (ej - 1) + 1:6 * ej) &
                        + Ke (6 * (i - 1) + 1:6 * i, 6 * (j - 1) + 1:6 * j)
                end do
            end do
            deallocate (Ke)
        end do
        
    end function assemble_tangent_dynamic
    
    ! compute global internal force vector
    function assemble_internal_force (mesh) result (Fint)
    
        implicit none
        
        type (ElementMesh)                                                      :: mesh
        double precision, dimension (6, mesh%NoNodes)                           :: Fint
        
        type (LineElement)  :: element
        integer             :: j, i, ei
        double precision, dimension (:,:), allocatable :: Fj
                
        Fint = 0.0D0
        
        do j = 1, mesh%NoElements
            element = mesh%Elements (j)
            allocate (Fj (6, element%NoNodes))
            Fj = internal_forces (mesh%Coordinates, mesh%Positions, element)
            do i = 1, element%NoNodes
                ei = element%Nodes (i)
                Fint (:, ei) = Fint (:, ei) + Fj (:, i)
            end do
            deallocate (Fj)
        end do
    
    end function assemble_internal_force
    
    ! compute global external force vector
    function assemble_external_force (mesh, Q) result (Fext)
    
        implicit none
        
        type (ElementMesh)                                                      :: mesh
        double precision, dimension (6, mesh%NoNodes)                           :: Q
        double precision, dimension (6, mesh%NoNodes)                           :: Fext
        
        type (LineElement)  :: element
        integer             :: j, i, ei
        double precision, dimension (:,:), allocatable :: Fj
                
        Fext = 0.0D0
        
        do j = 1, mesh%NoElements
            element = mesh%Elements (j)
            allocate (Fj (6, element%NoNodes))
            Fj = external_forces (mesh%Coordinates, element)
            do i = 1, element%NoNodes
                ei = element%Nodes (i)
                Fext (:, ei) = Fext (:, ei) + Fj (:, i)
            end do
            deallocate (Fj)
        end do
        
        Fext = Fext + Q
    
    end function assemble_external_force
    
    ! compute global external force vector
    subroutine update_stress_strain (mesh, theta)
    
        implicit none
        
        type (ElementMesh)                                                      :: mesh
        double precision, dimension (:,:)                                       :: theta
        
        type (LineElement)  :: element
        integer             :: j
        
        do j = 1, mesh%NoElements
            element = mesh%Elements (j)
            call curvature (mesh%Coordinates, mesh%Positions, element, theta(:, element%Nodes))
            mesh%Elements (j) = element
        end do
        
    end subroutine update_stress_strain
    
    
end module

