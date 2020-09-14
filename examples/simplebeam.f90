program example

	use solver
	use mesher
	use mesh_object
	
	implicit none
	
	character (len=100)							:: LineIn = ''
	character (len=10)							:: string = ''
	integer										:: lenLineIn, startValue, endValue
	double precision							:: length, nel_d
	
	integer										:: Nele, Nno, order, Ngauss
	type (ElementMesh) 							:: mesh
	double precision, allocatable 				:: X (:,:), U (:,:)
	double precision, dimension (6, 6)			:: C
	logical, allocatable						:: DOF (:,:)
	double precision, allocatable				:: dU (:,:), Q (:,:)
	double precision, allocatable				:: p (:,:,:), f (:,:,:)
	double precision, allocatable				:: om (:,:,:)
	double precision, allocatable				:: rot (:,:,:,:)
	double precision, allocatable				:: R (:,:)
	double precision, dimension (6), parameter	:: material = (/ 1.0D1, 1.0D1, 1.0D1, 1.0D1, 2.0D1, 1.0D1 /)
	integer, parameter							:: MAXITER = 20
	double precision, parameter					:: TOLER = 1D-5
	integer										:: i, j, Niter, info
	
	! Determine the length of the data string
    ! that is to be read in the next section.
	call Get_Environment_Variable('CONTENT_LENGTH', string)
	read (string, *) lenLineIn

	! Read the data from the html form into the
	! variable LineIn, a single character at a time.
	do i = 1, lenLineIn
		read (*, ADVANCE='NO', FMT='(A)') LineIn(i:i)
	end do
	
	! Locate and read the value of 'length' from LineIn (name in html form)
	startValue = index (LineIn, 'length=') + 7
	endValue = startValue + index(LineIn(startValue:), '&') - 2
	read (LineIn(startValue:endValue), *) length
	
	! Locate and read the value of 'nel' from LineIn
	startValue = index(LineIn,'nel=') + 4
	read (LineIn(startValue:), *) nel_d
	
	! Send a header to the browser, identifying
	! the type of information that will be sent.
	write (*, '("Content-type: text/html",//)')
	
	! Begin with html
	write (*, '(1X,"<html><body>")')
	
	! =================================================
	! MESH
	Nele = int (nel_d)
	order = 1
	Ngauss = order
	Nno = Nele * order + 1
	mesh%Nno = Nno
	mesh%Nele = Nele
	mesh%order = 1
	call lmsh (length, mesh)
	do i = 1, Nele
		do j = 1, Ngauss
			rot (i, j, 1, :)= (/ 1.0D1, 0.0D1, 0.0D1 /)
			rot (i, j, 2, :)= (/ 0.0D1, 1.0D1, 0.0D1 /)
			rot (i, j, 3, :)= (/ 0.0D1, 0.0D1, 1.0D1 /)
		end do
	end do
	
	allocate (X (3, Nno), U (3, Nno), DOF (6, Nno), dU (6, Nno), Q (6, Nno))
	allocate (p (Nele, 6, Ngauss), f (Nele, 6, Ngauss), om (Nele, 3, Ngauss))
	allocate (rot (Nele, Ngauss, 3, 3), R (6, Nno))
	
	! =================================================
	! ELASTIC MODULI MATRIX
	C = 0.0D1  ! add material
	do i = 1, 6
		C (i, i) = material (i)
	end do
	
	
	! =================================================
	! DATA INITIALIZATION
	U = 0.0D1
	p = 0.0D1
	f = 0.0D1
	om = 0.0D1
	
	! =================================================
	! BOUNDARY CONDITIONS
	DOF = .TRUE.
	DOF (:, 1) = .FALSE.
	dU = 0.0D1
	Q = 0.0D1
	
	! =================================================
	! FORCE CONTROL ROUTINE
	do j = 0, 10
		if (j > 0) then
			Q (5, Nno) = Q (5, Nno) + 25 / 10
			call newton_iter (mesh%ele, mesh%X0, U, C, DOF, dU, Q, p, rot, om, f, R, TOLER, MAXITER, 'RSD', Niter, info, .FALSE.)
		end if
		
		X (1, :) = mesh%X0 (1, :) + U (1, :)
		X (2, :) = mesh%X0 (2, :) + U (2, :)
		X (3, :) = mesh%X0 (3, :) + U (3, :)
				
	end do
		
	! =================================================
	! HTML OUTPUT
	
	write (*, '(1X, "<p>Success!</p>")')
	!call htmlmatrix (X, 3, Nno)
	write (*, '(1X,"</html></body>")')
	
end program example

subroutine htmlmatrix (A, ndim, n)
		
	implicit none
	
	double precision, dimension (ndim,n), intent(in) 	:: A
	integer 											:: i, ndim, n
	character (len = 80)								:: arrayfmt
	
	write (arrayfmt, '(A, i5, A)') '("<p>", ', ndim, 'f5.1, "</p>")'
	
	do i = 1, size(A (1, :))
		write (*, arrayfmt) A (:, i)
	end do
		
end subroutine htmlmatrix
