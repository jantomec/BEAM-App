program example

	use solver
	
	implicit none
	
	character (len=100)							:: LineIn = ''
	character (len=10)							:: string = ''
	integer										:: lenLineIn, startValue, endValue
	double precision							:: length, nel_d
	integer, dimension (1, 2)					:: ele
	double precision, dimension (3, 2)			:: X0, U
	double precision, dimension (6, 6)			:: C
	logical, dimension (6, 2)					:: DOF
	double precision, dimension (6, 2)			:: dU, Q, res
	double precision, dimension (1, 6, 1)		:: p, f
	double precision, dimension (1, 3, 1)		:: om
	double precision, dimension (1, 1, 3, 3)	:: rot
	integer, parameter							:: MAXITER = 20
	double precision, parameter					:: TOLER = 1D-5
	integer										:: i, Niter, errck

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
	
	! Write the html results page to the browser,
	! with the sum of the two numbers.
	write (*, '(1X,"<html><body>")')
	
	ele (1, :) = (/ 1, 2 /)
	X0 (:, 1) = (/ 0.0D1, 0.0D1, 0.0D1 /)
	X0 (:, 2) = (/ 1.0D1, 2.0D1, 3.0D1 /)
	U = 0.0D1
	C = 0.0D1
	do i = 1, 6
		C (i, i) = 1.0D1
	end do
	DOF (:, 1) = .FALSE.
	DOF (:, 2) = .TRUE.
	dU = 0.0D1
	Q = 0.0D1
	Q (2, 2) = 0.001
	p = 0.0D1
	f = 0.0D1
	om = 0.0D1
	rot = 0.0D1
	do i = 1, 3
		rot (1, 1, i, i) = 1
	end do
	
	call newton_iter (ele, X0, U, C, DOF, dU, Q, p, rot, om, f, res, TOLER, MAXITER, 'RSD', Niter, errck, .FALSE.)
	
	write (*, '(1X, "<p>Success! The length is given below:</p>", E12.4)') length
	write (*, '(1X, "<p>Success! The no of elements is given below:</p>", E12.4)') nel
	write (*, '(1X,"</html></body>")')
	
end program example
