! Files and folders management routines
!
! author ....... Jan Tomec
! copyright .... Copyright 2020, Project THREAD - University of Rijeka, Faculty of Civil Engineering
! credits ...... Jan Tomec, Gordan Jelenić
! license ...... GPL
! version ...... 1.0.0
! maintainer ... Jan Tomec
! email ........ jan.tomec@gradri.uniri.hr
! status ....... Development
! date ......... 22/10/2020
!
! ------------------------------------------------------------------------------

module files
    
    use mesh_objects
	
	implicit none
	
	private
	
	public :: removeFolder, createFolder, writeResults
	
	contains
	
    subroutine removeFolder (folderName)
	
		implicit none
        
        character (len = *), intent (in) :: folderName
        character (len = 255)            :: wdir_long
        character (len = :), allocatable :: wdir
        character (len = 255)            :: command
		integer                          :: st
        
		! Get current working directoy
        st = getcwd (wdir_long)
        wdir = trim (wdir_long)
        
        ! Remove old folder with results if it exists
        write (command, '("if exist", X, A, "\", A, X, "rmdir", X, "/S", X, "/Q", X, A, "\" A)') wdir, folderName, wdir, folderName
        call execute_command_line (command, .TRUE.)
		
	end subroutine removeFolder
    
    subroutine createFolder (folderName)
	
		implicit none
        
        character (len = *), intent (in) :: folderName
        character (len = 255)            :: wdir_long
        character (len = :), allocatable :: wdir
        character (len = 255)            :: command
		integer                          :: st
        
		! Get current working directoy
        st = getcwd (wdir_long)
        wdir = trim (wdir_long)
        
        ! Create folder for results
        write (command, '("mkdir", X, A, "\" A)') wdir, folderName
        call execute_command_line (command, .TRUE.)
		
	end subroutine createFolder
    
    subroutine writeResults(mesh, R, j, folderName)
        
        implicit none
        
        type (ElementMesh),                            intent (in) :: mesh
        double precision, dimension (6, mesh%NoNodes), intent (in) :: R
        integer,                                       intent (in) :: j
        character (len = *),                           intent (in) :: folderName
        
        character (len = *), parameter   :: fname_format = '("step", I0.3, ".dat")'
        character (len = 255)            :: fname_long
        character (len = :), allocatable :: fname
        character (len = :), allocatable :: relpath
        integer                          :: i
                
        ! Write new results file name to variable
        write (fname_long, fname_format) j
        fname = trim (fname_long)
        
        ! Open new file with designated name
        relpath = folderName//'/'//fname
        open (unit = 12, file = relpath, status = 'unknown', action = 'write')
        
        ! Write header line
        write (12, '("node", 19X, "X1", 19X, "X2", 19X, "X3", 19X, "U1", 19X, "U2", 19X, "U3", 19X, "R1", 19X, "R2", 19X, &
            "R3", 19X, "M1", 19X, "M2", 19X, "M3")')
        
        ! Write for each node
        do i = 1, mesh%NoNodes
            write (12, '(I0.4, X, ES20.11E3, X, ES20.11E3, X, ES20.11E3, X, ES20.11E3, X, ES20.11E3, X, ES20.11E3, X, ES20.11E3, &
                   X, ES20.11E3, X, ES20.11E3, X, ES20.11E3, X, ES20.11E3, X, ES20.11E3)') i, &
                   mesh%Positions (1, i), mesh%Positions (2, i), mesh%Positions (3, i), &
                   mesh%Displacements (1, i), mesh%Displacements (2, i), mesh%Displacements (3, i), &
                   R (1, i), R (2, i), R (3, i), R (4, i), R (5, i), R (6, i)
        end do
        
        ! Close file
        close(12)
        
    end subroutine
	
end module