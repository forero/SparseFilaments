
PROGRAM MAIN

	implicit none

	integer, parameter :: char_len = 200
	integer :: i, n_lines
	integer :: n_first_ignore, n_last_ignore, x_column !  the first/last several lines to ignore. The last  
	integer(kind=8) :: id
	character(len=char_len) :: inputfile, outputfile
	character(len=1000) str
	double precision :: x, y, z, temp(100)

	i = iargc()
	
	if(i .ne. 2) then
		print*
		print*, "ERROR! Number of arg must be 2."
		print*
		print*, "USAGE: "
		print*
		print*, "  ./post_proc inputfile outputfile"
		print*
		stop
	endif

	call getarg(1,inputfile)
	call getarg(2,outputfile)

	n_first_ignore = 37
	n_last_ignore = 1
	x_column = 5
	
	call count_line_number(inputfile, n_lines)
	
!	write(*,*) "n_lines = ", n_lines
	
	open(unit=1,file=inputfile)
	open(unit=2,file=outputfile)
	do i = 1, n_first_ignore
		read(1,*) str
	enddo
	
	do i = n_first_ignore+1, n_lines-n_last_ignore
		read(1,*) temp(1:x_column-1), x,y,z, str
		write(2,'(3(f14.7,1x))') x,y,z
	enddo
	
	close(1)
	close(2)
	
contains

	SUBROUTINE count_line_number (file_name, line_number)
		INTEGER :: line_number
		CHARACTER(LEN=char_len) :: file_name, inline
    
		OPEN(UNIT=1456,NAME=file_name,ERR=2)
    
		line_number = 0
		DO While(1 .EQ. 1)
			READ(1456,*,ERR=3,END=100) inline
			line_number = line_number + 1
		ENDDO

2		WRITE(*,*) "Error occurs when opening the file ", file_name
		CLOSE(1456)
		STOP

3		WRITE(*,*) "Error occurs when couting the lines of the file ", file_name
		CLOSE(1456)
		STOP

100		CLOSE(1456)
	END SUBROUTINE count_line_number
	
	SUBROUTINE read_in (file_name, n_column, n_lines, total_data)
		CHARACTER(LEN=char_len), INTENT(IN) :: file_name
		INTEGER, INTENT(IN)  :: n_column
		INTEGER, INTENT(OUT) :: n_lines
		DOUBLE PRECISION, ALLOCATABLE :: total_data(:,:)
		INTEGER :: i
		CALL count_line_number(file_name, n_lines)
		ALLOCATE(total_data(n_lines, n_column))
		OPEN(UNIT=4321,FILE=file_name)
		DO i = 1, n_lines
			READ(4321, *) total_data(i,1:n_column)
		ENDDO
		CLOSE(4321)
	END SUBROUTINE read_in 
	
END PROGRAM MAIN
