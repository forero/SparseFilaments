
program main
	implicit none

	integer, parameter :: char_len = 100
	integer :: iarg, i, m, n, n_lines
	character(len=100) :: datafile, outputfile, str
	double precision :: x, y, z, xmin, xmax, ymin, ymax, zmin, zmax
	double precision, allocatable :: coord_data(:,:)

	iarg = iargc()
	if((iarg .ne. 4) .and. (iarg .ne. 8)) then
		print*
		print*, "ERROR! Number of arg must be 4 or 8."
		print*
		print*, "USAGE: "
		print*
		print*, "  ./select_obj datafile outputfile xmin xmax ymin ymax zmin zmax"
		print*, "  ./select_obj datafile outputfile min max"
		print*
		print*, "     datafile: 3-d coordinates of a list of objects."
		print*, "     xmin, max, ymin, ...: ranges of x, y, z"
		stop
	endif

	call getarg(1,datafile)
	call getarg(2,outputfile)
	call getarg(3,str)
	read(str,*,err=100) xmin
	call getarg(4,str)
	read(str,*,err=100)  xmax
	if(iarg .eq. 4) then
		ymin = xmin; ymax = xmax;
		zmin = xmin; zmax = xmax;
		else
		call getarg(5,str)
		read(str,*,err=100) ymin
		call getarg(6,str)
		read(str,*,err=100) ymax
		call getarg(7,str)
		read(str,*,err=100) zmin
		call getarg(8,str)
		read(str,*,err=100) zmax
	endif
	goto 101
	
100	print*, 'ERROR OCCURS DURING READ!'
	stop
	
101	continue
!	write(*,*) "xmin, xmax, ymin, ymax, zmin, zmax = ", xmin, xmax, ymin, ymax, zmin, zmax 
!	datafile = "../testdata1"
!	pairfile = "../pairs_1.0.txt"
!	outputfile = "../pairs_1.0_xyz.txt"
	call read_in(datafile, 3, n_lines, coord_data)
!	write(*,*) "Reading in ", n_lines, " lines data from ", trim(adjustl(datafile))

!	do i = 1, n_lines
!		write(*,'(3(f14.7,1x))') coord_data(i,1:3)
!	enddo

!	write(*,*) "Totally ", n_lines, " pairs"

	open (unit=1, file=datafile)
	open (unit=2, file=outputfile)
	do i = 1, n_lines
		read(1,*) x, y, z
!		write(*,*) x, y, z
		if ((x.le.xmax) .and. (x.ge.xmin) &
			.and. (y.le.ymax) .and. (y.ge.ymin) &
			.and. (z.le.zmax) .and. (z.ge.zmin)) then
			write(2,'(3(f14.7),1x)') x, y, z
		endif
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
end program main
