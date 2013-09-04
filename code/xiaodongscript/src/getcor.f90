
program main
	implicit none

	integer, parameter :: char_len = 100
	integer :: i, i1, i2, j,  n_lines
	character(len=100) :: datafile, pairfile, outputfile, str
	double precision :: x, y, z
	double precision, allocatable :: coord_data(:,:)

	i = iargc()

	if(i .ne. 3) then
		print*
		print*, "ERROR! Number of arg must be 3."
		print*
		print*, "USAGE: "
		print*
		print*, "  ./getcor datafile pair_index_file pair_coord_file"
		print*
		print*, "     datafile: 3-d coordinates of a list of objects."
		print*, "     pair_index_file: index of objects pairs (start from 0). "	
		print*, "     pair_coord_file: coordinates of objects pairs (output, generated based on datafile and pari_indexfile)"
		stop
	endif

	call getarg(1,datafile)
	call getarg(2,pairfile)
	call getarg(3,outputfile)

!	datafile = "../testdata1"
!	pairfile = "../pairs_1.0.txt"
!	outputfile = "../pairs_1.0_xyz.txt"
	call read_in(datafile, 3, n_lines, coord_data)
!	write(*,*) "Reading in ", n_lines, " lines data from ", trim(adjustl(datafile))

!	do i = 1, n_lines
!		write(*,'(3(f14.7,1x))') coord_data(i,1:3)
!	enddo

	call count_line_number(pairfile, n_lines)
!	write(*,*) "Totally ", n_lines, " pairs"

	open (unit=1, file=pairfile)
	open (unit=2, file=outputfile)
	do i = 1, n_lines
		read(1,*) i1, i2
		i1 = i1+1; i2 = i2+1; !from c-convention to fortran
		if(i1<i2) cycle !ngl-beta double-counted pairs. reserve half 
		if(coord_data(i1,3) > coord_data(i2,3)) then !keep all vecotrs pointing 'up' (from low-z to high-z)
			j=i1; i1=i2;i2=j
		endif
		x = coord_data(i1,1)
		y = coord_data(i1,2)
		z = coord_data(i1,3)
		write(2,'(3(f14.7),1x)') x, y, z
		x = coord_data(i2,1)
		y = coord_data(i2,2)
		z = coord_data(i2,3)
		write(2,'(3(f14.7),1x)') x, y, z
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

