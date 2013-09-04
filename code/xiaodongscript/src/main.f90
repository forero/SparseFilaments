!####################################
!This module does statistical an
!####################################


program main
use ap_tools
use ap_cosmo_funs
use smooth
	implicit none

	integer :: i,j,k, ix, iy, iz, n, m
	integer :: di,imin,imax,mini,maxi
	integer, allocatable :: selected_list(:)
	character(len=char_len) :: xyz_mass_file_name
	real(dl) :: x, y, z, time1, time2, rho, drhodx, drhody, drhodz
	real(dl) :: test_array(1000), test_2d_array(3,100)
	real(dl), allocatable :: ref_1d_array(:)
	integer, allocatable :: index_array(:)

	!Test the get_z function
	if (1 .eq. 2) then
		gb_omegam = 0.27
		gb_h = 0.73
		gb_w = -1.0
		call de_tools_init()
		call de_calc_comovr()
!		do i = 1, 100
!			print *, de_zdata(i), de_comovr_data(i), comov_r(de_zi(i))
!		enddo
!		do i = 1, 100
!			z = i * 0.01d0
!			print *, z, (comov_r(z) - de_get_comovr(z))/comov_r(z)
!		enddo
		z = 0.53
		call cpu_time(time1)
		do i = 1, 1
			x = de_get_comovr(z)
		enddo
		print *, x
		call cpu_time(time2)
		print *, time2 - time1
		stop
	endif
	
	!Test the bubble sort
	if (1 .eq. 2) then
		call random_seed()
		call random_number(test_2d_array)
		n = size(test_2d_array,2)
		m = size(test_2d_array,1)
		allocate(ref_1d_array(n))
		ref_1d_array = 0
		print *, 'array not sorted:'
		do i = 1, n
			write(*,'(i4,1x,<m>(e14.7,1x))'), i, test_2d_array(1:m,i)	
			do j = 1, m
				ref_1d_array(i) = ref_1d_array(i) &
					+ test_2d_array(j,i)**2.0
			enddo
		enddo
		call bubble_2d(test_2d_array,ref_1d_array)
		print *, 'array sorted:'
		do i = 1, n
			write(*,'(i4,1x,<m+1>(e14.7,1x))'), i, test_2d_array(1:m,i),	 ref_1d_array(i)
			x = 0
			do j = 1, m
				x = x + test_2d_array(j,i)**2.0
			enddo
			write(*,*) x
		enddo
		stop
	endif	
	
	xyz_mass_file_name = '../../../data/output/MD_halos_20w/MD_halos_20w__xyz_mass__om_1p0_RSD.dat'
	
!	use_sd_unit_len = .false.
	
!	print*, 'sd_unit_len = ', sd_unit_len   
	num_in_x = 50
	call do_cell_initialize(xyz_mass_file_name)
	
!	! interactive test of cell_mat			
!	if(1.eq.2) then
!101		read(*,*) ix,iy,iz
!		if(ix .le. n_cellx .and. iy .le. n_celly .and. iz.le.n_cellz) then	
!			print*, 'ix,iy,iz = ', ix,iy,iz
!			print*, 'list of halos:', cell_mat(ix,iy,iz)%list(1:cell_mat(ix,iy,iz)%halo_num)
!			write(*,'(A,e14.7)'), ' total mass:', mass_mat(ix,iy,iz)
!			call cell_pos(ix,iy,iz,x,y,z)
!			print*, 'cell_position = ', x,y,z
!			print*, 'cell rho = ', cell_rho(ix,iy,iz)
!		else
!			print*, 'ix,iy,iz overflow!'
!		endif
!		goto 101
!	endif
	
	! interactive test of ilist
!	do while(1 .eq. 2) 
!		read(*,*) i, di, imin, imax
!		call ilist(i,di,imin,imax,mini,maxi)
!		print*, mini,maxi
!	enddo
		
	!select halos in cells nearby
	ix = 4; iy = 9; iz = 6;
	call cell_pos(ix,iy,iz,x,y,z)
	x=1100;y=1200;z=1050;
	call cpu_time(time1)
	do i = 1, 10
		call nb_list(x,y,z,200,rho,drhodx,drhody,drhodz)
	enddo
	call cpu_time(time2)
	print *, 'cpu_time = ', time2-time1
   	print *, 'position:', x,y,z
   	print *, 'rho, drho: ', rho, drhodx, drhody, drhodz
   	stop
!   	call cpu_time(time1)
!  	do i = 1, 100000
!		call nb_select(x,y,z,300,selected_list=selected_list)
!	enddo
!	call cpu_time(time2)
!	print*, 'cpu_time = ', time2-time1
!	print*, selected_list
!	print*, 'size of selected = ', size(selected_list,1)

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

	!call getarg(1,datafile)
	!call getarg(2,pairfile)
	!call getarg(3,outputfile)

end program main

