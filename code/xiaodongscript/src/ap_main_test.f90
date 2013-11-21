
program ap_main_test

!use mpi
use ap_tools
use ap_chisq
use ap_grids
use ap_structure_count

	implicit none

! Common variables
	integer :: i,j,k, ix,iy,iz, i1,i2,i3
	real(dl) :: x,y,z
	character(len=char_len) :: outputfile
	type(chisq_settings) :: cs

! Variables for Test Polynomial fit
	integer, parameter :: ndat = 100, norder = 10
	real(dl) :: Xarray(ndat), Yarray(ndat), Yfit(ndat), A(norder+1)
! Variables for Test binned_quan
	integer, parameter :: num_quan = 100000, num_bin = 1
	real(dl) :: quan_val_list(num_quan), quan_r_list(num_quan), rmin, rmax, &
		quan_val_av_list(num_bin), quan_val_er_list(num_bin), r_av_list(num_bin), &
		quan_val_var_list(num_bin), quan_medval_list(num_bin)
		
! Variables for Test strucount
		
!Test strucount
	if(.true.) then
		! chisq settings
		cs%smnum = 30
		cs%num_in_x = 100
		cs%check_boundary = .true.
		cs%print_info = .true.
		cs%use_num_density = .true.
		cs%cb_adjust_ratio = 1.0_dl
		cs%remov_dist_ratio = 0.0_dl
		cs%nbins_rhoav = 4
		use_fixmd = .true.
		gb_fixmd = 40.0_dl
		! halo data file
		halo_data_file = '../../../data/input/NewHR3LC00_xyzge0_rge500.dat'
		local_pecu_v = .true.
		! initialize
		call cosmo_funs_init()
		call read_in_halo_data()
		call init_halo_info()
		! test struc count
		call strucount(RSD=0, AP=0, cs=cs, pixelsize=40.0_dl)
		
		stop
	endif
	
! Test cell_index	
	if(.false.) then		
		do ix=1,gb_n_cellx
		do iy=1,gb_n_celly
		do iz=1,gb_n_cellz
			call cell_pos(ix,iy,iz,x,y,z)
			call cell_index(i1,i2,i3,x,y,z)
			call random_number(x)
			if(x .le. 0.0001) then
				print *, ix,i1,iy,i2,iz,i3
			endif
			if(i1.ne.ix .or. i2.ne.iy .or. i3.ne.iz) then
				print *, 'ERROR!!! : ', i1,ix,i2,iy,i3,iz
			endif
		enddo
		enddo
		enddo
	endif
	
! Test Polynomial fit		
	if(.false.) then
		do i = 1, ndat
			Xarray(i) = 0.1*(i-1)
			Yarray(i) = dsin(Xarray(i))
			print *, Yarray(i), sin(Xarray(i)), Yarray(i) - sin(Xarray(i))
		enddo
	
		call poly_fit(Xarray,Yarray,A,ndat,norder)
		do i = 1, ndat
			Yfit(i) = poly(Xarray(i),A,norder)
		enddo

	! 	output Xarray, Y, Yfit to files
		outputfile = 'Test/X.dat'
		call output_1d(outputfile, Xarray, ndat)
		outputfile = 'Test/Y.dat'
		call output_1d(outputfile, Yarray, ndat)
		outputfile = 'Test/Yfit.dat'
		call output_1d(outputfile, Yfit, ndat)
	endif
		
! Test binned_quan 
	if(.false.) then
		call random_number(quan_val_list)
		quan_val_list = quan_val_list**2.0
		quan_r_list = quan_val_list
		call find_min_max(quan_r_list, num_quan, rmin, rmax)
	!	do i = 1, num_quan
	!		quan_r_list(i) = rando
	!	enddo
		outputfile = 'Test/RandomQuans.txt'
		call output_1d(outputfile, quan_val_list, num_quan)
		call binned_quan(quan_val_list, quan_r_list, rmin, rmax, num_quan, num_bin, &
	  		quan_val_av_list, quan_val_er_list, r_av_list, quan_val_var_list, quan_medval_list)
	
		do i = 1, num_bin
			print *, 'Infor for the ', i, 'th bin:'
			print *, '  <quan> = ', real(quan_val_av_list(i)), ' +/- ', real(quan_val_er_list(i))
			print *, '  <r>    = ', real(r_av_list(i))
			print *, ' medium  = ', real(quan_medval_list(i))
			print *
		enddo
	endif
end program ap_main_test
