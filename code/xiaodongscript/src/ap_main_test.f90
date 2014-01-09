
program ap_main_test

!use mpi
use ap_tools
use ap_chisq
use ap_grids
use ap_structure_count

	implicit none

! Common variables
	integer :: i,j,k, ix,iy,iz, i1,i2,i3, RSD,AP
	real(dl) :: x,y,z
	character(len=char_len) :: outputfile,str1,str2
	type(chisq_settings) :: cs
	integer, parameter :: numbins = 3
	real(dl) :: time1, time2

! Variables for Test Polynomial fit
	integer, parameter :: ndat = 100, norder = 10
	real(dl) :: Xarray(ndat), Yarray(ndat), Yfit(ndat), A(norder+1)
! Variables for Test binned_quan
	integer, parameter :: num_quan = 100000, num_bin = 1
	real(dl) :: quan_val_list(num_quan), quan_r_list(num_quan), rmin, rmax, &
		quan_val_av_list(num_bin), quan_val_er_list(num_bin), r_av_list(num_bin), &
		quan_val_var_list(num_bin), quan_medval_list(num_bin)
		
! Variables for Test strucount
	real(dl) :: pixelsize
	real(dl) :: chisqcrA,chisqmpA,nowom
	real(dl) :: binnedcrratio(numbins),binnedcrratioer(numbins),binnedmpratio(numbins),binnedmpratioer(numbins)
	integer, parameter :: numom = 3
	real(dl), parameter :: ommin = 0.06, ommax = 0.46

!Test strucount
	if(.true.) then
		call cpu_time(time1)
		call random_seed()
		! chisq settings
		cs%smnum = 30
		cs%num_in_x = 100
		cs%check_boundary = .true.
		use_num_density = .true.
		cs%cb_adjust_ratio = 1.0_dl
		cs%remov_dist_ratio = 0.0_dl
		cs%nbins_rhoav = numbins
		gb_use_fixmd = .true.
		gb_fixmd = 15.0_dl
		gb_do_seg_cut = .false.

		! halo data
		str1 = "NewCBPHR3LC00_zge0_rge500" !'NewHR3LC00_xyzge0_rge500'!'MDR1_shelldata_sep15'
		halo_data_file = '../../../data/input/'//trim(adjustl(str1))//'.dat'
		local_pecu_v = .true.
		write(str2,*) cs%nbins_rhoav
		outputfile = 'Test/'//trim(adjustl(str1))//'_nbins'//trim(adjustl(str2))//'.dat'

		! print info
		cs%print_info = .false.

		! initialize
		call cosmo_funs_init()
		call read_in_halo_data()
		call init_halo_info()
		print *
		! test struc count
		open(unit=604582,file=outputfile)
		do i = 1, numom
		do RSD = 0, 1
			nowom = ommin + (ommax-ommin) / dble(numom-1) * dble(i-1)
			write(teststr,'(f5.2)') nowom
			print *, '##################################'
			if(RSD.eq.0) then
				write(str1,*) "NO RSD; "
				teststr = 'Test/'//trim(adjustl(teststr))//'_NORSD_'
			else
				write(str1,*) 'With RSD; '
				teststr = 'Test/'//trim(adjustl(teststr))//'_WithRSD_'				
			endif
			write(str2,'(A,f5.2)') 'Omegam = ', real(nowom)
			print *, trim(adjustl(str1))//' '//trim(adjustl(str2))
			print *, '##################################'
			call init_AP_cosmo(nowom,gb_w,gb_h,1,cs%print_info)

			pixelsize = 15.0_dl;
			print *, '  SIZE of PIXEL = ', real(pixelsize)
			print *
!			call strucount(RSD,1,cs,pixelsize,chisqcrA,chisqmpA,binnedcrratio,binnedcrratioer,binnedmpratio,binnedmpratioer)
			call strucount(RSD,1,cs,pixelsize)
			write(604582,'(i2,e14.7,<2+3+numbins*4>(e14.7,1x))') pixelsize,RSD,gb_omegam,gb_w,chisqcrA,chisqmpA,&
				binnedcrratio,binnedcrratioer,binnedmpratio,binnedmpratioer
			cs%print_info = .false.
			print *
		enddo
		enddo
		close(604582)
		call cpu_time(time2)
		print *, 'Time used: ', time2
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
			Yarray(i) = sin(Xarray(i))
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
