!####################################
! scan chisqs
!####################################


program ap_main

!use mpi
use ap_chisq
use ap_structure_count

	implicit none

	type(chisq_settings) :: cs
	logical :: auto_num_in_x
	character(len=char_len) :: inputfile, outputname, &
		rhofile, rhoRSDfile, dltfile, dltRSDfile, ndltfile, ndltRSDfile, str, str2
	integer :: i, j, n, num_omw, tmpnum
	real(dl), allocatable :: om_w_list(:,:)
	! chis
	real(dl), allocatable :: rhochi(:), rhoRSDchi(:), rhodfchi(:), rhoRSDdfchi(:), &
		dltchi(:), dltRSDchi(:), dltdfchi(:), dltRSDdfchi(:), &
		ndltchi(:), ndltRSDchi(:), ndltdfchi(:), ndltRSDdfchi(:)
!	real(dl) :: rhochi(gb_numwei), rhoRSDchi(gb_numwei)
	real(dl) :: ommin, ommax, wmin, wmax, dom, dw, om, w, time1, time2, timebegin, timeend
	! tmp variables used for testing
	real(dl) :: drmean, drvar, ratio
	real(dl), allocatable :: tmp(:), drho_mu_data(:), ddlt_mu_data(:), dnormed_dlt_mu_data(:)
	real(dl) :: chisq1, chisq2, chisq3
	integer :: i_LC, datalab
	character(len=char_len) :: mustr, nbliststr
	
	i = iargc()
	
	if(i .ne. 2) then
		print*
		print*, "ERROR! Number of arg must be 2."
		print*
		print*, "USAGE: "
		print*
		print*, "  ./post_proc inputfile outputfile"
		print*
!		stop
	endif

	call getarg(1,inputfile)
	call getarg(2,outputname)
	
	open(unit=1,file=inputfile)
	
	read(1,*) cs%smnum, cs%num_in_x, auto_num_in_x, cs%check_boundary, cs%print_info, cs%cb_adjust_ratio, &
		cs%remov_dist_ratio, cs%use_num_density, cs%nbins_rhoav, cs%use_intpl_rho, cs%numdrop, gb_num_changenuminx, gb_amp_changenuminx, num_omw
	
	write(*,*) 'Settings of the programm:'
	write(*,*) '    smnum             = ', cs%smnum
	if(.not.auto_num_in_x) &
		write(*,*) '    num_in_x          = ', cs%num_in_x
	write(*,*) '    auto_num_in_x     = ', auto_num_in_x
	write(*,*) '    check_boundary    = ', cs%check_boundary
	write(*,*) '    print_info        = ', cs%print_info
	write(*,*) '    cb_adjust_ratio   = ', real(cs%cb_adjust_ratio)
	write(*,*) '    remov_dist_ratio  = ', real(cs%remov_dist_ratio)
	write(*,*) '    use_num_density   = ', cs%use_num_density
	write(*,*) '    nbins_rhoav       = ', cs%nbins_rhoav
	write(*,*) '    use_intpl_rho     = ', cs%use_intpl_rho
	write(*,'(A,i4)')    '     # of cellchange   : ', gb_num_changenuminx
	write(*,'(A,e14.7)') '     amp of cellchange : ', gb_amp_changenuminx	
	write(*,*) '    num of drop       = ', cs%numdrop
	write(*,*) '    num of om/w       = ', num_omw

	allocate(cs%dropval(cs%numdrop),  cs%lowdropvalratio(cs%numdrop),  cs%highdropvalratio(cs%numdrop), &
		 cs%dropdval(cs%numdrop), cs%lowdropdvalratio(cs%numdrop), cs%highdropdvalratio(cs%numdrop))	
	
	print *
	do i = 1, cs%numdrop
		read(1,*) cs%dropval(i), cs%lowdropvalratio(i),  cs%highdropvalratio(i), &
			cs%dropdval(i), cs%lowdropdvalratio(i), cs%highdropdvalratio(i)
		cs%lowdropvalratio(i) = cs%lowdropvalratio(i) / 100.0_dl
		cs%highdropvalratio(i) = cs%highdropvalratio(i) / 100.0_dl
		cs%lowdropdvalratio(i) = cs%lowdropdvalratio(i) / 100.0_dl
		cs%highdropdvalratio(i) = cs%highdropdvalratio(i) / 100.0_dl
		write(*,*) 'Droppint pixels settings ', i,':'
		write(*,*) '    Density:            ', cs%dropval(i), real(cs%lowdropvalratio(i)),  real(cs%highdropvalratio(i))
		write(*,*) '    Density Gradient:   ', cs%dropdval(i), real(cs%lowdropdvalratio(i)), real(cs%highdropdvalratio(i))
	enddo

	! read in the om_w_list
	allocate(om_w_list(2,num_omw))
	print *
	do i = 1, num_omw
		read(1,*) om_w_list(1:2,i)
		if(gbtp) print*, '  read in omegam/w = ', real(om_w_list(1:2,i))
	enddo
	close(1)

	i_LC = 0
	write(LClabelChar,*) i_LC

!BEGIN TESTING	
	if(.true.) then
		datalab = 3
		if(datalab .eq. 1) then
			halo_data_file = '../../data/input/MD_fullhalos_z0_100w.dat'; 
			local_pecu_v = .false.
			mustr = 'chisqlike_LookAtRSD/MDR1muinfo_'
		elseif(datalab .eq. 2) then
			halo_data_file = '../../data/input/HR3LCQua'//trim(adjustl(LClabelChar))//'_57w_600to1787.dat'; 
			local_pecu_v = .false.
			mustr = 'chisqlike_LookAtRSD/HR3LCQua0muinfo_'
		elseif(datalab .eq. 3) then
			halo_data_file = '../../data/input/NewHR3LC00.dat'
			local_pecu_v = .true.
			mustr = 'chisqlike_LookAtRSD/NewHR3LC00muinfo_'
		else
			print *, 'Wrong datalab (must be 1 or 2 or 3)!'
			stop
		endif
	endif	
!END TESTING	

	dotsbe = .false.
	
	if(.not.dotsbe) then
!CCC 	Testing of several MDR1 datasets
!		halo_data_file = '../../data/input/MDR1_shelldata_sep15.dat'
!		halo_data_file = '../../data/input/MDR1_shelldata_sep10.dat'
!		halo_data_file = '../../data/input/MDR1_shelldata_sep8.dat'
!		halo_data_file = '../../data/input/NewHR3LC00_xyzge0.dat'
!		halo_data_file = '../../data/input/NewHR3LC00.dat'			
		halo_data_file = '../../data/input/NewHR3LC00_xyzge0_rge500.dat'
!		halo_data_file = '../../data/input/NewCBPHR3LCx.dat'
		local_pecu_v = .true.
	endif
!	halo_data_file = '../../data/input/NewHR3LC00_xyzge0_rge500.dat'
!	halo_data_file = '../../data/input/HR3LC0_1eighth.dat'
!	halo_data_file = '../../data/input/HR3LCQua'//trim(adjustl(LClabelChar))//'_57w_600to1787.dat'; local_pecu_v = .true.
!	halo_data_file = '../../data/input/HR3LC0_NMC_600to1787.dat'
!	halo_data_file = '../../data/input/HR3LCQua'//trim(adjustl(LClabelChar))//'_57w_1200to1787.dat'

	if(use_fixmd) then
		print *
		print *, 'Applying fixed max_dist = ', gb_fixmd
		print *
	endif
	
	if(use_seg) then
		print *
		print *, 'Applying seg_cut_dist = ', seg_cut_dist
		print *
	endif
	
	
	
!	call cosmo_funs_init()
!	call read_in_halo_data()
!	call init_halo_info()

	allocate(rhochi(cs%numdrop),rhoRSDchi(cs%numdrop),rhodfchi(cs%numdrop),rhoRSDdfchi(cs%numdrop), &
		dltchi(cs%numdrop),dltRSDchi(cs%numdrop),dltdfchi(cs%numdrop),dltRSDdfchi(cs%numdrop), &
		ndltchi(cs%numdrop),ndltRSDchi(cs%numdrop),ndltdfchi(cs%numdrop),ndltRSDdfchi(cs%numdrop))

	str = outputname
	rhofile = trim(adjustl(str))//'_rho.txt'
	rhoRSDfile = trim(adjustl(str))//'_rhoRSD.txt'
	dltfile = trim(adjustl(str))//'_dlt.txt'
	dltRSDfile = trim(adjustl(str))//'_dltRSD.txt'
	ndltfile = trim(adjustl(str))//'_ndlt.txt'
	ndltRSDfile = trim(adjustl(str))//'_ndltRSD.txt'
	
	print *
	print *, 'Calculating chisqs...'	
	print *, 'Results saved in the following file:'
	print *, '     ', trim(adjustl(rhofile))
	print *, '     ', trim(adjustl(rhoRSDfile))
	print *, '     ', trim(adjustl(dltfile))
	print *, '     ', trim(adjustl(dltRSDfile))
	print *, '     ', trim(adjustl(ndltfile))
	print *, '     ', trim(adjustl(ndltRSDfile))
	print *

	open(unit=1,file=rhofile)
	open(unit=2,file=rhoRSDfile)
	open(unit=3,file=dltfile)
	open(unit=4,file=dltRSDfile)
	open(unit=5,file=ndltfile)
	open(unit=6,file=ndltRSDfile)
	
	call cpu_time(timebegin)
	time1 = timebegin

	print *, 'Scanning with both RSD, no RSD...'
	do i = 1, num_omw
		om = om_w_list(1,i)
		w  = om_w_list(2,i)

		if(gbtp) then
			write(*,'(A,i8,1x,f8.4,1x,f8.4,1x,f7.3,A)') &
				'  Step, om, w, ratio = ', i, real(om), real(w), real(i/(num_omw+0.0)), '...'
		endif
		
		cs%has_RSD = .false.
		
!		call gf_mldprho_mlchi2s(om, w, h_dft, cs, rhochi, rhodfchi, dltchi, dltdfchi, &
!			ndltchi, ndltdfchi, cs%numdrop, calc_comvr = .true.)
		call gf_mldprho_chi2s(om, w, h_dft, cs, rhochi, rhodfchi, cs%numdrop, calc_comvr = .true.)

		write(1,'(<2+cs%numdrop+cs%numdrop>(e14.7,1x))') om, w, rhochi(1:cs%numdrop), rhodfchi(1:cs%numdrop)
		write(3,'(<2+cs%numdrop+cs%numdrop>(e14.7,1x))') om, w, dltchi(1:cs%numdrop), dltdfchi(1:cs%numdrop)
		write(5,'(<2+cs%numdrop+cs%numdrop>(e14.7,1x))') om, w,ndltchi(1:cs%numdrop),ndltdfchi(1:cs%numdrop)
		
		print *, 'Chisqs without RSD...'
		write(*,'(<2+cs%numdrop+cs%numdrop>(f7.2,1x))') om, w, rhochi(1:cs%numdrop), rhodfchi(1:cs%numdrop)
		write(*,'(<2+cs%numdrop+cs%numdrop>(f7.2,1x))') om, w, dltchi(1:cs%numdrop), dltdfchi(1:cs%numdrop)
		write(*,'(<2+cs%numdrop+cs%numdrop>(f7.2,1x))') om, w,ndltchi(1:cs%numdrop),ndltdfchi(1:cs%numdrop)
		
		cs%has_RSD = .true.
		
!		call gf_mldprho_mlchi2s(om, w, h_dft, cs, rhoRSDchi, rhoRSDdfchi, dltRSDchi, dltRSDdfchi, &
!			ndltRSDchi, ndltRSDdfchi, cs%numdrop, calc_comvr = .true.)
		call gf_mldprho_chi2s(om, w, h_dft, cs, rhoRSDchi, rhoRSDdfchi, cs%numdrop, calc_comvr = .true.)
		
		if(gbtp) then
			write(*,'(A,2f10.4,A)') '  Mean d_sep / r_sm_sphere = ', &
				(gbtotvol/num_halo)**(1.0/3.0), est_sm_sphe_r(gbtotvol, num_halo, cs%smnum), '  ...'
		endif
			
		write(2,'(<2+cs%numdrop+cs%numdrop>(e14.7,1x))') om, w, rhoRSDchi(1:cs%numdrop), rhoRSDdfchi(1:cs%numdrop)
		write(4,'(<2+cs%numdrop+cs%numdrop>(e14.7,1x))') om, w, dltRSDchi(1:cs%numdrop), dltRSDdfchi(1:cs%numdrop)
		write(6,'(<2+cs%numdrop+cs%numdrop>(e14.7,1x))') om, w,ndltRSDchi(1:cs%numdrop),ndltRSDdfchi(1:cs%numdrop)

		print *, 'Chisqs with RSD...'
		write(*,'(<2+cs%numdrop+cs%numdrop>(f7.2,1x))') om, w, rhoRSDchi(1:cs%numdrop), rhoRSDdfchi(1:cs%numdrop)
		write(*,'(<2+cs%numdrop+cs%numdrop>(f7.2,1x))') om, w, dltRSDchi(1:cs%numdrop), dltRSDdfchi(1:cs%numdrop)
		write(*,'(<2+cs%numdrop+cs%numdrop>(f7.2,1x))') om, w,ndltRSDchi(1:cs%numdrop),ndltRSDdfchi(1:cs%numdrop)
		
		call cpu_time(time2)
		if(gbtp) then
			write(*,'(A,f10.4)') '  Time used in this step: ', time2-time1
			time1 = time2
			print *
		endif
	enddo
	
	call cpu_time(timeend)
	
	print *, 'Total time:     ',  timeend - timebegin
	print *, 'Total time / n: ', (timeend - timebegin) / (num_omw+0.0)

	deallocate(om_w_list, rhochi,rhoRSDchi,rhodfchi,rhoRSDdfchi,dltchi,&
		dltRSDchi,dltdfchi,dltRSDdfchi,ndltchi,ndltRSDchi,ndltdfchi,ndltRSDdfchi)

	close(1); close(2); 
end program ap_main
