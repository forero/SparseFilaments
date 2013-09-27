!####################################
!This module does statistical an
!####################################

!		integer :: smnum, num_in_x
!		logical :: check_boundary = .true.
!		logical :: print_info = .false.
!		real(dl) :: cb_adjust_ratio = 1.0_dl, remov_dist_ratio = 1.0_dl
!		logical :: use_num_density = .true.
!		integer :: nbins_rhoav = 10
!		logical :: use_intpl_rho = .true.
!		logical :: has_RSD = .true.
!		integer, allocatable :: nbins_list(:)

! List of abondoned settings:

! gb_chisq_method
! cs%nbins_rhoav
! cs%nbins_list(:)

program ap_main

use mpi
use ap_chisq

	implicit none

	type(chisq_settings) :: cs
	character(len=char_len) :: inputfile, outputname, &
		rhofile, rhoRSDfile, deltafile, deltaRSDfile, ndeltafile, ndeltaRSDfile, str, str2
	integer :: i, j, n, num_omw
	real(dl), allocatable :: om_w_list(:,:), rho_chisqlist(:), delta_chisqlist(:), ndelta_chisqlist(:),  &
			rhoRSD_chisqlist(:),deltaRSD_chisqlist(:),ndeltaRSD_chisqlist(:)
	real(dl) :: ommin, ommax, wmin, wmax, dom, dw, om, w, time1, time2, timebegin, timeend
	logical :: auto_num_in_x, scan1, scan2	
	
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

	print *
	if (.not. gb_chisq_initied) then
		call cosmo_funs_init()
		call read_in_halo_data()
		call init_halo_info()
		gb_chisq_initied = .true.
	endif

	!settings of smnum
	if(auto_num_in_x) then
		cs%num_in_x = floor(est_num_in_x(cs%smnum)+0.5)
		print *, 'Automatically determine num_in_x = ', cs%num_in_x, '...'
	endif	
	
!	stop
	
	str = outputname
	rhofile = trim(adjustl(str))//'_rho.txt'
	rhoRSDfile = trim(adjustl(str))//'_rhoRSD.txt'
!	deltafile = trim(adjustl(str))//'_delta.txt'
!	deltaRSDfile = trim(adjustl(str))//'_deltaRSD.txt'
!	ndeltafile = trim(adjustl(str))//'_ndelta.txt'
!	ndeltaRSDfile = trim(adjustl(str))//'_ndeltaRSD.txt'
	
	print *
	print *, 'Calculating chisqs...'	
	print *, 'Results saved in the following file:'
	print *, '     ', trim(adjustl(rhofile))
	print *, '     ', trim(adjustl(rhoRSDfile))
!	print *, '     ', trim(adjustl(deltafile))
!	print *, '     ', trim(adjustl(deltaRSDfile))
!	print *, '     ', trim(adjustl(ndeltafile))
!	print *, '     ', trim(adjustl(ndeltaRSDfile))
	print *
	
	open(unit=1,file=rhofile)
	open(unit=2,file=rhoRSDfile)
!	open(unit=3,file=deltafile)
!	open(unit=4,file=deltaRSDfile)
!	open(unit=5,file=ndeltafile)
!	open(unit=6,file=ndeltaRSDfile)
	
	allocate(rho_chisqlist(cs%numdrop),delta_chisqlist(cs%numdrop),ndelta_chisqlist(cs%numdrop), &
		rhoRSD_chisqlist(cs%numdrop),deltaRSD_chisqlist(cs%numdrop),ndeltaRSD_chisqlist(cs%numdrop))
		
	call cpu_time(timebegin)
	time1 = timebegin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	scan1 = .false. ! only +RSD scan
	if(scan1) then
		print *, 'Warning!!!! Scan of only with RSD!'
		do i = 1, num_omw
			om = om_w_list(1,i)
			w  = om_w_list(2,i)

			if(gbtp) &
				write(*,'(A,i8,1x,f8.4,1x,f8.4,1x,f7.3,A)') '  Step, om, w, ratio = ', i, real(om), real(w), real(i/(num_omw+0.0)), '...'

			cs%has_RSD = .true.		
			call gf_mldprho_chi2s(om, w, h_dft, cs, rhoRSD_chisqlist, calc_comvr = .true.)		
			write(2,'(<2+cs%numdrop>(e14.7,1x))') om, w, rhoRSD_chisqlist(1:cs%numdrop)
				
			if(gbtp) then
				write(*,'(A,2f10.4,A)') '  Mean d_sep / r_sm_sphere = ', &
					(gbtotvol/num_halo)**(1.0/3.0), est_sm_sphe_r(gbtotvol, num_halo, cs%smnum), '  ...'
			endif
		
			call cpu_time(time2)
			if(gbtp) then
				write(*,'(A,f10.4)') '  Time used in this step: ', time2-time1
				write(*,'(A,<cs%numdrop>(f12.7,1x))') '  Chisqs of rho (RSD): ', rhoRSD_chisqlist(1:cs%numdrop)
				time1 = time2
				print *
			endif
		enddo
		stop
	endif
	
	scan2 = .true. ! only no RSD scan
	if(scan2) then
		print *, 'Warning!!!! Scan of only no RSD!'
		do i = 1, num_omw
			om = om_w_list(1,i)
			w  = om_w_list(2,i)

			if(gbtp) &
				write(*,'(A,i8,1x,f8.4,1x,f8.4,1x,f7.3,A)') '  Step, om, w, ratio = ', i, real(om), real(w), real(i/(num_omw+0.0)), '...'
		
			cs%has_RSD = .false.
			call gf_mldprho_chi2s(om, w, h_dft, cs, rho_chisqlist, calc_comvr = .true.)
			write(1,'(<2+cs%numdrop>(e14.7,1x))') om, w, rho_chisqlist(1:cs%numdrop)

			if(gbtp) then
				write(*,'(A,2f10.4,A)') '  Mean d_sep / r_sm_sphere = ', &
					(gbtotvol/num_halo)**(1.0/3.0), est_sm_sphe_r(gbtotvol, num_halo, cs%smnum), '  ...'
			endif
		
			call cpu_time(time2)
			if(gbtp) then
				write(*,'(A,f10.4)') '  Time used in this step: ', time2-time1
				write(*,'(A,<cs%numdrop>(f12.7,1x))') '  Chisqs of rho:       ', rho_chisqlist(1:cs%numdrop)
				time1 = time2
				print *
			endif
		enddo
		stop
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	print *, 'Scanning with both RSD, no RSD...'
	do i = 1, num_omw
		om = om_w_list(1,i)
		w  = om_w_list(2,i)

		if(gbtp) &
			write(*,'(A,i8,1x,f8.4,1x,f8.4,1x,f7.3,A)') '  Step, om, w, ratio = ', i, real(om), real(w), real(i/(num_omw+0.0)), '...'
		
		cs%has_RSD = .false.
		
!		call gradient_chisqs(om, w, h_dft, cs, rho_chisqlist, delta_chisqlist, ndelta_chisqlist, calc_comvr = .true.)
		call gf_mldprho_chi2s(om, w, h_dft, cs, rho_chisqlist, calc_comvr = .true.)

		write(1,'(<2+cs%numdrop>(e14.7,1x))') om, w, rho_chisqlist(1:cs%numdrop)
!		write(3,'(<2+cs%numdrop>(e14.7,1x))') om, w, delta_chisqlist(1:cs%numdrop)
!		write(5,'(<2+cs%numdrop>(e14.7,1x))') om, w, ndelta_chisqlist(1:cs%numdrop)
		
		cs%has_RSD = .true.
		
!		call gradient_chisqs(om, w, h_dft, cs, rhoRSD_chisqlist, deltaRSD_chisqlist, ndeltaRSD_chisqlist, calc_comvr = .false.)
		call gf_mldprho_chi2s(om, w, h_dft, cs, rhoRSD_chisqlist, calc_comvr = .false.)		
		
		if(gbtp) then
			write(*,'(A,2f10.4,A)') '  Mean d_sep / r_sm_sphere = ', &
				(gbtotvol/num_halo)**(1.0/3.0), est_sm_sphe_r(gbtotvol, num_halo, cs%smnum), '  ...'
		endif
			
		write(2,'(<2+cs%numdrop>(e14.7,1x))') om, w, rhoRSD_chisqlist(1:cs%numdrop)
!		write(4,'(<2+cs%numdrop>(e14.7,1x))') om, w, deltaRSD_chisqlist(1:cs%numdrop)
!		write(6,'(<2+cs%numdrop>(e14.7,1x))') om, w, ndeltaRSD_chisqlist(1:cs%numdrop)	
		
		call cpu_time(time2)
		if(gbtp) then
			write(*,'(A,f10.4)') '  Time used in this step: ', time2-time1
			write(*,'(A,<cs%numdrop>(f12.7,1x))') '  Chisqs of rho:       ', rho_chisqlist(1:cs%numdrop)
			write(*,'(A,<cs%numdrop>(f12.7,1x))') '  Chisqs of rho (RSD): ', rhoRSD_chisqlist(1:cs%numdrop)
!			write(*,'(A,<6*cs%numdrop>(f12.7,1x))') '  Chisqs: ', rho_chisqlist(1:cs%numdrop), delta_chisqlist(1:cs%numdrop), &
!				ndelta_chisqlist(1:cs%numdrop), rhoRSD_chisqlist(1:cs%numdrop), ndeltaRSD_chisqlist(1:cs%numdrop)
			time1 = time2
			print *
		else
			if(time2 - time1 > 30) then
				write(*,'(4x,A,i8,1x,f12.7,1x,f6.3,1x,f7.3,A)') ' step, om, w, ratio = ', i, real(om), real(w), real(i/(n+0.0)), '...'
!			write(*,'(4x,<2+num_nbins>(f8.2,1x))') om, w, rho_chisqlist(1:num_nbins)
!			write(*,'(4x,<2+num_nbins>(f8.2,1x))') om, w, delta_chisqlist(1:num_nbins)
!			write(*,'(4x,<2+num_nbins>(f8.2,1x))') om, w, ndelta_chisqlist(1:num_nbins)
				time1 = time2
			endif
		endif


	enddo
	
	call cpu_time(timeend)
	
	print *, 'Total time:     ',  timeend - timebegin
	print *, 'Total time / n: ', (timeend - timebegin) / (n+0.0)

	close(1); close(2); close(3); close(4); close(5); close(6)
end program ap_main
