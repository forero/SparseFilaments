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

!use mpi
use ap_chisq

	implicit none

	type(chisq_settings) :: cs
	character(len=char_len) :: inputfile, outputname, &
		rhofile, rhoRSDfile, deltafile, deltaRSDfile, ndeltafile, ndeltaRSDfile, str, tmpstr
	integer :: i, j, n, num_omw
	real(dl), allocatable :: om_w_list(:,:), rho_chisqlist(:), delta_chisqlist(:), ndelta_chisqlist(:),  &
			rhoRSD_chisqlist(:),deltaRSD_chisqlist(:),ndeltaRSD_chisqlist(:)
	real(dl) :: ommin, ommax, wmin, wmax, dom, dw, om, w, time1, time2, timebegin, timeend
	logical :: calc_om, calc_w

!-----------------------------------------
!	Basic settings
!----!			  	!---------
	cs%smnum 		= 10
	cs%num_in_x 		= 40
	cs%check_boundary 	= .true.
	cs%print_info		= .false.
	cs%cb_adjust_ratio	= 1.0_dl
	cs%remov_dist_ratio	= 1.0_dl
	cs%use_num_density	= .true.
	cs%nbins_rhoav		= 10
	cs%use_intpl_rho	= .true.
	gb_num_changenuminx	= 15
	gb_amp_changenuminx	= 0.04_dl
!----!			  	!---------
!-----------------------------------------

	write(*,*) 'Settings of the programm:'
	write(*,*) '    smnum             = ', cs%smnum
	write(*,*) '    num_in_x          = ', cs%num_in_x
	write(*,*) '    check_boundary    = ', cs%check_boundary
	write(*,*) '    print_info        = ', cs%print_info
	write(*,*) '    cb_adjust_ratio   = ', real(cs%cb_adjust_ratio)
	write(*,*) '    remov_dist_ratio  = ', real(cs%remov_dist_ratio)
	write(*,*) '    use_num_density   = ', cs%use_num_density
	write(*,*) '    nbins_rhoav       = ', cs%nbins_rhoav
	write(*,*) '    use_intpl_rho     = ', cs%use_intpl_rho
	write(*,'(A,i4)')    '     # of cellchange   : ', gb_num_changenuminx
	write(*,'(A,e14.7)') '     amp of cellchange : ', gb_amp_changenuminx	

!-----------------------------------------
!	Settings of droppings...
!----!			  	!---------
	cs%numdrop		= 20
!----!			  	!---------	
!-----------------------------------------

	allocate(cs%dropval(cs%numdrop),  cs%lowdropvalratio(cs%numdrop),  cs%highdropvalratio(cs%numdrop), &
		 cs%dropdval(cs%numdrop), cs%lowdropdvalratio(cs%numdrop), cs%highdropdvalratio(cs%numdrop))	
		 
	do i = 1, cs%numdrop
		cs%dropval(i) 		= .true. 
		cs%lowdropvalratio(i)	= 0.00_dl
		cs%dropdval(i)		= .true. 
		cs%lowdropdvalratio(i)	= 0.01_dl
	enddo

	cs%highdropvalratio(1:cs%numdrop) = (/0.3_dl, 0.4_dl, 0.45_dl, 0.5_dl, 0.55_dl, &
		0.6_dl, 0.65_dl, 0.7_dl, 0.75_dl, 0.8_dl, 0.82_dl, 0.84_dl, 0.86_dl, 0.88_dl, 0.9_dl, 0.92_dl, 0.94_dl, 0.96_dl, 0.97_dl, 0.98_dl/)
	cs%highdropdvalratio(1:cs%numdrop) = cs%highdropvalratio(1:cs%numdrop)

	print *
	write(*,*) '    num of drop       = ', cs%numdrop
	do i = 1, cs%numdrop
		write(*,*) 'Droppint pixels settings ', i,':'
		write(*,*) '    Density:            ', cs%dropval(i), real(cs%lowdropvalratio(i)),  real(cs%highdropvalratio(i))
		write(*,*) '    Density Gradient:   ', cs%dropdval(i), real(cs%lowdropdvalratio(i)), real(cs%highdropdvalratio(i))
	enddo

!-----------------------------------------
!	Settings of om/w
!----!			  	!---------
	num_omw = 48
	calc_om = .true.
	ommin 	= 0.22 
	ommax 	= 0.3
	calc_w  = .false.
	wmin  	= -1.2_dl 
	wmax 	= -0.85_dl	
!----!			  	!---------
!-----------------------------------------
	
	print *
	write(*,*) '    num of om/w       = ', num_omw
	
	if(calc_om) then
		outputname 	= 'om__smnum__'
		write(str,*) cs%smnum
		outputname 	= outputname//trim(adjustl(str))//'numpixels__'
		write(str,*) cs%num_in_x
		outputname 	= outputname//trim(adjustl(str))//'cubic'
		
		! read in the om_w_list
		allocate(om_w_list(2,num_omw))
		print * 
		print *, 'List of om/w:'
		dom 	= (ommax-ommin) / (num_omw-1.0_dl)
		do i = 1, num_omw
			om_w_list(1,i) = ommin + dom * (i-1)
			om_w_list(2,i) = w_dft
			if(gbtp) print*, '  omegam/w = ', real(om_w_list(1:2,i))
		enddo

		print *
		if (.not. gb_chisq_initied) then
			call cosmo_funs_init()
			call read_in_halo_data()
			call init_halo_info()
			gb_chisq_initied = .true.
		endif
	
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
		print *, 'Total time / n: ', (timeend - timebegin) / (num_omw+0.0)

		close(1); close(2); close(3); close(4); close(5); close(6)
		deallocate(om_w_list)
	endif
	
	
	if(calc_w) then
		print *
		outputname 	= 'w__smnum__'
		write(str,*) cs%smnum
		outputname 	= outputname//trim(adjustl(str))//'numpixels__'
		write(str,*) cs%num_in_x
		outputname 	= outputname//trim(adjustl(str))//'cubic'
		
		! read in the om_w_list
		allocate(om_w_list(2,num_omw))
		print *, 'List of om/w:'
		dw 	= (wmax-wmin) / (num_omw-1.0_dl)
		do i = 1, num_omw
			om_w_list(1,i) = om_dft
			om_w_list(2,i) = wmin + dw * (i-1)
			if(gbtp) print*, '  omegam/w = ', real(om_w_list(1:2,i))
		enddo

		print *
		if (.not. gb_chisq_initied) then
			call cosmo_funs_init()
			call read_in_halo_data()
			call init_halo_info()
			gb_chisq_initied = .true.
		endif
	
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
		print *, 'Total time / n: ', (timeend - timebegin) / (num_omw+0.0)

		close(1); close(2); close(3); close(4); close(5); close(6)
		deallocate(om_w_list)
	endif
end program ap_main
