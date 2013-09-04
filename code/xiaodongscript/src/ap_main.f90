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
program ap_main

use mpi
use ap_chisq

	implicit none

	type(chisq_settings) :: cs
	character(len=char_len) :: inputfile, outputname, &
		rhofile, rhoRSDfile, deltafile, deltaRSDfile, ndeltafile, ndeltaRSDfile, str
	integer :: num_nbins
	integer*2 :: tmpint 
	integer :: i, n, num_om, num_w
	real(dl), allocatable :: rho_chisqlist(:), delta_chisqlist(:), ndelta_chisqlist(:), om_w_list(:,:)
	real(dl) :: ommin, ommax, wmin, wmax, dom, dw, om, w, time1, time2
	logical :: auto_num_in_x
	
	
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
		cs%remov_dist_ratio, cs%use_num_density, cs%nbins_rhoav, cs%use_intpl_rho, num_nbins


	read(1,*) gb_dropval, gb_dropvalratio(1:2), gb_dropdval, gb_dropdvalratio(1:2)
	allocate(cs%nbins_list(num_nbins))
	read(1,*) cs%nbins_list(1:num_nbins)		 
	write(*,*) 'Settings of the programm:'
	tmpint = cs%smnum
	write(*,*) '    smnum             = ', tmpint
	if(.not.auto_num_in_x) then
		tmpint = cs%num_in_x
		write(*,*) '    num_in_x          = ', tmpint
	endif
	write(*,*) '    auto_num_in_x     = ', auto_num_in_x
	write(*,*) '    check_boundary    = ', cs%check_boundary
	write(*,*) '    print_info        = ', cs%print_info
	write(*,*) '    cb_adjust_ratio   = ', real(cs%cb_adjust_ratio)
	write(*,*) '    remov_dist_ratio  = ', real(cs%remov_dist_ratio)
	write(*,*) '    use_num_density   = ', cs%use_num_density
	tmpint = cs%nbins_rhoav
	write(*,*) '    nbins_rhoav       = ', tmpint 
	write(*,*) '    use_intpl_rho     = ', cs%use_intpl_rho
	tmpint = num_nbins
	write(*,*) '    num_nbins         = ', tmpint
	write(*,'(A,<num_nbins>(i4,2x))') '     list of nbins     : ', cs%nbins_list	
	write(*,*) 'Droppint pixels settings: '
	write(*,*) '    Density:            ', gb_dropval, real(gb_dropvalratio(1:2))
	write(*,*) '    Density Gradient:   ', gb_dropdval, real(gb_dropdvalratio(1:2))
	n = 0
	do while(.true.)
		read(unit=1,end=100,fmt=*) str
		n = n + 1
	enddo
100	allocate(om_w_list(2,n))
	!print *, '# of om_w = ', n
	
	! read in the om_w_list
	close(1)
	open(unit=1,file=inputfile)
	read(1,*) str
	read(1,*) str
	read(1,*) str
	do i = 1, n
		read(1,*) om_w_list(1:2,i)
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

	

	str = outputname
	rhofile = trim(adjustl(str))//'_rho.txt'
	rhoRSDfile = trim(adjustl(str))//'_rhoRSD.txt'
	deltafile = trim(adjustl(str))//'_delta.txt'
	deltaRSDfile = trim(adjustl(str))//'_deltaRSD.txt'
	ndeltafile = trim(adjustl(str))//'_ndelta.txt'
	ndeltaRSDfile = trim(adjustl(str))//'_ndeltaRSD.txt'
	
	print *
	print *, 'Calculating chisqs...'	
	print *, 'Results saved in the following file:'
	print *, '     ', trim(adjustl(rhofile))
	print *, '     ', trim(adjustl(rhoRSDfile))
	print *, '     ', trim(adjustl(deltafile))
	print *, '     ', trim(adjustl(deltaRSDfile))
	print *, '     ', trim(adjustl(ndeltafile))
	print *, '     ', trim(adjustl(ndeltaRSDfile))

	open(unit=1,file=rhofile)
	open(unit=2,file=rhoRSDfile)
	open(unit=3,file=deltafile)
	open(unit=4,file=deltaRSDfile)
	open(unit=5,file=ndeltafile)
	open(unit=6,file=ndeltaRSDfile)
	
	allocate(rho_chisqlist(num_nbins),delta_chisqlist(num_nbins),ndelta_chisqlist(num_nbins))
		
	call cpu_time(time1)
	
	do i = 1, n
		om = om_w_list(1,i)
		w  = om_w_list(2,i)
		
		cs%has_RSD = .false.
		
		call gradient_chisqs(om, w, h_dft, cs, rho_chisqlist, delta_chisqlist, ndelta_chisqlist)

		call cpu_time(time2)
		if(time2 - time1 > 30) then
			write(*,'(4x,A,i8,1x,f6.2,1x,f6.2,1x,f6.3,A)') ' step, om, w, ratio = ', i, real(om), real(w), real(i/(n+0.0)), '...'
!			write(*,'(4x,<2+num_nbins>(f8.2,1x))') om, w, rho_chisqlist(1:num_nbins)
!			write(*,'(4x,<2+num_nbins>(f8.2,1x))') om, w, delta_chisqlist(1:num_nbins)
!			write(*,'(4x,<2+num_nbins>(f8.2,1x))') om, w, ndelta_chisqlist(1:num_nbins)
			time1 = time2
		endif
		write(1,'(<2+num_nbins>(e14.7,1x))') om, w, rho_chisqlist(1:num_nbins)
		write(3,'(<2+num_nbins>(e14.7,1x))') om, w, delta_chisqlist(1:num_nbins)
		write(5,'(<2+num_nbins>(e14.7,1x))') om, w, ndelta_chisqlist(1:num_nbins)
		
		cs%has_RSD = .true.
		
		
		call gradient_chisqs(om, w, h_dft, cs, rho_chisqlist, delta_chisqlist, ndelta_chisqlist)

				
		write(2,'(<2+num_nbins>(e14.7,1x))') om, w, rho_chisqlist(1:num_nbins)
		write(4,'(<2+num_nbins>(e14.7,1x))') om, w, delta_chisqlist(1:num_nbins)
		write(6,'(<2+num_nbins>(e14.7,1x))') om, w, ndelta_chisqlist(1:num_nbins)	
		
!		if(time2 - time1 > 10) then
!			write(*,*) ' step, om, w = ', i, real(om), real(w), '...'
!			write(*,'(4x,<2+num_nbins>(f8.2,1x))') om, w, rho_chisqlist(1:num_nbins)
!			write(*,'(4x,<2+num_nbins>(f8.2,1x))') om, w, delta_chisqlist(1:num_nbins)
!			write(*,'(4x,<2+num_nbins>(f8.2,1x))') om, w, ndelta_chisqlist(1:num_nbins)
!			time1 = time2
!		endif
	enddo
	

	close(1); close(2); close(3); close(4); close(5); close(6)
end program ap_main
