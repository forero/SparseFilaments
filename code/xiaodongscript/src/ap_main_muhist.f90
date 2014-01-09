

program ap_main_muhist

!use mpi
use ap_tools
use ap_chisq
use ap_grids

implicit none

	integer :: numarg, i, i_LC
	character(len=char_len) :: str, inputfile, outputfile
	integer, parameter :: numbins = 10
	type(mustatype) :: mustat, NORSDmustat(numbins), RSDmustat(numbins)
	
	halo_data_file = '../../../data/input/HR3LC0_57w_600to1787.dat'
!	halo_data_file = '../../../data/input/HR3halo100w_pos1000.dat' ! Successful if use 100-150 distance cut, RSD 
	call cosmo_funs_init()
	call read_in_halo_data()
	call init_halo_info()	
	
	numarg = iargc()
	if(numarg .ne. 1) then
		print*
		print*, "ERROR! Number of arg must be 1."
		print*
		print*, "USAGE: "
		print*
		print*, "  ./post_proc inputfile"
		print*
!		stop
	endif

	call getarg(1,inputfile)
	open(unit=1,file=inputfile)
	read(1,*) mustat%numbins, i_LC, mindr, maxdr, omegam, w
	
	AP = 1
	gb_omegam 	= omegam
	gb_w 		= w
	gb_h 		= h
	if(calc_cos) then
		call de_calc_comovr()
		do i = 1, num_halo
			halo_info(i)%r_AP(AP) = de_get_comovr(halo_info(i)%z_real)
			halo_info(i)%r_AP_RSD(AP) = de_get_comovr(halo_info(i)%z_obs)
		enddo
	endif
	
	call tpCor(mindr,maxdr,mu1,mu2,countmu1,countmu2,diffmu,chisq,mustat)
	
	
	
	
