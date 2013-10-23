
program ap_main_test

!use mpi
use ap_tools
use ap_chisq
use ap_grids

	implicit none

	integer, parameter :: n = 10000, numchisq = 10000
	integer :: i
	real(dl) :: array(n), weight(n), chisq1(numchisq), chisq2(numchisq), mean, var

!	call random_seed()
	call cosmo_funs_init()
	call read_in_halo_data()
	call init_halo_info()
	call init_pixels(RSD=0, AP=0, rl_num_in_x=50.0_dl, print_info=.true., use_num_density=.true., do_xyzmassinit=.true.)
	stop

	weight = 1.0_dl
	do i = 1, n / 2
		weight(i) = 0.5_dl
	enddo
!	weight = weight * 10.0_dl

	do i = 1, numchisq	
		call random_number(array)
		chisq1(i) = chisq_of_mu_data2(array)
		chisq2(i) = weighted_chisq(array, weight, n)
	enddo
	
!	print *, chisq1
!	print *, chisq2
	
	call get_mean_var(chisq1, mean, var)
	print *, mean, var

	call get_mean_var(chisq2, mean, var)
	print *, mean, var	
end program ap_main_test
