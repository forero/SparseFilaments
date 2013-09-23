
program ap_main_test

use mpi
use ap_tools
use ap_chisq

	implicit none

	integer, parameter :: np = 920, ndata = 2000
	real(dl) :: chisqlist(ndata), chisqmean, chisqvar
	real(dl) :: mudata(ndata), chisq
	integer :: i

!	do i = 2, 20
!		print *, 'nbins = ', i

!	i = 5
	call random_chisqs(np, ndata, chisqlist, chisqmean, chisqvar, 5, .true.)
!	enddo

!	call random_seed()
!	call random_number(mudata)
!	do i = 1, ndata
!		mudata(i) = mudata(i)*2 - 1.0
!	enddo
!	chisq = chisq_of_mu_data_shift(mudata, 5, 100, .true.)
end program ap_main_test
