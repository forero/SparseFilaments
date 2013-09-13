
program ap_main_test

use mpi
use ap_tools
use ap_chisq

	implicit none

	integer, parameter :: nmean = 100, np = 92000, ndata = 2000, nbins = 5
	real(dl) :: chisqlist(ndata), chisqmean, chisqvar, meanlist(nmean), sqvarlist(nmean)
	integer :: i
	real(dl) :: meanmean, varmean, meansqrtvar

	call random_chisqs(np, ndata, chisqlist, chisqmean, chisqvar, nbins, .true.)
	
	goto 100

	do i = 1, nmean
		call random_chisqs(np, ndata, chisqlist, chisqmean, chisqvar, nbins, .false.)
	  	if(mod(i,max(nmean/5,1)).eq.0) then
			print*, '  Step ', i, ':'
			print*, '    Mean, var, sqrt(var) = ', real(chisqmean), real(chisqvar), real(sqrt(chisqvar))
		endif
		meanlist(i) = chisqmean
		sqvarlist(i) = sqrt(chisqvar)
	enddo
	
	call get_mean_var(meanlist, meanmean, varmean)
	print *, 'Meanmean, varmean, sqrt(varmean) =   ', real(meanmean), real(varmean), real(sqrt(varmean))
	call get_mean_var(sqvarlist,meansqrtvar)
	print *, 'Mean of sqrt(chisqvar)             : ', real(meansqrtvar)
	
100	continue	
end program ap_main_test
