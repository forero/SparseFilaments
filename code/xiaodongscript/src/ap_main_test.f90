
program ap_main_test

use mpi
use ap_tools
use ap_chisq

	implicit none

	integer, parameter :: narray = 10000
	integer :: i, ix, iy, iz, di,istart,iend,nbclist(narray)
	
	ix = 0
	iy = 0
	iz = 0
	
	istart = 1
	
	do i = 1, 1000000
		do di = 0, 3
			call nbcell_list(ix,iy,iz,di,narray,istart,iend,nbclist)
			istart = iend+1
		enddo
	enddo
	
	print *, 'Finally, iend = ', iend

end program ap_main_test
