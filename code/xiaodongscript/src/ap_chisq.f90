
module ap_chisq

use ap_grad_fields
implicit none

	logical :: gb_chisq_initied = .false.
	integer, public :: gb_chisq_method = 3
	integer, public :: gb_num_changenuminx = 10
	real(dl), public :: gb_amp_changenuminx = 0.1_dl
	! 1: (<abs(mu)>-0.5)**2.0 / var(abs(mu))
	! 2: histogram, nbins.
contains


  !------------------------------------------
  ! calculating chisq based
  !------------------------------------------
	  subroutine random_chisqs(np, ndata, chisqlist, chisqmean, chisqvar, nbins, printinfo)
	  	! DUMMy
	  	integer :: np, ndata, nbins
	  	real(dl) :: chisqlist(ndata), chisqmean, chisqabsmean, chisqvar
	  	logical :: printinfo 
	  	! LOCAL
	  	integer :: i, idata, j, n
	  	real(dl) :: mudata(np), x
	  	
	  	if(printinfo) print *, 'Totally ', ndata, 'steps...'
	  	call random_seed()
		do idata = 1, ndata
		  	do i = 1, np
		  		call random_number(x)
		  		mudata(i) = x + x - 1.0
		  	enddo
		  	chisqlist(idata) = 0
!		  	j = 0
!		  	do n = 2, 5
!			  	chisqlist(idata) = chisqlist(idata) + chisq_of_mu_data(mudata, n, nbins)
!			 	j = j + 1
!			enddo
!			chisqlist(idata) = chisqlist(idata) / (j+0.0)
!			chisqlist(idata) = chisq_of_mu_data(mudata, nbins)
!			chisqlist(idata) = chisq_of_mu_data2(mudata)	
			chisqlist(idata) = chisq_of_mu_data_shift(mudata, nbins, 3, .false.)			
		  	if(mod(idata,max(ndata/20,1)).eq.0.and.printinfo) then
			  	print *, '  step, chisq = ', idata, chisqlist(idata)
			endif
		  enddo
		  call get_mean_var(chisqlist, chisqmean, chisqvar)
		  if(printinfo) then
		  	print *, 'Mean, var, sqrt(var) = ', real(chisqmean), real(chisqvar), real(sqrt(chisqvar))
		  endif
	  end subroutine random_chisqs

  !------------------------------------------
  ! calculating chisq based
  !------------------------------------------
	  subroutine gf_mldprho_chi2s(omegam, w, h, cs, rho_chisqlist, rho_dfchisqlist, numchisq, calc_comvr)
		! dummy
		real(dl), intent(in) :: omegam, w, h
		integer, intent(in) :: numchisq
		type(chisq_settings), intent(inout) :: cs
		real(dl), intent(out) :: rho_chisqlist(numchisq), rho_dfchisqlist(numchisq)
		logical, intent(in) :: calc_comvr 
		! local
		integer :: i, j, RSD, AP
		real(dl), allocatable :: changenuminxlist(:), tmpchisqlist(:), tmpdfchisqlist(:)

		if (size(rho_chisqlist).ne.cs%numdrop) then
			print *, 'ERROR (gf_mldprho_chi2s)! size of rho_chisqlist mismatches with ', cs%numdrop
			stop
		endif

		if (.not. gb_chisq_initied) then
			call cosmo_funs_init()
			call read_in_halo_data()
			call init_halo_info()
			gb_chisq_initied = .true.
		endif

		! take the first cosmology as the input cosmology
		AP = 1
		gb_omegam 	= omegam
		gb_w 		= w
		gb_h 		= h
		if(calc_comvr) then
			if (cs%print_info) then
				print *, 'Estimating cosmology omegam/w = ', real(omegam), real(w)
			endif
			call de_calc_comovr()
			do i = 1, num_halo
				halo_info(i)%r_AP(AP) = de_get_comovr(halo_info(i)%z_real)
				halo_info(i)%r_AP_RSD(AP) = de_get_comovr(halo_info(i)%z_obs)
			enddo
		endif

		if(cs%has_RSD) then
			RSD = 1
		else
			RSD = 0
		endif

!		cs%print_info = .true.
		
		!get the chisqs
		allocate(changenuminxlist(2*gb_num_changenuminx+1), tmpchisqlist(cs%numdrop), tmpdfchisqlist(cs%numdrop))
		changenuminxlist(1) = 0.0_dl
		do i = 1, gb_num_changenuminx
			changenuminxlist(i+1) = gb_amp_changenuminx*i/(gb_num_changenuminx+0.0)
			changenuminxlist(2*gb_num_changenuminx-i+2) = -1.0_dl *  gb_amp_changenuminx*i/(gb_num_changenuminx+0.0)
		enddo
!			if(gbtp) then
!				print *, 'Using changenuminxlist: ', real(changenuminxlist)
!			endif
		rho_chisqlist = 0.0_dl
		do i = 1, size(changenuminxlist)
!			print *, 'Calculating many chisqs (with pixels rescaling): i = ', i
			call gd_mldprho_chi2s(RSD, AP, cs, changenuminxlist(i), tmpchisqlist, tmpdfchisqlist, cs%numdrop)
			if(gbtp) then
!				write(*,'(A,f9.4,A,<cs%numdrop>(f13.7,1x))'), '   chisqs at changenuminx = ', real(changenuminxlist(i)), ':', &
!					 real(tmpchisqlist(1:cs%numdrop))
			endif
			do j = 1, cs%numdrop
				rho_chisqlist(j) = rho_chisqlist(j) + tmpchisqlist(j) / (size(changenuminxlist) + 0.0)
				rho_dfchisqlist(j) = rho_dfchisqlist(j) + tmpdfchisqlist(j) / (size(changenuminxlist) + 0.0)				
			enddo
		enddo
		deallocate(changenuminxlist, tmpchisqlist, tmpdfchisqlist)
	end subroutine gf_mldprho_chi2s

  !------------------------------------------
  ! calculating chisq based
  !------------------------------------------
	  subroutine gf_mdwei_chi2s(omegam, w, h, cs, rho_chisqlist, calc_comvr)
		! dummy
		real(dl), intent(in) :: omegam, w, h
		type(chisq_settings), intent(inout) :: cs
		real(dl), intent(out) :: rho_chisqlist(gb_numwei)
		logical, intent(in) :: calc_comvr 
		! local
		integer :: i, j, RSD, AP, gb_numwei
		real(dl), allocatable :: changenuminxlist(:), tmpchisqlist(:)

		if (size(rho_chisqlist).ne.cs%numdrop) then
			print *, 'ERROR (gf_mldprho_chi2s)! size of rho_chisqlist mismatches with ', cs%numdrop
			stop
		endif

		if (.not. gb_chisq_initied) then
			call cosmo_funs_init()
			call read_in_halo_data()
			call init_halo_info()
			gb_chisq_initied = .true.
		endif

		! take the first cosmology as the input cosmology
		AP = 1
		gb_omegam 	= omegam
		gb_w 		= w
		gb_h 		= h
		if(calc_comvr) then
			if (cs%print_info) then
				print *, 'Estimating cosmology omegam/w = ', real(omegam), real(w)
			endif
			call de_calc_comovr()
			do i = 1, num_halo
				halo_info(i)%r_AP(AP) = de_get_comovr(halo_info(i)%z_real)
				halo_info(i)%r_AP_RSD(AP) = de_get_comovr(halo_info(i)%z_obs)
			enddo
		endif

		if(cs%has_RSD) then
			RSD = 1
		else
			RSD = 0
		endif

!		cs%print_info = .true.
		
		!get the chisqs
		allocate(changenuminxlist(2*gb_num_changenuminx+1), tmpchisqlist(gb_numwei))
		changenuminxlist(1) = 0.0_dl
		do i = 1, gb_num_changenuminx
			changenuminxlist(i+1) = gb_amp_changenuminx*i/(gb_num_changenuminx+0.0)
			changenuminxlist(2*gb_num_changenuminx-i+2) = -1.0_dl *  gb_amp_changenuminx*i/(gb_num_changenuminx+0.0)
		enddo
!			if(gbtp) then
!				print *, 'Using changenuminxlist: ', real(changenuminxlist)
!			endif
		rho_chisqlist = 0.0_dl
		do i = 1, size(changenuminxlist)
!			print *, 'Calculating many chisqs (with pixels rescaling): i = ', i
			call gd_mdwei_chi2s(RSD, AP, cs, changenuminxlist(i), tmpchisqlist, gb_numwei)
			if(gbtp) then
!				write(*,'(A,f9.4,A,<cs%numdrop>(f13.7,1x))'), '   chisqs at changenuminx = ', real(changenuminxlist(i)), ':', &
!					 real(tmpchisqlist(1:cs%numdrop))
			endif
			do j = 1, gb_numwei
				rho_chisqlist(j) = rho_chisqlist(j) + tmpchisqlist(j) / (size(changenuminxlist) + 0.0)
			enddo
		enddo
		deallocate(changenuminxlist, tmpchisqlist)
	end subroutine gf_mdwei_chi2s


end module ap_chisq
