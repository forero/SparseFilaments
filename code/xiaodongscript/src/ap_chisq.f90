
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
	  subroutine gradient_chisqs(omegam, w, h, cs, rho_chisqlist, delta_chisqlist, ndelta_chisqlist, calc_comvr)
		real(dl) :: omegam, w, h
		type(chisq_settings) :: cs
		integer :: i, n, RSD, AP, nbins
		real(dl), allocatable :: drho_mu_data(:), ddelta_mu_data(:), dndelta_mu_data(:)
		real(dl), intent(out) :: rho_chisqlist(:), delta_chisqlist(:), ndelta_chisqlist(:)
		logical :: calc_comvr
		real(dl), allocatable :: changenuminxlist(:), tmpchisqlist(:,:)
		real(dl) :: tmp1, tmp2, tmp3

		if(.not.allocated(cs%nbins_list)) then
			print *, 'ERROR! nbins_list as a mandatory setting must be provided!'
			stop
		endif
		
		n = size(cs%nbins_list)
		if(size(rho_chisqlist).ne.n .or.size(delta_chisqlist).ne.n .or.size(ndelta_chisqlist).ne.n) then
			print *, 'ERROR! size of chisq_list must match with the size of cs%nbins_list: ', n
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
		if(gb_chisq_method .eq. 1) then
			call grid_rho_delta_list(RSD, AP, cs, drho_mu_data, ddelta_mu_data, dndelta_mu_data)
			rho_chisqlist(1) = chisq_of_mu_data2(drho_mu_data)
			delta_chisqlist(1) = chisq_of_mu_data2(ddelta_mu_data)
			ndelta_chisqlist(1) = chisq_of_mu_data2(dndelta_mu_data)
		elseif(gb_chisq_method .eq. 2) then
			call grid_rho_delta_list(RSD, AP, cs, drho_mu_data, ddelta_mu_data, dndelta_mu_data)
			do i = 1, n
				nbins = cs%nbins_list(i)
				rho_chisqlist(i) = chisq_of_mu_data(drho_mu_data,nbins)
				delta_chisqlist(i) = chisq_of_mu_data(ddelta_mu_data,nbins)
				ndelta_chisqlist(i) = chisq_of_mu_data(dndelta_mu_data,nbins)
			enddo
		elseif(gb_chisq_method .eq. 3) then
			allocate(changenuminxlist(2*gb_num_changenuminx+1), tmpchisqlist(3,2*gb_num_changenuminx+1))
			changenuminxlist(1) = 0.0_dl
			do i = 1, gb_num_changenuminx
				changenuminxlist(i+1) = gb_amp_changenuminx*i/(gb_num_changenuminx+0.0)
				changenuminxlist(2*gb_num_changenuminx-i+2) = -1.0_dl *  gb_amp_changenuminx*i/(gb_num_changenuminx+0.0)
			enddo
!			if(gbtp) then
!				print *, 'Using changenuminxlist: ', real(changenuminxlist)
!			endif
			rho_chisqlist = 0.0_dl
			delta_chisqlist = 0.0_dl
			ndelta_chisqlist = 0.0_dl
			do i = 1, size(changenuminxlist)
				call grid_rho_delta_list2(RSD, AP, cs, changenuminxlist(i), drho_mu_data, ddelta_mu_data, dndelta_mu_data)
				tmpchisqlist(1,i) = chisq_of_mu_data2(drho_mu_data)
				tmpchisqlist(2,i) = chisq_of_mu_data2(ddelta_mu_data)
				tmpchisqlist(3,i) = chisq_of_mu_data2(dndelta_mu_data)
				if(gbtp) then
					write(*,'(A,f7.4,A,3f11.7)'), '   chisqs at changenuminx = ', real(changenuminxlist(i)), ':', &
					 real(tmpchisqlist(1:3,i))
				endif
				rho_chisqlist(1) = rho_chisqlist(1) + tmpchisqlist(1,i) / (size(changenuminxlist) + 0.0)
				delta_chisqlist(1) = delta_chisqlist(1) + tmpchisqlist(2,i) / (size(changenuminxlist) + 0.0)
				ndelta_chisqlist(1) = ndelta_chisqlist(1) + tmpchisqlist(3,i) / (size(changenuminxlist) + 0.0)
			enddo
			!print *
			write(*,'(A,<2*gb_num_changenuminx+1>(f9.3,1x))')  '  List of chisqs based on rho:' , real(tmpchisqlist(1,1:size(changenuminxlist)))
			call get_mean_var(tmpchisqlist(1,1:size(changenuminxlist)), tmp1, tmp2)
			print *, 'Mean, var, sqrt(var) of chisq = ', real(tmp1), real(tmp2), real(sqrt(tmp2))
			!print *
			deallocate(changenuminxlist, tmpchisqlist)
		endif
	end subroutine gradient_chisqs

end module ap_chisq
