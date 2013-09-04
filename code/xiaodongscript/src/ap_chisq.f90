
module ap_chisq

use ap_grad_fields
implicit none

	logical :: gb_chisq_initied = .false.
	type :: chisq_settings
		integer :: smnum, num_in_x
		logical :: check_boundary = .true.
		logical :: print_info = .false.
		real(dl) :: cb_adjust_ratio = 1.0_dl, remov_dist_ratio = 1.0_dl
		logical :: use_num_density = .true.
		integer :: nbins_rhoav = 10
		logical :: use_intpl_rho = .true.
		logical :: has_RSD = .true.
		integer, allocatable :: nbins_list(:)
	end type
	

contains


  !------------------------------------------
  ! calculating chisq based
  !------------------------------------------
	  subroutine gradient_chisqs(omegam, w, h, cs, rho_chisqlist, delta_chisqlist, ndelta_chisqlist)
		real(dl) :: omegam, w, h
		type(chisq_settings) :: cs
		integer :: i, j, n, RSD, AP, nbins
		real(dl), allocatable :: drho_mu_data(:), ddelta_mu_data(:), dndelta_mu_data(:)
		real(dl), intent(out) :: rho_chisqlist(:), delta_chisqlist(:), ndelta_chisqlist(:)

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
		if (cs%print_info) then
			print *, 'Estimating cosmology om, w = ', omegam, w
		endif
		call de_calc_comovr()
		do i = 1, num_halo
			halo_info(i)%r_AP(AP) = de_get_comovr(halo_info(i)%z_real)
			halo_info(i)%r_AP_RSD(AP) = de_get_comovr(halo_info(i)%z_obs)
		enddo

		if(cs%has_RSD) then
			RSD = 1
		else
			RSD = 0
		endif

		! get the mu_data
		call grid_rho_delta_list(RSD, AP, cs%smnum, cs%num_in_x, cs%print_info, cs%check_boundary, cs%cb_adjust_ratio, cs%remov_dist_ratio, &
			cs%use_num_density, cs%nbins_rhoav, cs%use_intpl_rho, &
			drho_mu_data, ddelta_mu_data, dndelta_mu_data)
		
		!get the chisqs
		do i = 1, n
			nbins = cs%nbins_list(i)
			rho_chisqlist(i) = chisq_of_mu_data(drho_mu_data,nbins)
			delta_chisqlist(i) = chisq_of_mu_data(ddelta_mu_data,nbins)
			ndelta_chisqlist(i) = chisq_of_mu_data(dndelta_mu_data,nbins)
		enddo
	end subroutine gradient_chisqs

end module ap_chisq
			
