


!####################################
!This module does smooth
!####################################
module ap_settings_init
use ap_cosmo_funs

	implicit none
	
!--------------------------------------------------------------------
!--------------------------------------------------------------------	
!--  settings of the programme   ------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------

	!----------------------------------
	! ( 1 ). Halo dat file
	!----------------------------------
	
		!character(char_len), parameter :: halo_data_file = '../../../data/input/MD_fullhalos_z0_max160w.dat'
		!character(char_len), parameter :: halo_data_file = '../../data/input/MD_fullhalos_z0_1000-1600_shell.dat'
		character(char_len), parameter :: halo_data_file = '../../data/input/MD_fullhalos_z0_100w.dat'
		!character(char_len), parameter :: halo_data_file = '../../data/input/MD_fullhalos_z0_330w.dat'
	
	!----------------------------------
	! ( 2 ). Name of the project.	
	!----------------------------------
	
		character(char_len), parameter :: halo_data_name = 'MD_halos_160w_1000-1600-shell'
	
		character(char_len), parameter :: output_dir = '../../data/output/'//trim(adjustl(halo_data_name))//'/'
		
	!----------------------------------
	! ( 3 ). Real cosmological parameters	
	!----------------------------------	
		
		real(dl), parameter :: om_dft = 0.27_dl, h_dft = 0.70_dl, w_dft = -1.0_dl
	
	!----------------------------------
	! ( 4 ). Assumed cosmologies.
	!----------------------------------	
		
		integer, parameter :: num_assumed_cosmo = 3
		
		real(dl), parameter :: assumed_cosmo(3,num_assumed_cosmo) = &
			(/ 0.5_dl,  w_dft, h_dft, & 
			   0.75_dl, w_dft, h_dft, &
			   1.0_dl,  w_dft, h_dft   /)
			   
		character(char_len), parameter :: ascos_name1 = 'om_0p5'
		character(char_len), parameter :: ascos_name2 = 'om_0p75'
		character(char_len), parameter :: ascos_name3 = 'om_0p1'

		character(char_len), parameter :: assumed_cosmo_name(num_assumed_cosmo) = (/ ascos_name1, ascos_name2, ascos_name3 /)
	

!--------------------------------------------------------------------
!--------------------------------------------------------------------	
!--  End of settings of the programme   -----------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
	
	type :: point
		real(dl) :: mass, x,y,z,vx,vy,vz,vlos
		real(dl) :: r,theta, phi
		real(dl) :: z_real, z_obs,  r_RSD
		real(dl) :: r_AP(num_assumed_cosmo), r_AP_RSD(num_assumed_cosmo)
	end type

	real(dl), allocatable :: halo_data(:,:)
	integer :: num_halo
	
	type(point), allocatable :: halo_info(:)

contains

!--------------------------------------------------------------------
!--------------------------------------------------------------------	
!--  settings of the programme   ------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
	!----------------------------------
	! ( 5 ). Assumed cosmologies.
	!----------------------------------
  		!------------------------------------------
		! check whether a point has boundary effect
		!------------------------------------------
		logical function has_boundary_effect(x,y,z,max_dist,rmin,rmax,cb_adjust_ratio)
			real(dl) :: x,y,z,max_dist,rmin,rmax,cb_adjust_ratio,r,dr
			dr = max_dist*cb_adjust_ratio
			r = sqrt(x*x+y*y+z*z)
			if (r<rmin+dr.or.r>rmax-dr) then
				has_boundary_effect = .true.; return
			endif
			if(x<dr.or.y<dr.or.z<dr) then
				has_boundary_effect = .true.; return
			endif
			has_boundary_effect=.false.
		end function has_boundary_effect
		
	!----------------------------------
	! ( 6 ). Volume function.
	!----------------------------------		
		!------------------------------------------
		! volume between r1 and r2
		!------------------------------------------
		real(dl) function vol_fun(r1, r2)
			real(dl) :: r1, r2
			vol_fun = (4.0*const_pi/3.0) * (r2**3.0-r1**3.0) / 8.0
		end function vol_fun

	
!--------------------------------------------------------------------
!--------------------------------------------------------------------	
!--  End of settings of the programme   -----------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------

  !------------------------------------------
  ! read in the halo data
  !------------------------------------------
	subroutine read_in_halo_data()
		real(dl) :: shell_rmin = 500.0, shell_rmax = 800.0
		real(dl) :: rsq1, rsq2, x,y,z,rsq
		real(dl), allocatable :: tmp(:,:)
		integer :: i, j, n, num

		print *, 'Read in the halo data from ', trim(adjustl(halo_data_file)), '...'		

		! directly read in from ready-data file.
		if(allocated(halo_data)) &
			deallocate(halo_data)
		call read_in_revfmt(halo_data_file, 7, num_halo, halo_data)
		return
		! no need to do shell by ourselves.
		
		rsq1 = shell_rmin*shell_rmin
		rsq2 = shell_rmax * shell_rmax
		call read_in(halo_data_file, 8, n, tmp)
		num = 0
		do i = 1, n
			x = tmp(i,1)
			y = tmp(i,2)
			z = tmp(i,3)
			rsq = x*x + y*y + z*z	
			if(rsq > rsq1 .and. rsq < rsq2) num=num+1
		enddo
		num_halo = num
		if(allocated(halo_data)) &
			deallocate(halo_data)		
		allocate(halo_data(7,num_halo))
		j = 1
		do i = 1, n
			x = tmp(i,1)
			y = tmp(i,2)
			z = tmp(i,3)
			rsq = x*x + y*y + z*z	
			if(rsq > rsq1 .and. rsq < rsq2) then
				halo_data(1:3,j) = 2.0_dl*tmp(i,1:3)
				halo_data(4:6,j) = tmp(i,4:6)
				halo_data(7,j) = tmp(i,7)
				j = j + 1
			endif
		enddo
		deallocate(tmp)
	end subroutine read_in_halo_data


  !------------------------------------------
  ! initializing the halo_info array
  !------------------------------------------
	subroutine init_halo_info()
		integer :: i, i_cosmo
		
		if(.not. allocated(halo_data)) then
			print *, 'ERROR! halo_data not allocated. Unable to initialize halo_info!'
			stop
		endif
		
		print *, 'Initializing the halo_info...'
		print *, '	num of halo = ', num_halo
		! initialize the comoving r array
		gb_omegam = om_dft; gb_w = w_dft; gb_h = h_dft
		call de_calc_comovr()
		if(allocated(halo_info)) &
			deallocate(halo_info)
		allocate(halo_info(num_halo))
		do i = 1, num_halo
			halo_info(i)%x = halo_data(1,i)
			halo_info(i)%y = halo_data(2,i)
			halo_info(i)%z = halo_data(3,i)
			halo_info(i)%vx = halo_data(4,i)
			halo_info(i)%vy = halo_data(5,i)
			halo_info(i)%vz = halo_data(6,i)
			halo_info(i)%mass = halo_data(7,i)
			call getSC(halo_info(i)%x,halo_info(i)%y,halo_info(i)%z,halo_info(i)%r,halo_info(i)%theta,halo_info(i)%phi)
			halo_info(i)%vlos = (halo_info(i)%x*halo_info(i)%vx + halo_info(i)%y*halo_info(i)%vy + halo_info(i)%z*halo_info(i)%vz) / halo_info(i)%r 
			! redshifts
			halo_info(i)%z_real = get_z(halo_info(i)%r)
			halo_info(i)%z_obs = halo_info(i)%z_real + halo_info(i)%vlos / const_c
			halo_info(i)%r_RSD = de_get_comovr(halo_info(i)%z_obs)
		enddo
		
		do i_cosmo = 1, num_assumed_cosmo
			gb_omegam 	= assumed_cosmo(1,i_cosmo)
			gb_w 		= assumed_cosmo(2,i_cosmo)
			gb_h 		= assumed_cosmo(3,i_cosmo)
			call de_calc_comovr()
			do i = 1, num_halo
				halo_info(i)%r_AP(i_cosmo) = de_get_comovr(halo_info(i)%z_real)
				halo_info(i)%r_AP_RSD(i_cosmo) = de_get_comovr(halo_info(i)%z_obs)
			enddo
		enddo
		
		print *, 'Creating outpu_dir = ', trim(adjustl(output_dir)), '...'
		call system('mkdir -p '//trim(adjustl(output_dir)))
	end subroutine init_halo_info

  !------------------------------------------
  ! self-introduction of a point
  !------------------------------------------	
	subroutine point_speak(p)
		integer :: i
		type(point) :: p
		print *, 'Information of the point:'
		print *, '  x,y,z = ', p%x, p%y, p%z
		print *, '  vx,vy,vz = ', p%vx, p%vy, p%vz
		print *, '  r,theta,phi = ', p%r, p%theta, p%phi
		print *, '  los speed =', p%vlos
		print *, '  cosmo redshift = ', p%z_real
		print *, '  apparent redshift = ', p%z_obs
		print *, '  apparent comov r  = ', p%r_RSD
		do i = 1, num_assumed_cosmo
			print *, '  Info for cosmology ', i, ':'
			print *, '    Inferred r (no   RSD) = ', p%r_AP(i)  
			print *, '    Inferred r (with RSD) = ', p%r_AP_RSD(i)
		enddo
	end subroutine point_speak

  !------------------------------------------
  ! function estimating mean radius
  !  of the smoothing sphere
  !------------------------------------------	
	real(dl) function est_sm_sphe_r(tot_vol, tot_num, smooth_num)
		real(dl) :: tot_vol
		integer :: tot_num, smooth_num
		est_sm_sphe_r = (3.0 * (tot_vol/(tot_num+0.0)) * smooth_num) / (4.0*const_pi)
		est_sm_sphe_r = est_sm_sphe_r**(1.0_dl/3.0_dl)
	end function est_sm_sphe_r
		

end module ap_settings_init
	



