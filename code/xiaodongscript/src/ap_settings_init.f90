
! List of settings 

! ap_smooth.f90
!	logical, public :: fixgridrange = .false.
!	logical, public :: use_fixmd = .true.
!	logical, public :: use_seg = .false.
!	logical, public :: use_realrange = .false.


!####################################
!This module does smooth
!####################################
module ap_settings_init
use ap_cosmo_funs

	implicit none

	real(dl) :: gbgridxmin, gbgridxmax, gbgridymin, gbgridymax, gbgridzmin, gbgridzmax, gbgridrmin, gbgridrmax
	real(dl) :: gbxmin, gbxmax, gbymin, gbymax, gbzmin, gbzmax, gbrmin, gbrmax
	real(dl) :: gbrealxmin, gbrealxmax, gbrealymin, gbrealymax, gbrealzmin, gbrealzmax
	real(dl) :: gbvratio = 1.0_dl
	logical  :: local_pecu_v = .true.
	logical  :: dotsbe !TESTING
!--------------------------------------------------------------------
!--------------------------------------------------------------------	
!--  settings of the programme   ------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------

	!----------------------------------
	! ( 1 ). Halo dat file
	!----------------------------------
	
!		character(char_len), parameter :: halo_data_file = '../../../data/input/MD_fullhalos_z0_max160w.dat'
!		character(char_len), parameter :: halo_data_file = '../../data/input/MD_fullhalos_z0_1000-1600_shell.dat'
!		character(char_len), parameter :: halo_data_file = '../../data/input/MD_fullhalos_z0_100w.dat'
!		character(char_len), parameter :: halo_data_file = '../../data/input/HR3halo100w_pos1000.dat'		
!		character(char_len), parameter :: halo_data_file = '../../data/input/HR3halo9w_pos1000.dat'		
!		character(char_len), parameter :: halo_data_file = '../../data/input/HR3halo100w_pos1000_times2.dat'	
!		character(char_len), parameter :: halo_data_file = '../../data/input/MD_fullhalos_z0_20wBox.dat'
!		character(char_len), parameter :: halo_data_file = '../../data/input/MD_fullhalos_z0_20wBox_add800.dat'		
!		character(char_len), parameter :: halo_data_file = '../../data/input/MD_fullhalos_z0_100wBox_add800.dat'
!		character(char_len), parameter :: halo_data_file = '../../data/input/HR2_20w_add800.dat'			
!		character(char_len), parameter :: halo_data_file = '../../data/input/HR3halo_20wbox_add800.dat'
!		character(char_len), parameter :: halo_data_file = '../../data/input/HR3_71w_800to2300box.dat'
!		character(char_len), parameter :: halo_data_file = '../../data/input/HR3_71w_500to2000box.dat'
!		character(char_len), parameter :: halo_data_file = '../../data/input/HR3_71w_500to2000boxA.dat'		
!		character(char_len), parameter :: halo_data_file = '../../data/input/HR3_71w_500to2000boxB.dat'		
!		character(char_len), parameter :: halo_data_file = '../../data/input/HR3_170w_0to2000box.dat'
!		character(char_len), parameter :: halo_data_file = '../../data/input/HR3_500to2000shellB.dat'		
!		character(char_len), parameter :: halo_data_file = '../../data/input/HR3_500to2000shellA.dat'		

!	In order to check the boundary effect, we test a series of boxes...
!		character(char_len), parameter :: halo_data_file = '../../data/input/HR3_1p7Gpc_100w_boxA.dat'		
!		character(char_len), parameter :: halo_data_file = '../../data/input/HR3_2p7_1p35_1p35_box.dat'
!		character(char_len), parameter :: halo_data_file = '../../data/input/HR3_2p1_100w_fan.dat'

!	HR3 LightCone Data
		character(char_len) :: LClabelChar != '0'
		character(char_len) :: halo_data_file != '../../data/input/HR3LC'//trim(adjustl(LClabelChar))//'_57w_600to1787.dat'
!		character(char_len), parameter :: halo_data_file = '../../data/input/'
!		character(char_len), parameter :: halo_data_file = '../../data/input/'
				
	!----------------------------------
	! ( 2 ). Name of the project.	
	!----------------------------------
	
		character(char_len), parameter :: halo_data_name = 'MD_halos_160w_1000-1600-shell'
	
		character(char_len), parameter :: output_dir = '../../data/output/'!//trim(adjustl(halo_data_name))//'/'
		
	!----------------------------------
	! ( 3 ). Real cosmological parameters	
	!----------------------------------	

!		real(dl), parameter :: om_dft = 0.50_dl, h_dft = 0.72_dl, w_dft = -1.0_dl		
		real(dl), parameter :: om_dft = 0.26_dl, h_dft = 0.72_dl, w_dft = -1.0_dl
!		real(dl), parameter :: om_dft = 0.27_dl, h_dft = 0.7_dl, w_dft = -1.0_dl
	
	!----------------------------------
	! ( 4 ). Assumed cosmologies.
	!----------------------------------	
		
		integer, parameter :: num_assumed_cosmo = 2
		
		real(dl), parameter :: assumed_cosmo(3,num_assumed_cosmo) = &
			(/ 0.15_dl,  w_dft, h_dft,  0.50_dl,  w_dft, h_dft /)
						   
		character(char_len), parameter :: ascos_name1 = 'om_0p05'
		character(char_len), parameter :: ascos_name2 = 'om_1p0'
		character(char_len), parameter :: ascos_name3 = 'om_0p1'

		character(char_len), parameter :: assumed_cosmo_name(num_assumed_cosmo) = (/ ascos_name1, ascos_name2 /)
	


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
		logical function has_boundary_effect(x,y,z,max_dist,cb_adjust_ratio)
			!dummy
			real(dl) :: x,y,z,max_dist,cb_adjust_ratio
			!local 
			real(dl) :: r,dr, extra_correct = 15.0
			if(max_dist .ge. 1.0e5) then
!				print *, 'Encountered max_dist ', max_dist
				has_boundary_effect = .true.; return
			endif
			dr = max_dist*cb_adjust_ratio + extra_correct
!			dr = 30.0 ! TESTING			
			r = sqrt(x*x+y*y+z*z)
!			dr = 300.0;
			if (abs(r-gbrmin)<dr.or.abs(r-gbrmax)<dr) then
				has_boundary_effect = .true.; return
			endif
			if(abs(x-gbxmin)<dr.or.abs(y-gbymin)<dr.or.abs(z-gbzmin)<dr &
				.or.abs(x-gbxmax)<dr.or.abs(y-gbymax)<dr.or.abs(z-gbzmax)<dr) then
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
		real(dl), allocatable :: tmp(:)
		character(len=char_len) :: outputfile = '~/tmp.txt'
		
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
		
		open(unit=1987,file=outputfile)
		do i = 1, num_halo
			read(1987,*) halo_info(i)%z_real
		enddo
		close(1987)
		
		do i = 1, num_halo
			halo_info(i)%x = halo_data(1,i)
			halo_info(i)%y = halo_data(2,i)
			halo_info(i)%z = halo_data(3,i)
			halo_info(i)%vx = halo_data(4,i)*gbvratio
			halo_info(i)%vy = halo_data(5,i)*gbvratio
			halo_info(i)%vz = halo_data(6,i)*gbvratio
			halo_info(i)%mass = halo_data(7,i)
			call getSC(halo_info(i)%x,halo_info(i)%y,halo_info(i)%z,halo_info(i)%r,halo_info(i)%theta,halo_info(i)%phi)
			halo_info(i)%vlos = (halo_info(i)%x*halo_info(i)%vx + halo_info(i)%y*halo_info(i)%vy + halo_info(i)%z*halo_info(i)%vz) / halo_info(i)%r 
			! redshifts
!			halo_info(i)%z_real = get_z(halo_info(i)%r)
			if(local_pecu_v) then
				halo_info(i)%vx = (1.0_dl + halo_info(i)%z_real) * halo_info(i)%vx 
				halo_info(i)%vx = (1.0_dl + halo_info(i)%z_real) * halo_info(i)%vx
				halo_info(i)%vx = (1.0_dl + halo_info(i)%z_real) * halo_info(i)%vx
				halo_info(i)%vlos = (1.0_dl + halo_info(i)%z_real) * halo_info(i)%vlos
			endif
			halo_info(i)%z_obs = halo_info(i)%z_real + halo_info(i)%vlos / const_c 
			halo_info(i)%r_RSD = de_get_comovr(halo_info(i)%z_obs)
			if(mod(i,int(num_halo/15)) .eq. 1) then
				print *, 'Check r:', i, de_get_comovr(halo_info(i)%z_obs) - ( halo_info(i)%r + halo_info(i)%vlos / Hz(halo_info(i)%z_real) * gb_h )
			endif
		enddo
		
!		allocate(tmp(num_halo))
!		do i =1, num_halo
!			tmp(i) = halo_info(i)%z_real
!		enddo
!		call output_1d(outputfile,tmp,num_halo)
		
		gbrealxmin = halo_info(1)%x
		gbrealxmax = halo_info(1)%x
		gbrealymin = halo_info(1)%y
		gbrealymax = halo_info(1)%y
		gbrealzmin = halo_info(1)%z
		gbrealzmax = halo_info(1)%z
		
		do i = 2, num_halo
			gbrealxmin = min(gbrealxmin, halo_info(i)%x)
			gbrealxmax = max(gbrealxmax, halo_info(i)%x)
			gbrealymin = min(gbrealymin, halo_info(i)%y)
			gbrealymax = max(gbrealymax, halo_info(i)%y)
			gbrealzmin = min(gbrealzmin, halo_info(i)%z)
			gbrealzmax = max(gbrealzmax, halo_info(i)%z)
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
		
!BEGIN TESTING
		if(.false.) then
!			open(unit=1,file='chisqlike_LookAtRSD/MDR1info.dat')
			open(unit=113,file='chisqlike_LookAtRSD/HR3LCQua0_ComovPecuV_info.dat')
!			open(unit=1,file='chisqlike_LookAtRSD/HR3LCQua0info.dat')
			do i = 1, num_halo
				write(113,'(11(e14.7,1x))') halo_info(i)%x, halo_info(i)%y, halo_info(i)%z, halo_info(i)%r, halo_info(i)%r_RSD-halo_info(i)%r, halo_info(i)%vx, halo_info(i)%vy, halo_info(i)%vz, halo_info(i)%vlos, halo_info(i)%z_real, halo_info(i)%z_obs-halo_info(i)%z_real
			enddo
			close(113)
!			stop
		endif
!END TESTING		
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
	



