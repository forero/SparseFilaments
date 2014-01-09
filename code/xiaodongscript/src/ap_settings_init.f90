
! List of settings 


!####################################
!This module does smooth
!####################################
module ap_settings_init
use ap_cosmo_funs

	implicit none

	real(dl), parameter :: vol_frac = 0.5_dl ! half sphere

	real(dl) :: gbgridxmin, gbgridxmax, gbgridymin, gbgridymax, gbgridzmin, gbgridzmax, gbgridrmin, gbgridrmax
	real(dl) :: gbxmin, gbxmax, gbymin, gbymax, gbzmin, gbzmax, gbrmin, gbrmax
	real(dl) :: gbnoRSDxmin, gbnoRSDxmax, gbnoRSDymin, gbnoRSDymax, gbnoRSDzmin, gbnoRSDzmax, gbnoRSDrmin, gbnoRSDrmax
	real(dl) :: gbrealxmin, gbrealxmax, gbrealymin, gbrealymax, gbrealzmin, gbrealzmax, gbrealrmin, gbrealrmax
	real(dl) :: gbvratio = 1.0_dl
	logical  :: local_pecu_v = .true.
	logical  :: dotsbe !TESTING
	character(len=char_len), public :: basestr, basestr_om, basestr_om_cube_RSD

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
		
		integer, parameter :: num_assumed_cosmo = 1
		
		real(dl), parameter :: assumed_cosmo(3,num_assumed_cosmo) = &
			(/ 0.16_dl,  w_dft, h_dft /)
						   
		character(char_len), parameter :: ascos_name1 = 'om_0p05'
		character(char_len), parameter :: ascos_name2 = 'om_1p0'
		character(char_len), parameter :: ascos_name3 = 'om_0p1'

		character(char_len) :: assumed_cosmo_name(num_assumed_cosmo)
	
!--------------------------------------------------------------------
!--------------------------------------------------------------------	
!--  End of settings of the programme   -----------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
	
	type :: point
		real :: mass, x,y,z
		real :: vx,vy,vz,vlos
		real :: vxth,vyth,vzth,vlosth
		logical :: vth_er = .false. ! whether there is error due to boundary touching
		real :: r,theta, phi
		real :: z_real, z_obs,  r_RSD
		real :: r_AP(num_assumed_cosmo), r_AP_RSD(num_assumed_cosmo), r_AP_RSD_VCOR(num_assumed_cosmo)
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
			real(dl) :: r,dr, extra_correct = 8.0!15.0
			if(max_dist .ge. 1.0e5) then
!				print *, 'Encountered max_dist ', max_dist
				has_boundary_effect = .true.; return
			endif
			dr = max_dist*cb_adjust_ratio + extra_correct
			r = sqrt(x*x+y*y+z*z)
			if(r<gbrmin.or.r>gbrmax.or.r>gbgridrmax) then
				has_boundary_effect = .true.; return
			endif
			if (abs(r-gbrmin)<dr.or.abs(r-gbrmax)<dr.or.abs(r-gbgridrmin)<dr.or.abs(r-gbgridrmax)<dr) then
				has_boundary_effect = .true.; return
			endif
			if(abs(x-gbxmin)<dr.or.abs(y-gbymin)<dr.or.abs(z-gbzmin)<dr &
				.or.abs(x-gbxmax)<dr.or.abs(y-gbymax)<dr.or.abs(z-gbzmax)<dr &
				.or.abs(x-gbgridxmin)<dr.or.abs(y-gbgridymin)<dr.or.abs(z-gbgridzmin)<dr &
				.or.abs(x-gbgridxmax)<dr.or.abs(y-gbgridymax)<dr.or.abs(z-gbgridzmax)<dr) then
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
			vol_fun = (4.0*const_pi/3.0) * (r2**3.0-r1**3.0) * vol_frac
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
		real(dl) :: shell_rmin = 500.0, shell_rmax = 800.0, rsq, r, rfrac, randomx, randomy
		real(dl) :: rsq1, rsq2, x,y,z, minrsq=0.0!1000.0**2.0
		real(dl), allocatable :: tmp(:,:)
		integer :: i, j, n, num

		write(*,'(A,A)') '   (read_in_halo_data) Read in the halo data from: ', trim(adjustl(halo_data_file))	
		
		! directly read in from ready-data file.
!		if(allocated(halo_data)) &
!			deallocate(halo_data)
!		call read_in_revfmt(halo_data_file, 7, num_halo, halo_data)

		! Do selections:::
		write(*,'(A,f12.5)') '   (read_in_halo_data) Applying distance cut: r > ', sqrt(minrsq)
		call read_in(halo_data_file, 7, n, tmp)
		num = 0
		do i = 1, n
			if(tmp(i,1)**2.0+tmp(i,2)**2.0+tmp(i,3)**2.0 > minrsq) num=num+1
		enddo
		num_halo = num
		if(allocated(halo_data)) deallocate(halo_data)		
		allocate(halo_data(7,num_halo))
		j = 1
!		print *, 'TESTING!!! EMBEDING AN ARTIFICAL DROPPING OF R!!!'
		do i = 1, n
			rsq = tmp(i,1)**2.0+tmp(i,2)**2.0+tmp(i,3)**2.0
			r = sqrt(rsq)
			rfrac = (r-500.0) / 1800.0
			call random_number(randomx)
			call random_number(randomy)
			! rfrac increases with distance; at large distance r will decreases (up to a fraction of 30% at r=1800).
!			if((randomy<0.7.or.(randomy>0.7.and.rfrac<randomx)).and.rsq > minrsq) then
			if(.true.) then
				halo_data(1:7,j) = tmp(i,1:7)
				j = j + 1
			endif
		enddo
		num_halo = j-1
!		if(j.ne.num+1) then
!			print *, 'ERROR (read_in_halo_data)!!! # mismatches: ', j, num_halo+1
!			stop
!		endif
		deallocate(tmp)
	end subroutine read_in_halo_data


  !------------------------------------------
  ! initializing the halo_info array
  !------------------------------------------
	subroutine init_halo_info()
		integer :: i, i_cosmo, n
		real(dl), allocatable :: tmp(:)
		character(len=char_len) :: outputfile = '~/tmp.txt'
		!Testing variables
		real(dl) :: testdiff
		integer, parameter :: binnum = 200
		integer :: nowbinnum
		real(dl) :: binnednum(binnum), binnedminmass(binnum), binnedmaxmass(binnum), &
			binnedminred(binnum), binnedmaxred(binnum)
		character(len=char_len) :: binstatfile = 'Test/halo_binned_stat.txt'
		
		if(.not. allocated(halo_data)) then
			print *, 'ERROR (init_halo_info)! halo_data not allocated. Unable to initialize halo_info!'
			stop
		endif
		
		print *, '  (init_halo_info) Initializing the halo_info...'
		print *, '                   # of halos: ', num_halo
		! initialize the comoving r array
		gb_omegam = om_dft; gb_w = w_dft; gb_h = h_dft
		write(*,'(A,f6.2,A,f6.2,A,f6.2)') '                    Default Cosmology: Omegam/w/h = ', gb_omegam,' /',gb_w,' /',gb_h
		call de_calc_comovr()
		if(allocated(halo_info)) &
			deallocate(halo_info)
		allocate(halo_info(num_halo))
		
!		open(unit=1987,file=outputfile)
!		do i = 1, num_halo
!			read(1987,*) halo_info(i)%z_real
!		enddo
!		close(1987)
		
		testdiff = 0.0_dl
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
			! use some interpolation method to get z (very fast, accurate enough)
			halo_info(i)%z_real = de_zfromintpl(dble(halo_info(i)%r)) 
			if(local_pecu_v) then
				halo_info(i)%vx = (1.0_dl + halo_info(i)%z_real) * halo_info(i)%vx 
				halo_info(i)%vy = (1.0_dl + halo_info(i)%z_real) * halo_info(i)%vy
				halo_info(i)%vz = (1.0_dl + halo_info(i)%z_real) * halo_info(i)%vz
				halo_info(i)%vlos = (1.0_dl + halo_info(i)%z_real) * halo_info(i)%vlos
			endif
			halo_info(i)%z_obs = halo_info(i)%z_real + halo_info(i)%vlos / const_c 
			halo_info(i)%r_RSD = de_get_comovr(dble(halo_info(i)%z_obs))
			if(mod(i,int(num_halo/15)) .eq. 1) then
				print *, '                        Check z/r:', i, real(halo_info(i)%z_real-get_z(dble(halo_info(i)%r))), &
				 real(de_get_comovr(dble(halo_info(i)%z_obs))-(halo_info(i)%r+halo_info(i)%vlos/Hz(dble(halo_info(i)%z_real))*gb_h))
			endif
		enddo
		
		gbrealxmin = halo_info(1)%x
		gbrealxmax = halo_info(1)%x
		gbrealymin = halo_info(1)%y
		gbrealymax = halo_info(1)%y
		gbrealzmin = halo_info(1)%z
		gbrealzmax = halo_info(1)%z
		gbrealrmin = halo_info(1)%r 
		gbrealrmax = halo_info(1)%r
		
		do i = 2, num_halo
			gbrealxmin = min(gbrealxmin, halo_info(i)%x)
			gbrealxmax = max(gbrealxmax, halo_info(i)%x)
			gbrealymin = min(gbrealymin, halo_info(i)%y)
			gbrealymax = max(gbrealymax, halo_info(i)%y)
			gbrealzmin = min(gbrealzmin, halo_info(i)%z)
			gbrealzmax = max(gbrealzmax, halo_info(i)%z)
			gbrealrmin = min(gbrealrmin, halo_info(i)%r)
			gbrealrmax = max(gbrealrmax, halo_info(i)%r)
		enddo			

		if(.false.) then
			binnednum  	=  0.0 
			binnedminmass 	=  logzero 
			binnedmaxmass 	= -logzero
			binnedminred  	=  logzero 
			binnedmaxred  	= -logzero
			do i = 1, num_halo
				nowbinnum = min(int((halo_info(i)%r-gbrealrmin)/((gbrealrmax-gbrealrmin)/dble(binnum))+1),binnum)
				binnednum(nowbinnum) = binnednum(nowbinnum) +1
				binnedminmass(nowbinnum) = min(binnedminmass(nowbinnum),halo_info(i)%mass)
				binnedmaxmass(nowbinnum) = max(binnedmaxmass(nowbinnum),halo_info(i)%mass)
				binnedminred(nowbinnum) = min(binnedminred(nowbinnum),halo_info(i)%z_real)
				binnedmaxred(nowbinnum) = max(binnedmaxred(nowbinnum),halo_info(i)%z_real)
			enddo
			open(unit=1,file=binstatfile)
			write(*,*) 'TESTING! Outputing binned statistics to ', trim(adjustl(binstatfile))
			do i = 1, binnum
				write(1,'(i3,1x,7(e14.7,1x))') i, gbrealrmin + (gbrealrmax-gbrealrmin)/dble(binnum)*dble(i-1), &
					gbrealrmin + (gbrealrmax-gbrealrmin)/dble(binnum)*dble(i), &
					binnednum(i),binnedminmass(i),binnedmaxmass(i),binnedminred(i),binnedmaxred(i)
			enddo
			close(1)
			stop
		endif
		
		write(*,'(A,2f9.3)') '                      xrange: ', real(gbrealxmin), real(gbrealxmax)
		write(*,'(A,2f9.3)') '                      yrange: ', real(gbrealymin), real(gbrealymax)
		write(*,'(A,2f9.3)') '                      zrange: ', real(gbrealzmin), real(gbrealzmax)
		write(*,'(A,2f9.3)') '                      rrange: ', real(gbrealrmin), real(gbrealrmax)
		write(*,'(A,e14.7)') '                      # density of sample = ', dble(num_halo) / vol_fun(gbrealrmin,gbrealrmax)
		print *, '                   Calculating postions of wrong cosmologies:'
		do i_cosmo = 1, num_assumed_cosmo
			gb_omegam 	= assumed_cosmo(1,i_cosmo)
			gb_w 		= assumed_cosmo(2,i_cosmo)
			gb_h 		= assumed_cosmo(3,i_cosmo)
			write(*,'(A,i2,A,f6.2,A,f6.2,A,f6.2)') '                      Cosmology ', i_cosmo, &
				': Omegam/w/h = ', gb_omegam,' /',gb_w,' /',gb_h
			call de_calc_comovr()
			do i = 1, num_halo
				halo_info(i)%r_AP(i_cosmo) = de_get_comovr(dble(halo_info(i)%z_real))
				halo_info(i)%r_AP_RSD(i_cosmo) = de_get_comovr(dble(halo_info(i)%z_obs))
				halo_info(i)%r_AP_RSD_VCOR(i_cosmo) = halo_info(i)%r_AP_RSD(i_cosmo)
			enddo
		enddo
		
		call system('mkdir -p '//trim(adjustl(output_dir)))
	end subroutine init_halo_info

  !------------------------------------------
  ! self-introduction of a point
  !------------------------------------------	
	subroutine point_speak(p)
		integer :: i
		type(point) :: p
		real(dl) :: ratio
		print *, '  (point_speak) Information of the point:'
		print *, '                x,y,z = ', p%x, p%y, p%z
		print *, '                vx,vy,vz = ', p%vx, p%vy, p%vz
		print *, '                r,theta,phi = ', p%r, p%theta, p%phi
		print *, '                los speed =', p%vlos
		print *, '                cosmo redshift = ', p%z_real
		print *, '                apparent redshift = ', p%z_obs
		print *, '                apparent comov r  = ', p%r_RSD
		ratio = p%r_RSD/p%r
		write(*,'(A,3f12.5)') '                 apparent xyz: ', p%x*ratio, p%y*ratio, p%z*ratio
		do i = 1, num_assumed_cosmo
			print *, '                Info for cosmology ', i, ':'
			print *, '                  Inferred r (no   RSD) = ', p%r_AP(i)  
			print *, '                  Inferred r (with RSD) = ', p%r_AP_RSD(i)
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

  !------------------------------------------
  ! use given parameters to init AP cosmology
  !------------------------------------------		
	subroutine init_AP_cosmo(omegam,w,h,AP,printinfo)
		! Dummy
		real(dl), intent(in) :: omegam, w, h
		integer, intent(in) :: AP
		logical, intent(in) :: printinfo
		! Local
		integer :: i
		gb_omegam 	= omegam
		gb_w 		= w
		gb_h 		= h
		if(printinfo) then
			write(*,'(3x,A,i3,A,f7.3,f7.3)') '(init_AP_cosmo) Initializing ', AP, 'th AP cosmology omegam/w = ', &
				real(omegam), real(w)
		endif
		call de_calc_comovr()
		do i = 1, num_halo
			halo_info(i)%r_AP(AP) = de_get_comovr(dble(halo_info(i)%z_real))
			halo_info(i)%r_AP_RSD(AP) = de_get_comovr(dble(halo_info(i)%z_obs))
			halo_info(i)%r_AP_RSD_VCOR(AP) = halo_info(i)%r_AP_RSD(AP)
			halo_info(i)%vth_er = .false.
		enddo
	end subroutine init_AP_cosmo		

end module ap_settings_init
	



