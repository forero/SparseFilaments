

!####################################
!This module does smooth
!####################################
module ap_grids
use ap_mu_tools
use ap_settings_init
use ap_smooth

	implicit none
	
! settings of the programe	
	integer, parameter :: CIC_order = 1 
	! 1: standard CIC (8 nearest cells);
	! 2: 27 nearest
	
	logical, parameter :: cell_smooth = .false.
	! whether to smooth the cells?
	
! nuisance settings	
!	character(len=char_len) :: rho_gradient_file_name
	
	type :: pixel
!		real(dl) :: x, y, z, r
		real(dl) :: mass = 0.0_dl!, rho = 0.0_dl ! if cs%use_num_density = .true., then mass is number, rho is number density
		logical :: bdeffect = .false., nbbdeffect = .false.
		logical :: has_mass = .false.
!		logical :: has_diff(30) = .false.
!		real(dl) :: diff(3,30)
	end type

	type(pixel), allocatable :: gbpixels(:,:,:)

	integer, parameter :: maxmubins = 50
	type :: mustatype !maximal 50 bins
		real(dl) :: weimuarray(maxmubins)
		integer :: nummuarray(maxmubins)
	end type

	real(dl) :: gbpixelmassmean
! 
!	integer :: gb_num_xyz_mass
!	real(dl), allocatable :: gb_xyz_list(:,:), gb_mass_list(:), gb_bf_mass_list(:), gb_r_list(:)
!	real(dl) :: gbxmin, gbxmax, gbymin, gbymax, gbzmin, gbzmax
!	real(dl) :: gbrmin, gbrmax, gbtotvol
!	real(dl) :: unit_len
!	integer :: gb_n_cellx, gb_n_celly, gb_n_cellz
!	real(dl) :: gbdeltax, gbdeltay, gbdeltaz
!	real(dl) :: gb_cell_vol
	
contains


  !------------------------------------------
  ! clean up allocatable arrays
  !------------------------------------------
	subroutine grid_clean_up()
		if(allocated(gbpixels)) then
			deallocate(gbpixels)
		endif
	end subroutine grid_clean_up
	
  !------------------------------------------
  ! initialize the xyz_mass array
  !------------------------------------------
	subroutine init_pixels(RSD, AP, rl_num_in_x, print_info, use_num_density, do_xyzmassinit)
		! DUMMY
		integer :: RSD, AP
		real(dl) :: rl_num_in_x
		logical :: print_info, use_num_density, do_xyzmassinit, test_print
		! Local
		integer :: i, ix,iy,iz, ix1,ix2, iy1,iy2, iz1,iz2, ipx,ipy,ipz
		real(dl) :: x,y,z,r, x1,y1,z1, x2,y2,z2, rsq,rsq1,rsq2, rminsq,rmaxsq, px,py,pz, mass, ratio, time1, time2, time3
		! Used for testing
		real(dl) :: sumratio
		integer :: icount, numnobd, numnonbbd, numzeromass

		if(print_info) then
		        write(*,'(A,i4,i4,A)'), '  Initializing grids of pixels with RSD, AP = ', RSD, AP, '  ...'
		endif
		
		call do_cell_init(RSD, AP, rl_num_in_x, print_info = .false., do_not_init_cellmat = .true.)

		! thresholds of x/y/z: pixels lies out of this range has no mass
		rminsq = gbrmin**2.0
		rmaxsq = gbrmax**2.0

		! Check boundary effect
		call grid_clean_up()
		allocate(gbpixels(gb_n_cellx, gb_n_celly, gb_n_cellz))
		do ix=1, gb_n_cellx
		do iy=1, gb_n_celly
		do iz=1, gb_n_cellz
			x2 = gbgridxmin + gbdeltax*ix 
			y2 = gbgridymin + gbdeltay*iy 
			z2 = gbgridzmin + gbdeltaz*iz
			if(x2>gbxmax .or. y2>gbymax .or. z2>gbzmax) then
				gbpixels(ix,iy,iz)%bdeffect = .true.
				cycle
			endif
			x1 = x2-gbdeltax
			y1 = y2-gbdeltay
			z1 = z2-gbdeltaz
			rsq1 = x1*x1 + y1*y1 + z1*z1
			rsq2 = x2*x2 + y2*y2 + z2*z2
			if(rsq1 .le. rminsq .or. rsq2 .ge. rmaxsq) then
				gbpixels(ix,iy,iz)%bdeffect = .true.
			endif
		enddo
		enddo
		enddo
		
		! Check neighbour-boundary-effect
		! By default nbbdeffect set as .true. if bdeffect is .true.
		do ix=1, gb_n_cellx
		do iy=1, gb_n_celly
		do iz=1, gb_n_cellz
			if(gbpixels(ix,iy,iz)%bdeffect) then
				gbpixels(ix,iy,iz)%nbbdeffect = .true.
				cycle
			endif			
			ix1=max(ix-1,1);iy1=max(iy-1,1);iz1=max(iz-1,1)
			ix2=min(ix+1,gb_n_cellx);iy2=min(iy+1,gb_n_celly);iz2=min(iz+1,gb_n_cellz);
			if(gbpixels(ix1,iy1,iz1)%bdeffect.or.gbpixels(ix2,iy2,iz2)%bdeffect) then
				gbpixels(ix,iy,iz)%nbbdeffect = .true.
			endif
		enddo
		enddo
		enddo
		
		if(use_num_density) &
			gb_mass_list = 1.0_dl
		
		! Calculating CIC mass for all pixels.
		call cpu_time(time2)
		do i = 1, gb_num_xyz_mass
			x= gb_xyz_list(1,i); y= gb_xyz_list(2,i); z= gb_xyz_list(3,i)
			r=gb_r_list(i)
			mass = gb_mass_list(i)

			ix=max(min(int((x-gbgridxmin)/gbdeltax)+1,gb_n_cellx),1)
			iy=max(min(int((y-gbgridymin)/gbdeltay)+1,gb_n_celly),1)
			iz=max(min(int((z-gbgridzmin)/gbdeltaz)+1,gb_n_cellz),1)

			ix2=int((x-gbgridxmin-gbdeltax/2.0)/gbdeltax+2)
			ix1=ix2-1
			iy2=int((y-gbgridymin-gbdeltay/2.0)/gbdeltay+2)
			iy1=iy2-1
			iz2=int((z-gbgridzmin-gbdeltaz/2.0)/gbdeltaz+2)
			iz1=iz2-1
			
			! if any of ix1/iy1/iz1 smaller than 1, assign this halo to this single pixel
!			if(ix1.le.0 .or. iy1.le.0 .or. iz1.le.0 .or. &
!				ix2.gt.gb_n_cellx .or. iy2.gt.gb_n_celly .or. iz2.gt.gb_n_cellz) then
!				gbpixels(ix,iy,iz)%mass = gbpixels(ix,iy,iz)%mass + mass
!				cycle
!			endif
			! if the cell has lowest r touches the boundary, assign this halo to this single pixel
!			if(gbpixels(ix1,iy1,iz1)%bdeffect .or. gbpixels(ix2,iy2,iz2)%bdeffect) then
!				gbpixels(ix,iy,iz)%mass = gbpixels(ix,iy,iz)%mass + mass
!				cycle
!			endif

			! o.w., do CIC 
!			sumratio = 0.0
			do ipx=max(1,ix1),min(ix2,gb_n_cellx)
			do ipy=max(1,iy1),min(iy2,gb_n_celly)
			do ipz=max(1,iz1),min(iz2,gb_n_cellz)
!				if(ix1<1 .or. iy1<1 .or. iz1<1 .or. ix2>gb_n_cellx .or. iy2>gb_n_celly .or. iz2>gb_n_cellz) cycle
				call cell_pos(ipx,ipy,ipz,px,py,pz)
				ratio = (gbdeltax-abs(x-px))/gbdeltax * (gbdeltay-abs(y-py))/gbdeltay * (gbdeltaz-abs(z-pz))/gbdeltaz
				gbpixels(ipx,ipy,ipz)%mass = gbpixels(ipx,ipy,ipz)%mass + mass*ratio
!				sumratio = sumratio+ratio
!				if(test_print) then
!					print *, 'ipx,ipy,ipz=',ipx,ipy,ipz
!					print *, 'px,py,pz=',real(px),real(py),real(pz)
!					print *, 'ratio=',ratio
!				endif
			enddo
			enddo
			enddo
		enddo
		call cpu_time(time3)

		! print info
		if(print_info) then	
			numnobd = 0
			numnonbbd = 0
			numzeromass = 0	
			do ix=1, gb_n_cellx
			do iy=1, gb_n_celly
			do iz=1, gb_n_cellz
				if(gbpixels(ix,iy,iz)%bdeffect.eq..false.) numnobd = numnobd + 1
				if(gbpixels(ix,iy,iz)%nbbdeffect.eq..false.) numnonbbd = numnonbbd + 1
				if(gbpixels(ix,iy,iz)%bdeffect.eq..false. .and. gbpixels(ix,iy,iz)%mass.eq.0.0_dl) &
					numzeromass = numzeromass + 1
			enddo
			enddo
			enddo
			print *, 'tot#/nobd#/nonbbd#/nullmass# = ', gb_n_cellx*gb_n_celly*gb_n_cellz, numnobd, numnonbbd, numzeromass
			print *, 'Time used in init_xyz_mass: ', time2-time1
			print *, 'Time used in init_pixels  : ', time3-time2
		endif
	end subroutine init_pixels
	

  !------------------------------------------
  ! write pixels info to file for investigation
  !------------------------------------------	
  	subroutine tpCor(mindr,maxdr,mu1,mu2,countmu1,countmu2,diffmu,numbins,mustat)
  		!Dummy
  		real(dl) :: mindr,maxdr, diffmu, mu1,mu2
  		integer :: countmu1, countmu2, numbins
  		type(mustatype) :: mustat
  		!Local
  		integer :: i, ix1,iy1,iz1, ix2,iy2,iz2, di, numpixel,numCor
  		real(dl) :: x1,y1,z1, x2,y2,z2, x,y,z, r, dx,dy,dz, dr,drsq,maxdrsq,mindrsq, mu, nmu1,nmu2, massvar, massmean
  		real(dl) :: rlnumbins, weimuarray(maxmubins)
  		integer :: nummuarray(maxmubins)
  		real(dl), allocatable :: massarray(:), rarray(:)
  		! testing variables
  		integer :: icount
  		real(dl) :: time1,time2
  		
		if(numbins > maxmubins) then
			print *, 'ERROR (tpCor)! maximal # of bins is ', maxmubins
			stop
		endif
		rlnumbins = numbins + 0.0_dl
		nummuarray = 0
		weimuarray = 0.0_dl
  		
  		di = int(maxdr/(gbdeltax)+1.0)
  		mindrsq = mindr**2.0
  		maxdrsq = maxdr**2.0
  		numCor = 0
  		countmu1 = 0; mu1=0
  		countmu2 = 0; mu2=0
  		icount = 0
  		
  		call cpu_time(time1)
  		
  		! We use very strict condition: skip boundary pixels (any index is 1 or gb_n_cell),
  		!  also skip any pixels whose neighbour has boundary effect 
!  		do ix1=1,gb_n_cellx
!  		do iy1=1,gb_n_celly
!  		do iz1=1,gb_n_cellz
  		do ix1=2,gb_n_cellx-1
  		do iy1=2,gb_n_celly-1
  		do iz1=2,gb_n_cellz-1
  			if(gbpixels(ix1,iy1,iz1)%nbbdeffect) cycle
  			call cell_pos(ix1,iy1,iz1,x1,y1,z1)
!  			do ix2=ix1, min(ix1+di,gb_n_cellx)
!  			do iy2=max(iy1-di,1), min(iy1+di,gb_n_celly)
!  			do iz2=max(iz1-di,1), min(iz1+di,gb_n_cellz)
  			do ix2=ix1, min(ix1+di,gb_n_cellx-1)
  			do iy2=max(iy1-di,2), min(iy1+di,gb_n_celly-1)
  			do iz2=max(iz1-di,2), min(iz1+di,gb_n_cellz-1)

 				if(ix2 .eq. ix1) then
 					if(iy2.lt.iy1) cycle
  					if(iy2.eq.iy1 .and. iz2.le.iz1) cycle
  				endif
	  			if(gbpixels(ix2,iy2,iz2)%nbbdeffect .or. gbpixels(ix1,iy1,iz1)%mass.eq.0.0_dl&
	  				.or. gbpixels(ix2,iy2,iz2)%mass.eq.0.0_dl) cycle
	  			call cell_pos(ix2,iy2,iz2,x2,y2,z2)
				dx=x2-x1; dy=y2-y1; dz=z2-z1;
				drsq = dx*dx+dy*dy+dz*dz
				if(drsq.gt.maxdrsq .or. drsq.le.mindrsq) cycle
	  			numCor = numCor + 1
				x=(x1+x2)*0.5; y=(y1+y2)*0.5; z=(z1+z2)*0.5;
				dr = sqrt(drsq)
				r = sqrt(x*x+y*y+z*z)
				mu = (x*dx+y*dy+z*dz)/dr/r
!				massvar = 1.0/(gbpixels(ix1,iy1,iz1)%mass*gbpixels(ix2,iy2,iz2)%mass)
!				massvar = 1.0*(gbpixels(ix1,iy1,iz1)%mass*gbpixels(ix2,iy2,iz2)%mass)
!				massvar = sqrt(gbpixels(ix1,iy1,iz1)%mass*gbpixels(ix2,iy2,iz2)%mass)
!				massvar = sqrt(gbpixels(ix1,iy1,iz1)%mass*gbpixels(ix2,iy2,iz2)%mass)
!				massvar = gbpixels(ix1,iy1,iz1)%mass + gbpixels(ix2,iy2,iz2)%mass
!				massvar = 1.0/(gbpixels(ix1,iy1,iz1)%mass*gbpixels(ix2,iy2,iz2)%mass)!*massmean**2.0
				massvar = 1.0/sqrt(gbpixels(ix1,iy1,iz1)%mass*gbpixels(ix2,iy2,iz2)%mass)!*massmean**2.0
				if(abs(mu).le.0.5) then
					countmu1 = countmu1 + 1
					mu1 = mu1+massvar
				else
					countmu2 = countmu2 + 1
					mu2 = mu2+massvar
				endif
				! 
				i = int(abs(mu) * rlnumbins) + 1
				i = min(i,numbins)
				weimuarray(i) = weimuarray(i) + massvar
				nummuarray(i) = nummuarray(i) + 1
  			enddo
  			enddo
  			enddo
		enddo
		enddo
		enddo
		
		mustat%weimuarray = weimuarray
		mustat%nummuarray = nummuarray
		
		call cpu_time(time2)

!		chisq = ((mu1 - (mu1+mu2)/2.0) / (mu1*sqrt(countmu1+0.0)/(countmu1+0.0))) ** 2.0 + &
!				((mu2 - (mu1+mu2)/2.0) / (mu2*sqrt(countmu2+0.0)/(countmu2+0.0))) ** 2.0
		nmu1 = mu1/dble(countmu1); nmu2 = mu2/dble(countmu2)
		diffmu = abs(nmu2-nmu1)/((nmu2+nmu1)*0.5)
!		print *, 'Time used in Cor: ', time2-time1
!		print *, '# of Cor = ', numCor
!		print *, '# of mu1/mu2 = ', countmu1, countmu2
!		print *, 'mu1, mu2 = ', mu1, mu2
!		print *, 'normalized mu1/mu2 = ', real(nmu1), real(nmu2)
!		print *, 'ratio of diff: ', diffmu
	end subroutine tpCor

  !------------------------------------------
  ! diff mu for 2p Cor
  !------------------------------------------	
	subroutine tpCormu(omegam, w, h, rl_num_in_x, refDx, drmin, drmax, diffmu, RSDdiffmu, calc_cos, numbins, mustat, RSDmustat)
		!Dummy
		real(dl), intent(in) :: omegam, w, h, rl_num_in_x, refDx, drmin, drmax
		real(dl), intent(out) :: diffmu, RSDdiffmu
		logical :: calc_cos
		integer :: numbins
		type(mustatype) :: mustat, RSDmustat
		!Local
		integer :: i, AP, countmu1, countmu2
		real(dl) :: deltaratio, mu1, mu2
		logical :: avgr = .true.
		AP = 1
		gb_omegam 	= omegam
		gb_w 		= w
		gb_h 		= h
		if(calc_cos) then
			call de_calc_comovr()
			do i = 1, num_halo
				halo_info(i)%r_AP(AP) = de_get_comovr(dble(halo_info(i)%z_real))
				halo_info(i)%r_AP_RSD(AP) = de_get_comovr(dble(halo_info(i)%z_obs))
			enddo
		endif
		
		!No RSD consraint
		call init_pixels(RSD=0, AP=1, rl_num_in_x=rl_num_in_x, &
			print_info=.false., use_num_density=.true., do_xyzmassinit=.true.)
		if(avgr) then
			call do_grid_avg_mlist()
			call init_pixels(RSD=0, AP=1, rl_num_in_x=rl_num_in_x, &
				print_info=.false., use_num_density=.true., do_xyzmassinit=.true.)
		endif
		deltaratio = (gbxmax - gbxmin) / refDx
!		deltaratio = (gbgridxmax - gbgridxmin) / refDx
		call tpCor(drmin*deltaratio, drmax*deltaratio, mu1, mu2, countmu1, countmu2, diffmu, numbins, mustat)
		
		!With RSD consraint
		call init_pixels(RSD=1, AP=1, rl_num_in_x=rl_num_in_x, &
			print_info=.false., use_num_density=.true., do_xyzmassinit=.true.)
		if(avgr) then
			call do_grid_avg_mlist()
			call init_pixels(RSD=1, AP=1, rl_num_in_x=rl_num_in_x, &
				print_info=.false., use_num_density=.true., do_xyzmassinit=.true.)
		endif
		deltaratio = (gbxmax - gbxmin) / refDx
!		deltaratio = (gbgridxmax - gbgridxmin) / refDx
		call tpCor(drmin*deltaratio, drmax*deltaratio, mu1, mu2, countmu1, countmu2, RSDdiffmu, numbins, RSDmustat)
	end subroutine tpCormu

  !------------------------------------------
  ! multiple diff mu for 2p Cor
  !------------------------------------------
	subroutine tpCormlmu(omegam, w, h, cubemin, cubemax, numcube, refDx, drmin, drmax, diffmu, RSDdiffmu, &
		numbins, mustat, RSDmustat)
		!Dummy
		real(dl), intent(in) :: omegam, w, h, cubemin, cubemax, refDx, drmin, drmax
		integer, intent(in) :: numcube, numbins
		real(dl), intent(out) :: diffmu(numcube), RSDdiffmu(numcube)
		type(mustatype), intent(out) :: mustat(numcube), RSDmustat(numcube)
		!Local
		integer :: i
		real(dl) :: dcube
		
		dcube = (cubemax-cubemin) / (numcube-1.0)
		
		i=1
		call tpCormu(omegam, w, h, cubemin+(i-1)*dcube, refDx, drmin, drmax, diffmu(i), RSDdiffmu(i), .true., &
			numbins, mustat(i), RSDmustat(i))
		
		do i = 2, numcube
			call tpCormu(omegam, w, h, cubemin+(i-1)*dcube, refDx, drmin, drmax, diffmu(i), RSDdiffmu(i), .false., &
				numbins, mustat(i), RSDmustat(i))
		enddo
	end subroutine tpCormlmu
	
  !------------------------------------------
  ! do averaging along the distance
  !------------------------------------------	
	subroutine do_grid_avg_mlist()
		integer :: numpixel, ix1,iy1,iz1, i
		real(dl) :: x1,y1,z1
		real(dl), allocatable :: rarray(:), massarray(:)

		gb_mass_list = 1.0_dl

  		numpixel = 0
  		do ix1=2,gb_n_cellx-1
  		do iy1=2,gb_n_celly-1
  		do iz1=2,gb_n_cellz-1
  			if(gbpixels(ix1,iy1,iz1)%nbbdeffect) cycle
  			numpixel = numpixel + 1
  		enddo
  		enddo
  		enddo
  		
  		allocate(rarray(numpixel), massarray(numpixel))
  		i = 0
  		do ix1=2,gb_n_cellx-1
  		do iy1=2,gb_n_celly-1
  		do iz1=2,gb_n_cellz-1
  			if(gbpixels(ix1,iy1,iz1)%nbbdeffect) cycle
  			i = i + 1
  			call cell_pos(ix1,iy1,iz1,x1,y1,z1)
  			rarray(i) = sqrt(x1*x1+y1*y1+z1*z1)
  			massarray(i) = gbpixels(ix1,iy1,iz1)%mass
  		enddo
  		enddo
  		enddo
	
		call grid_avg_mlist(rarray, massarray, numpixel, .false., 10, 1)
	end subroutine do_grid_avg_mlist
	
  !------------------------------------------
  ! Reset the mass_list: normalized by avg;
  !------------------------------------------
	subroutine grid_avg_mlist(distance_list, quan_list, n, print_info, nbins_rhoav, nsm)
		!Dummy
		integer :: n, nbins_rhoav, nsm
		logical :: print_info
		real(dl) :: distance_list(n), quan_list(n), remov_dist_ratio
		!Local
		real(dl) :: boundary_rmin, boundary_rmax, delta_r, r, quan_av
		real(dl), allocatable :: quan_av_list(:), quan_er_list(:), binned_r_list(:)
		integer :: i_sm, i, i1, i2, i3, binned_i
		logical :: use_intpl_rho = .true.
	
!		call find_min_max(r_list,size(r_list),rmin,rmax)
		boundary_rmin = gbrmin
		boundary_rmax = gbrmax
		delta_r = (boundary_rmax - boundary_rmin) / (nbins_rhoav + 0.0)
  		
  		allocate(quan_av_list(nbins_rhoav), quan_er_list(nbins_rhoav), binned_r_list(nbins_rhoav))
  		do i_sm = 1, nsm
  			call binned_quan(quan_list, distance_list, boundary_rmin, boundary_rmax, size(quan_list), nbins_rhoav, &
		  		quan_av_list, quan_er_list, binned_r_list)
	  		do i = 1, size(gb_r_list)
  				r = gb_r_list(i)
  				binned_i = max(min(int((r-boundary_rmin)/delta_r)+1,nbins_rhoav),1)
  				if(use_intpl_rho) then
  					call ilist(binned_i, 1, 1, nbins_rhoav,i1,i3)
  					i2 = i1 + 1
  					quan_av = adaptp_intpl_vl(r,binned_r_list(i1),quan_av_list(i1),binned_r_list(i2),quan_av_list(i2),binned_r_list(i3),quan_av_list(i3))
  					gb_mass_list(i) = gb_mass_list(i) / quan_av
  				else
  					gb_mass_list(i) = gb_mass_list(i) / quan_av_list(binned_i)
  				endif
  			enddo
  		enddo
  	end subroutine grid_avg_mlist

!	Calcuate the chisq based on pixels
  !------------------------------------------
  ! diff mu for 2p Cor
  !------------------------------------------	
	subroutine CICgfchisq(omegam, w, h, cubemin, cubemax, numcube, calc_cos, &
		droparray, ndrop, NORSDchisq, RSDchisq)
		!Dummy
		real(dl), intent(in) :: omegam, w, h, cubemin, cubemax, droparray(ndrop)
		logical, intent(in) :: calc_cos
		integer, intent(in) :: numcube, ndrop
		real(dl), intent(out) :: NORSDchisq(ndrop), RSDchisq(ndrop)
		!Local
		real(dl) :: NORSDchisqarray(numcube), RSDchisqarray(numcube)
		real(dl), allocatable :: mudata(:), tmp(:), massdata(:),  sortedmudata(:)
		integer, allocatable :: ordermu(:)
		logical :: avgr = .true.
		integer :: AP, i, j, ix,iy,iz, n, maxnummu, m
		real(dl) :: x,y,z, dx,dy,dz, rl_num_in_x,dcube
		integer :: di = 2
		AP = 1
		gb_omegam 	= omegam
		gb_w 		= w
		gb_h 		= h
		if(calc_cos) then
			call de_calc_comovr()
			do i = 1, num_halo
				halo_info(i)%r_AP(AP) = de_get_comovr(dble(halo_info(i)%z_real))
				halo_info(i)%r_AP_RSD(AP) = de_get_comovr(dble(halo_info(i)%z_obs))
			enddo
		endif
		
		if(numcube .eq. 1) then
			dcube = 0.0_dl
		else
			dcube = (cubemax-cubemin)/(numcube-1.0)
		endif
!		print *, cubemin, cubemax, dcube
		
		NORSDchisq = 0.0
		!No RSD consraint
		do i = 1, numcube
			rl_num_in_x = cubemin + (i-1)*dcube
!			print *, 'rl_num_in_x: ', rl_num_in_x
			call init_pixels(RSD=0, AP=1, rl_num_in_x=rl_num_in_x, &
				print_info=.false., use_num_density=.true., do_xyzmassinit=.true.)
			if(avgr) then
				call do_grid_avg_mlist()
				call init_pixels(RSD=0, AP=1, rl_num_in_x=rl_num_in_x, &
					print_info=.false., use_num_density=.true., do_xyzmassinit=.true.)
			endif
			maxnummu = (gb_n_cellx-2)*(gb_n_celly-2)*(gb_n_cellz-2)
			allocate(mudata(maxnummu), massdata(maxnummu))
			n = 0
			do ix = 2+di, gb_n_cellx-1-di
			do iy = 2+di, gb_n_celly-1-di
			do iz = 2+di, gb_n_cellz-1-di
				if(gbpixels(ix-di,iy,iz)%nbbdeffect .or. &
	 				gbpixels(ix,iy-di,iz)%nbbdeffect  .or. &
					gbpixels(ix,iy,iz-di)%nbbdeffect  .or. &
					gbpixels(ix+di,iy,iz)%nbbdeffect  .or. &
					gbpixels(ix,iy+di,iz)%nbbdeffect  .or. &
					gbpixels(ix,iy,iz+di)%nbbdeffect ) cycle
				call cell_pos(ix,iy,iz,x,y,z)
				dx=gbpixels(ix+di,iy,iz)%mass - gbpixels(ix-di,iy,iz)%mass
				dy=gbpixels(ix,iy+di,iz)%mass - gbpixels(ix,iy-di,iz)%mass
				dz=gbpixels(ix,iy,iz+di)%mass - gbpixels(ix,iy,iz-di)%mass
				n = n+1
				mudata(n) = get_mu_of_gradient(x,y,z,dx,dy,dz)				
				massdata(n) = gbpixels(ix,iy,iz)%mass * sqrt(dx*dx+dy*dy+dz*dz)
			enddo
			enddo
			enddo
			allocate(sortedmudata(n),ordermu(n))
			do j = 1, n
				ordermu(j) = j
			enddo
!			print *, 'A'
			call QSort2(massdata(1:n),ordermu,n)
!			print *, 'B'
			do j = 1, n
				sortedmudata(j) = mudata(ordermu(j))
			enddo
			do j = 1, ndrop
				m = max(min(int(n*(1-droparray(j))+0.5), n),10)
!				print *, 'j,droarray(j),m,n= ',j,real(droparray(j)),m,n
				NORSDchisq(j) = NORSDchisq(j) + chisq_of_mu_data3(sortedmudata(1:m),m)
!				print *, NORSDchisq(j)
			enddo
			deallocate(mudata, massdata, sortedmudata, ordermu)
		enddo
		
		!With RSD consraint
		RSDchisq = 0.0
		do i = 1, numcube
			!No RSD consraint
			rl_num_in_x = cubemin + (i-1)*dcube
			call init_pixels(RSD=1, AP=1, rl_num_in_x=rl_num_in_x, &
				print_info=.false., use_num_density=.true., do_xyzmassinit=.true.)
			if(avgr) then
				call do_grid_avg_mlist()
				call init_pixels(RSD=1, AP=1, rl_num_in_x=rl_num_in_x, &
					print_info=.false., use_num_density=.true., do_xyzmassinit=.true.)
			endif
			maxnummu = (gb_n_cellx-2)*(gb_n_celly-2)*(gb_n_cellz-2)
			allocate(mudata(maxnummu), massdata(maxnummu))
			n = 0
			do ix = 2+di, gb_n_cellx-1-di
			do iy = 2+di, gb_n_celly-1-di
			do iz = 2+di, gb_n_cellz-1-di
				if(gbpixels(ix-di,iy,iz)%nbbdeffect .or. &
	 				gbpixels(ix,iy-di,iz)%nbbdeffect  .or. &
					gbpixels(ix,iy,iz-di)%nbbdeffect  .or. &
					gbpixels(ix+di,iy,iz)%nbbdeffect  .or. &
					gbpixels(ix,iy+di,iz)%nbbdeffect  .or. &
					gbpixels(ix,iy,iz+di)%nbbdeffect ) cycle
				if(gbpixels(ix+di,iy,iz)%mass.le.1.0d-50 .or. &
					gbpixels(ix,iy+di,iz)%mass.le.1.0d-50 .or.&
					gbpixels(ix,iy,iz+di)%mass.le.1.0d-50 .or.&
					gbpixels(ix-di,iy,iz)%mass.le.1.0d-50 .or.&
					gbpixels(ix,iy-di,iz)%mass.le.1.0d-50 .or.&
					gbpixels(ix,iy,iz-di)%mass.le.1.0d-50) cycle
				call cell_pos(ix,iy,iz,x,y,z)					
				dx=gbpixels(ix+di,iy,iz)%mass - gbpixels(ix-di,iy,iz)%mass
				dy=gbpixels(ix,iy+di,iz)%mass - gbpixels(ix,iy-di,iz)%mass
				dz=gbpixels(ix,iy,iz+di)%mass - gbpixels(ix,iy,iz-di)%mass
				n = n+1
				mudata(n) = get_mu_of_gradient(x,y,z,dx,dy,dz)
				massdata(n) = gbpixels(ix,iy,iz)%mass * sqrt(dx*dx+dy*dy+dz*dz)
			enddo
			enddo
			enddo
			allocate(sortedmudata(n),ordermu(n))
			do j = 1, n
				ordermu(j) = j
			enddo
			call QSort2(massdata(1:n),ordermu,n)
			do j = 1, n
				sortedmudata(j) = mudata(ordermu(j))
			enddo
			do j = 1, ndrop
				m = max(min(int(n*(1-droparray(j))+0.5), n),10)
				RSDchisq(j) = RSDchisq(j) + chisq_of_mu_data3(sortedmudata(1:m), m)
			enddo
			deallocate(mudata, massdata, sortedmudata, ordermu)
		enddo
!		deltaratio = (gbgridxmax - gbgridxmin) / refDx
!		call tpCor(drmin*deltaratio, drmax*deltaratio, mu1, mu2, countmu1, countmu2, diffmu, numbins, mustat)
	end subroutine CICgfchisq

end module ap_grids


