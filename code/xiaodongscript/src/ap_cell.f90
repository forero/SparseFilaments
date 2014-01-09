

!####################################
!This module does smooth
!####################################
module ap_cell
use ap_settings_init

	implicit none

!!! IMPORTANT SETTINGS
	logical :: use_num_density = .true.
	logical :: normed_by_density = .true.

!! IMPORTANT SETTINGS
	! Settings of rhocri estimation
	! Separte the data into <nbinrhocri> equal-volume bins, 
	!   get fractional rhocri at each bin, 
	!   then use <nfitrhocri> order polynomial to fit them 
	!    and estimate rhocri at every distance
	integer, parameter :: nbinrhocri = 6, nfitrhocri = 2

!!! IMPORTANT SETTINGS
! settings of ranges for grids.
	integer, public, parameter :: use_fixgridrange = 1 ! use the range from some given values
	real(dl), public :: fixgridxmin, fixgridxmax, fixgridymin, fixgridymax, fixgridzmin, fixgridzmax
	integer, public, parameter :: use_realrange = 2 ! use the range of real data (in right cosmology, without RSD)
	integer, public, parameter :: use_noRSDrange = 3 ! use the range of current cosmology, no RSD case
	integer, public, parameter :: use_currentrange = 4 ! use the range of current cosmology, current RSD case
	integer, public :: gridrangeset = use_noRSDrange 

	
! nuisance settings	
!	character(len=char_len) :: rho_gradient_file_name
	integer :: gb_num_xyz_mass
	real(dl), allocatable :: gb_xyz_list(:,:), gb_mass_list(:), gb_r_list(:), gb_vel_list(:,:) 
!	logical, allocatable :: gb_vcorect_er_list(:)
	real(dl) :: gbtotvol
	real(dl) :: unit_len
	integer :: gb_n_cellx, gb_n_celly, gb_n_cellz
	real(dl) :: gbdeltax, gbdeltay, gbdeltaz
	real(dl) :: gb_cell_vol

	type :: halo_list
		integer :: halo_num = 0
		integer, allocatable :: list(:)
		real :: rhodrhos(4) = -1.0
		real :: saddlerho 
		! rho at the nearby saddle point.
		! Coordinate of the nearby saddle point for cell labeled (ix,iy,iz) &
		!   will be (gbgridxmin + ix*gbdeltax, gbgridymin + iy*gbdeltay, gbgridzmin + iz*gbdeltaz)
		logical :: rhodrho_has_er 
		! Flag indicates that rho / drho has error. M
		!  May be due to three reasons. 
		!   (1) boundary effect;
		!   (2) # of in-sphere halo is too small; 
		!   (3) there is error in vel correction (only happens in with RSD case)
	end type
	
	! array containing which halo in which cell
	type(halo_list), allocatable :: gb_cell_mat(:,:,:)
	
contains


  !------------------------------------------
  ! clean up allocatable arrays
  !------------------------------------------
	subroutine smooth_clean_up()
		if(allocated(gb_xyz_list)) &
			deallocate(gb_xyz_list)
		if(allocated(gb_mass_list)) &
			deallocate(gb_mass_list)
		if(allocated(gb_cell_mat)) &
			deallocate(gb_cell_mat)
		if(allocated(gb_r_list)) &
			deallocate(gb_r_list)	
		if(allocated(gb_vel_list)) &
			deallocate(gb_vel_list)	
!		if(allocated(gb_vcorect_er_list)) &
!			deallocate(gb_vcorect_er_list)
	end subroutine

  !------------------------------------------
  ! initialize the xyz_mass array
  !------------------------------------------
	subroutine init_mult_lists(RSD, AP, print_info)
		integer :: RSD, AP, i
		real(dl) :: x, y, z, ratio, noRSDratio, r
		logical, intent(in) :: print_info

		if(print_info) write(*,'(A,i4,i4)'), '   (init_mult_list) Init lists of xyz/r/mass with RSD, AP = ', RSD, AP
		call smooth_clean_up()
		
		allocate(gb_xyz_list(3,num_halo), gb_mass_list(num_halo), gb_r_list(num_halo), &
			 gb_vel_list(3,num_halo))
		gb_num_xyz_mass = num_halo

		gbxmin = logzero; gbxmax = -logzero;
		gbymin = logzero; gbymax = -logzero;
		gbzmin = logzero; gbzmax = -logzero;
		gbrmin = logzero; gbrmax = -logzero;
		gbnoRSDxmin = logzero;	gbnoRSDxmax = -logzero
		gbnoRSDymin = logzero;	gbnoRSDymax = -logzero
		gbnoRSDzmin = logzero;	gbnoRSDzmax = -logzero
		gbnoRSDrmin = logzero;	gbnoRSDrmax = -logzero
		do i = 1, num_halo
			x = halo_info(i)%x
			y = halo_info(i)%y
			z = halo_info(i)%z
			r = halo_info(i)%r
			if    (RSD .eq. 0 .and. AP .eq. 0) then
				gb_r_list(i) = r 
			elseif(RSD .eq. 1 .and. AP .eq. 0) then
				gb_r_list(i) = halo_info(i)%r_RSD
			elseif(RSD .eq. 0 .and. AP .gt. 0) then
				gb_r_list(i) = halo_info(i)%r_AP(AP)
			elseif(RSD .eq. 1 .and. AP .gt. 0) then
				gb_r_list(i) = halo_info(i)%r_AP_RSD_VCOR(AP)
			else
				print *, 'ERROR! No such RSD, AP : ', RSD, AP
				stop
			endif
			ratio = gb_r_list(i) / r 
			gb_xyz_list(1,i) = x * ratio
			gb_xyz_list(2,i) = y * ratio
			gb_xyz_list(3,i) = z * ratio
			gbxmin = min(gbxmin,gb_xyz_list(1,i)); gbxmax = max(gbxmax,gb_xyz_list(1,i))
			gbymin = min(gbymin,gb_xyz_list(2,i)); gbymax = max(gbymax,gb_xyz_list(2,i))
			gbzmin = min(gbzmin,gb_xyz_list(3,i)); gbzmax = max(gbzmax,gb_xyz_list(3,i))
			gbrmin = min(gbrmin,gb_r_list(i)); gbrmax = max(gbrmax,gb_r_list(i));
			if(RSD.eq.1) then
				noRSDratio = halo_info(i)%r_AP(AP) / r
				gbnoRSDxmin = min(gbnoRSDxmin,x*noRSDratio)
				gbnoRSDxmax = max(gbnoRSDxmax,x*noRSDratio)
				gbnoRSDymin = min(gbnoRSDymin,y*noRSDratio)
				gbnoRSDymax = max(gbnoRSDymax,y*noRSDratio)
				gbnoRSDzmin = min(gbnoRSDzmin,z*noRSDratio)
				gbnoRSDzmax = max(gbnoRSDzmax,z*noRSDratio)
				gbnoRSDrmin = min(gbnoRSDrmin,halo_info(i)%r_AP(AP))
				gbnoRSDrmax = max(gbnoRSDrmax,halo_info(i)%r_AP(AP))
			endif
			if(use_num_density) then
				gb_mass_list(i) = 1.0
			else
				gb_mass_list(i) = halo_info(i)%mass
			endif
			gb_vel_list(1,i) = halo_info(i)%vx
			gb_vel_list(2,i) = halo_info(i)%vy
			gb_vel_list(3,i) = halo_info(i)%vz
		enddo
		!min/max of r

		if(RSD.eq.0) then
			gbnoRSDxmin=gbxmin; gbnoRSDxmax=gbxmax; 
			gbnoRSDymin=gbymin; gbnoRSDymax=gbymax; 
			gbnoRSDzmin=gbzmin; gbnoRSDzmax=gbzmax; 
			gbnoRSDrmin=gbrmin; gbnoRSDrmax=gbrmax; 
		endif
		
		if(normed_by_density) then
			if(print_info) write(*,'(29x,A)'), 'Normalizing mass according to <rho>(r)...'
			! This converts M_i into M_i / <rho>(r_i)
			call normn_mlist(print_info,nbinrhocri,1)
		endif
		
		! lists of errors in v correction (true, has error or false, no error)
!		allocate(gb_vcorect_er_list(num_halo))
!		gb_vcorect_er_list = .false.
		
		if(print_info) then
			write(*,'(2x,A,f12.5,2x,f12.5)')  ' (init_mult_list) xmin / xmax = ', gbxmin, gbxmax
			write(*,'(2x,A,f12.5,2x,f12.5)')  '                  ymin / ymax = ', gbymin, gbymax
			write(*,'(2x,A,f12.5,2x,f12.5)')  '                  zmin / zmax = ', gbzmin, gbzmax
			write(*,'(2x,A,f12.5,2x,f12.5)')  '                  rmin / rmax = ', gbrmin, gbrmax
		endif
	end subroutine init_mult_lists


  !------------------------------------------
  ! get the list of x,y,z at centers of cells
  !------------------------------------------
	subroutine get_cell_pos_list(pos_list)
		real(dl), allocatable :: pos_list(:,:)
		real(dl) :: x,y,z
		integer :: ix, iy, iz, i
		allocate(pos_list(3,gb_n_cellx*gb_n_celly*gb_n_cellz))
		i = 1
		do ix = 1, gb_n_cellx
		do iy = 1, gb_n_celly
		do iz = 1, gb_n_cellz
			call cell_pos(ix,iy,iz,x,y,z)
			pos_list(1:3,i) = (/x,y,z/)
			i = i + 1
		enddo
		enddo
		enddo
	end subroutine get_cell_pos_list
	
  !------------------------------------------
  ! initialize the Calculation
  !------------------------------------------
	subroutine do_cell_init(RSD, AP, rl_num_in_x, print_info, do_not_init_cellmat)
		! Dummy
		integer, intent(in) :: RSD, AP
		real(dl), intent(in) :: rl_num_in_x
		logical, optional, intent(in) :: print_info, do_not_init_cellmat
		! Local
		integer :: numset, i, j, k, ix, iy, iz, num_in_x
		integer(2), allocatable :: mark_num_halo(:,:,:)
		real(dl) :: x,y,z
		integer :: nownum
		integer, parameter :: max_incell_halo_num = 20000  !maximal length of array to save index of halos in one cell
		integer, parameter :: re_alo_num = 5
		real(dl) :: numpixel, numhashalo, maxhalonum, avghalonum
		
		if(print_info) print *, '  (do_cell_init begin) initializing Grid of Cells.'
		if(gridrangeset .eq. use_fixgridrange) then
			if(print_info) print *, '  (do_cell_init) Use fixed range of grids:'
			gbgridxmin = fixgridxmin; gbgridxmax = fixgridxmax
			gbgridymin = fixgridymin; gbgridymax = fixgridymax
			gbgridzmin = fixgridzmin; gbgridzmax = fixgridzmax
		elseif(gridrangeset .eq. use_realrange) then
			if(print_info) print *, '  (do_cell_init) Set range of grids from real data:'
			gbgridxmin = gbrealxmin; gbgridxmax = gbrealxmax;
			gbgridymin = gbrealymin; gbgridymax = gbrealymax;
			gbgridzmin = gbrealzmin; gbgridzmax = gbrealzmax;
		elseif(gridrangeset .eq. use_currentrange) then
			if(print_info) print *, '  (do_cell_init) Set range of grids from current cosmology:'
			gbgridxmin = gbxmin; gbgridxmax = gbxmax
			gbgridymin = gbymin; gbgridymax = gbymax
			gbgridzmin = gbzmin; gbgridzmax = gbzmax
		elseif(gridrangeset .eq. use_noRSDrange) then
			if(print_info) print *, '  (do_cell_init) Set range of grids from current cosmology NORSD range:'
			gbgridxmin = gbnoRSDxmin; gbgridxmax = gbnoRSDxmax
			gbgridymin = gbnoRSDymin; gbgridymax = gbnoRSDymax
			gbgridzmin = gbnoRSDzmin; gbgridzmax = gbnoRSDzmax
		else
			print *, 'ERROR (do_cell_init)!!! Wrong gridrangeset: ', gridrangeset
			stop
		endif

		unit_len = (gbgridxmax - gbgridxmin) / rl_num_in_x
		if(print_info) write(*,'(20x,A,f12.5)')  'Width of cell: ', unit_len

		gb_n_cellx = int( (gbgridxmax-gbgridxmin)/unit_len)
		gb_n_celly = int( (gbgridymax-gbgridymin)/unit_len)
		gb_n_cellz = int( (gbgridzmax-gbgridzmin)/unit_len)

		gbdeltax = unit_len; gbdeltay = unit_len; gbdeltaz = unit_len

		gbgridxmax = gbgridxmin + gbdeltax*dble(gb_n_cellx)
		gbgridymax = gbgridymin + gbdeltay*dble(gb_n_celly)
		gbgridzmax = gbgridzmin + gbdeltaz*dble(gb_n_cellz)
		
		gbgridrmax = sqrt(gbgridxmax**2.0+gbgridymax**2.0+gbgridzmax**2.0)

		if(print_info) then
			write(*,'(20x,A,f12.5,2x,f12.5)')  'Grid Range of x: ', gbgridxmin, gbgridxmax
			write(*,'(20x,A,f12.5,2x,f12.5)')  'Grid Range of y: ', gbgridymin, gbgridymax
			write(*,'(20x,A,f12.5,2x,f12.5)')  'Grid Range of z: ', gbgridzmin, gbgridzmax
		endif
		
		gb_cell_vol = gbdeltax * gbdeltay * gbdeltaz
		if(print_info) then
			write(*,'(20x,A,f12.5)')  'Effective # of cell in x direction: ', real(rl_num_in_x)
			write(*,'(20x,A,i3,A,i3,A,i3,A,i11)')  'Grids of cell        = ', &
				gb_n_cellx,'*',gb_n_celly,'*',gb_n_cellz,' = ',gb_n_cellx*gb_n_celly*gb_n_cellz
			write(*,'(20x,A,3(f10.5,A),f12.5)') 'Each cell has size ', real(gbdeltax),' *',real(gbdeltay), &
				' *',real(gbdeltaz), ' =', real(gb_cell_vol)
		endif

		if(present(do_not_init_cellmat)) then
			if(do_not_init_cellmat .eq. .true.) return
		endif
		
		!allocating gb_cell_mat 
		if(allocated(gb_cell_mat)) deallocate(gb_cell_mat)
		allocate(gb_cell_mat(gb_n_cellx,gb_n_celly,gb_n_cellz),mark_num_halo(gb_n_cellx,gb_n_celly,gb_n_cellz))	
		
		! Count how many halos in each cell		
		do i = 1, gb_num_xyz_mass
			x= gb_xyz_list(1,i); y= gb_xyz_list(2,i); z= gb_xyz_list(3,i)
			if(x.ge.gbgridxmax .or. y.ge.gbgridymax .or. z.ge.gbgridzmax &
			   .or. x.le.gbgridxmin .or. y.le.gbgridymin .or. z.le.gbgridzmin) then
				cycle
			endif
			ix=max(min(int((x-gbgridxmin)/gbdeltax)+1,gb_n_cellx),1)
			iy=max(min(int((y-gbgridymin)/gbdeltay)+1,gb_n_celly),1)
			iz=max(min(int((z-gbgridzmin)/gbdeltaz)+1,gb_n_cellz),1)
			gb_cell_mat(ix,iy,iz)%halo_num = gb_cell_mat(ix,iy,iz)%halo_num + 1
		enddo
		
		! allocate halo_list in pixels having halos
		numpixel = gb_n_cellx*gb_n_celly*gb_n_cellz
		numhashalo = 0.0; maxhalonum = 0.0; avghalonum = 0.0
		do ix=1,gb_n_cellx
		do iy=1,gb_n_celly
		do iz=1,gb_n_cellz
			mark_num_halo(ix,iy,iz) = 0
			if(gb_cell_mat(ix,iy,iz)%halo_num > 0) then
				if(gb_cell_mat(ix,iy,iz)%halo_num > max_incell_halo_num) then
					print *, 'ERROR! incell halo num overflows!'
					print *, 'ix,iy,iz,num = ', ix,iy,iz,gb_cell_mat(ix,iy,iz)%halo_num
					stop
				endif
				allocate(gb_cell_mat(ix,iy,iz)%list(gb_cell_mat(ix,iy,iz)%halo_num))
				numhashalo = numhashalo + 1.0_dl
				avghalonum = avghalonum + gb_cell_mat(ix,iy,iz)%halo_num
				maxhalonum = max(maxhalonum,dble(gb_cell_mat(ix,iy,iz)%halo_num))
			endif
		enddo
		enddo
		enddo
		avghalonum = avghalonum / dble(numhashalo)
		
		! saving halos into halo lists
		do i = 1, gb_num_xyz_mass
			x= gb_xyz_list(1,i); y= gb_xyz_list(2,i); z= gb_xyz_list(3,i)
			if(x.ge.gbgridxmax .or. y.ge.gbgridymax .or. z.ge.gbgridzmax &
			   .or. x.le.gbgridxmin .or. y.le.gbgridymin .or. z.le.gbgridzmin) then
				cycle
			endif
			ix=max(min(int((x-gbgridxmin)/gbdeltax)+1,gb_n_cellx),1)
			iy=max(min(int((y-gbgridymin)/gbdeltay)+1,gb_n_celly),1)
			iz=max(min(int((z-gbgridzmin)/gbdeltaz)+1,gb_n_cellz),1)

			if(gb_cell_mat(ix,iy,iz)%halo_num > 0) then
				mark_num_halo(ix,iy,iz) = mark_num_halo(ix,iy,iz)+1
				gb_cell_mat(ix,iy,iz)%list(mark_num_halo(ix,iy,iz)) = i
			endif
		enddo
		
		! Check matching of halo #
		do ix=1,gb_n_cellx
		do iy=1,gb_n_celly
		do iz=1,gb_n_cellz
			if(mark_num_halo(ix,iy,iz) .ne. gb_cell_mat(ix,iy,iz)%halo_num) then
				print *, 'ERROR (do_cell_init)!!! halo # mismatch: '
				print *, 'ix,iy,iz, #1, #2 = ',ix,iy,iz,mark_num_halo(ix,iy,iz),gb_cell_mat(ix,iy,iz)%halo_num
				stop
			endif
		enddo
		enddo
		enddo
		deallocate(mark_num_halo)
		
		if(print_info) then
			write(*,'(18x,i7,A,f5.2,A,i5,A,f7.3)') int(numhashalo+0.5),'(',&
				real(numhashalo/numpixel)*100.0,'%) cells have halo inside; maximal # is ',&
				 int(maxhalonum+0.5), ', mean num is ', real(avghalonum)
			print *, '  (do_cell_init done)'
		endif
	end subroutine do_cell_init	

  !------------------------------------------
  ! central position of the cell
  !------------------------------------------	
  	subroutine cell_pos(ix,iy,iz,x,y,z)
		integer, intent(in) :: ix, iy, iz
  		real(dl), intent(out) :: x,y,z
  		x = gbgridxmin + (ix-0.5)*gbdeltax
  		y = gbgridymin + (iy-0.5)*gbdeltay
  		z = gbgridzmin + (iz-0.5)*gbdeltaz
	end subroutine cell_pos
	
  !------------------------------------------
  ! central position of the cell
  !------------------------------------------	
  	subroutine sad_pos(ix,iy,iz,x,y,z)
		integer, intent(in) :: ix, iy, iz
  		real(dl), intent(out) :: x,y,z
  		x = gbgridxmin + ix*gbdeltax
  		y = gbgridymin + iy*gbdeltay
  		z = gbgridzmin + iz*gbdeltaz
	end subroutine sad_pos
	real(dl) function cell_sad_dist(ixcell,iycell,izcell,ixsad,iysad,izsad)
		integer, intent(in) :: ixcell,iycell,izcell,ixsad,iysad,izsad
  		real(dl) :: dx,dy,dz
  		dx = (ixsad-ixcell+0.5)*gbdeltax
  		dy = (iysad-iycell+0.5)*gbdeltay
  		dz = (izsad-izcell+0.5)*gbdeltaz
  		cell_sad_dist = sqrt(dx**2.0+dy**2.0+dz**2.0)
	end function cell_sad_dist
	
  !------------------------------------------
  ! cell index 
  !------------------------------------------	
  	subroutine cell_index(ix,iy,iz,x,y,z)
		integer, intent(out) :: ix, iy, iz
  		real(dl), intent(in) :: x,y,z
  		ix = int((x-gbgridxmin) / gbdeltax +1.0)
  		iy = int((y-gbgridymin) / gbdeltay +1.0)
  		iz = int((z-gbgridzmin) / gbdeltaz +1.0)
	end subroutine cell_index


  !------------------------------------------
  ! Reset the mass_list: normalized by mass density;
  !------------------------------------------
	subroutine normn_mlist(printinfo, nbins, nsm)
		!Dummy
		integer, intent(in) :: nbins, nsm
		logical, intent(in) :: printinfo
		!Local
		real(dl), allocatable :: quan_av_list(:),quan_er_list(:),binned_r_list(:),quan_var_list(:),&
			quan_medval_list(:),numinbin(:),redgelist(:)
		integer, parameter :: polyorder = nfitrhocri
		real(dl) :: r, rmin, rmax, polyrho, polycoef(polyorder+1)
		integer :: i_sm, i
		
		if(printinfo) then
			write(*,'(A,i4,A,i3,A,2f10.3)'), '   (normn_mlist) Estimating medium rho for ', nbins, &
				' bins. Polyorder = ', polyorder ,'. Distance range: ', real(gbrmin), real(gbrmax)
		endif
		
		allocate(quan_av_list(nbins), quan_er_list(nbins), binned_r_list(nbins), &
			quan_var_list(nbins), quan_medval_list(nbins), numinbin(nbins), redgelist(nbins+1))
  		do i_sm = 1, nsm
  			call eqvl_binned_quan(gb_mass_list, gb_r_list, gbnoRSDrmin, gbnoRSDrmax, num_halo, nbins, &
		  		quan_av_list, quan_er_list, binned_r_list, numinbinlist= numinbin, redgelist = redgelist)
!		  	call poly_fit(binned_r_list,quan_medval_list,polycoef,nbins,polyorder)
			write(*,'(A)'), '   (normn_mlist) ATTENTION: Using mean value rather than middle value to normalize the mass!'
			! <rho>(r) = MASS / Vol = <M_i>*# / Vol
			do i = 1, nbins
				quan_av_list(i) = quan_av_list(i) * numinbin(i) / vol_fun(redgelist(i),redgelist(i+1))
			enddo
			if(printinfo) then
				write(*,'(19x,A,6x,A,10x,A,12x,A)') ' bin     #-in-bin     distance','redshift','   range of r ','density'
				do i = 1, nbins
					write(*,'(20x,i3,2x,i10,4(2x,f12.3),4x,e14.7)') i,int(numinbin(i)),binned_r_list(i),&
						de_zfromintpl(binned_r_list(i)), redgelist(i), redgelist(i+1), &
						quan_av_list(i)
				enddo
			endif
		  	call poly_fit(binned_r_list,quan_av_list,polycoef,nbins,polyorder)
		  	if(printinfo) then
		  		write(*,'(A,<polyorder+1>(e14.7,1x))') '                  Coefficients of poly: ', polycoef
		  	endif
	  		do i = 1, num_halo
  				r = gb_r_list(i)
  				polyrho = poly(r,polycoef,polyorder)
  				gb_mass_list(i) = gb_mass_list(i) / polyrho
  			enddo
  		enddo
  	end subroutine normn_mlist	

		
  !------------------------------------------
  ! Reset the mass_list: normalized by avg;
  !------------------------------------------
	subroutine medrho_mlist(distance_list, quan_list, num_quan, printinfo, nbins, nsm)
		!Dummy
		integer, intent(in) :: num_quan, nbins, nsm
		logical, intent(in) :: printinfo
		real(dl), intent(in) :: distance_list(num_quan), quan_list(num_quan)
		!Local
		real(dl), allocatable :: quan_av_list(:), quan_er_list(:), binned_r_list(:), quan_var_list(:), quan_medval_list(:)
		integer, parameter :: polyorder = nfitrhocri
		real(dl) :: r, rmin, rmax, polyrho, polycoef(polyorder+1)
		integer :: i_sm, i
		
		if(printinfo) then
			write(*,'(A,i4,A,i3,A,2f10.3)'), '   (medrho_mlist) Estimating medium rho for ', nbins, &
				' bins. Polyorder = ', polyorder ,'. Distance range: ', real(gbrmin), real(gbrmax)
		endif
		
		allocate(quan_av_list(nbins), quan_er_list(nbins), binned_r_list(nbins), &
			quan_var_list(nbins), quan_medval_list(nbins))
  		do i_sm = 1, nsm
  			call find_min_max(distance_list, num_quan, rmin, rmax)
  			call eqvl_binned_quan(quan_list, distance_list, rmin, rmax, num_quan, nbins, &
		  		quan_av_list, quan_er_list, binned_r_list, quan_var_list, quan_medval_list)
!		  	call poly_fit(binned_r_list,quan_medval_list,polycoef,nbins,polyorder)
			write(*,'(A)') '   (medrho_mlist) ATTENTION: Using mean value rather than middle value to normalize the mass!'
		  	call poly_fit(binned_r_list,quan_av_list,polycoef,nbins,polyorder)
		  	if(printinfo) then
		  		write(*,'(A,<polyorder+1>(e14.7,1x))') '                  Coefficients of poly: ', polycoef
		  	endif
	  		do i = 1, num_halo
  				r = gb_r_list(i)
  				polyrho = poly(r,polycoef,polyorder)
  				gb_mass_list(i) = gb_mass_list(i) / polyrho
  			enddo
  		enddo
  	end subroutine medrho_mlist		
end module ap_cell	
	

