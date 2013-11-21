


!####################################
!This module does smooth
!####################################
module ap_smooth
use ap_settings_init

	implicit none

! settings of ranges for grids.
	logical, public :: fixgridrange = .false.
	real(dl), public :: fixgridxmin, fixgridxmax, fixgridymin, fixgridymax, fixgridzmin, fixgridzmax
	logical, public :: use_realrange = .false. ! use the range of real data

! settings of the programe	
	real(dl) :: dft_ra_ratio = 4.0d0 ! 4.0d0 TESTING
	
! using fixed radius smoothing kernel
	logical, public :: use_fixmd = .false.
	real(dl), public :: gb_fixmd = 40.0_dl

! using cuts in smooth kernel (ignore nearby halos)
	logical, public :: use_seg = .false.
	real(dl), public :: seg_cut_ratio = 0.0_dl !-0.05_dl
	real(dl), public :: seg_cut_dist = 60.0_dl
	
!	real(dl) :: sep_distance = 17.065
!	real(dl) :: sd_unit_len = sep_distance * (dble(smooth_num)**(1.0/3.0)) ! unit length used to dividing cells, ...
!	logical :: use_sd_unit_len = .false.

!	
	
! nuisance settings	
!	character(len=char_len) :: rho_gradient_file_name
	integer :: gb_num_xyz_mass
	real(dl), allocatable :: gb_xyz_list(:,:), gb_mass_list(:), gb_bf_mass_list(:), gb_r_list(:)
	real(dl) :: gbtotvol
	real(dl) :: unit_len
	integer :: gb_n_cellx, gb_n_celly, gb_n_cellz
	real(dl) :: gbdeltax, gbdeltay, gbdeltaz
	real(dl) :: gb_cell_vol

	type :: halo_list
		integer :: halo_num = 0
		integer, allocatable :: list(:)
		real(dl) :: rho
		logical :: has_bd_effect
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
		if(allocated(gb_bf_mass_list)) &
			deallocate(gb_bf_mass_list)
		if(allocated(gb_cell_mat)) &
			deallocate(gb_cell_mat)
		if(allocated(gb_r_list)) &
			deallocate(gb_r_list)			
	end subroutine

  !------------------------------------------
  ! initialize the xyz_mass array
  !------------------------------------------
	subroutine init_xyz_r_gb_mass_list(RSD, AP, gv_print_info)
		integer :: RSD, AP, i
		real(dl) :: x, y, z, ratio, r
		logical, optional :: gv_print_info
		logical :: print_info

		if(present(gv_print_info)) then
			if(gv_print_info) &
			        write(*,'(A,i4,i4,A)'), '  (init_xyz_r_gb_mass_list) Initializing grids of cells with RSD, AP = ', RSD, AP, '  ...'
		endif

		call smooth_clean_up()
		
!		print *, 'Before allocating gb_xyz_list...'
		allocate(gb_xyz_list(3,num_halo), gb_mass_list(num_halo), gb_bf_mass_list(num_halo), gb_r_list(num_halo))
!		print *, 'After allocating gb_xyz_list...'
		gb_num_xyz_mass = num_halo
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
				gb_r_list(i) = halo_info(i)%r_AP_RSD(AP)
			else
				print *, 'ERROR! No such RSD, AP : ', RSD, AP
				stop
			endif
			ratio = gb_r_list(i) / r 
			gb_xyz_list(1,i) = x * ratio
			gb_xyz_list(2,i) = y * ratio
			gb_xyz_list(3,i) = z * ratio
			gb_mass_list(i) = halo_info(i)%mass
			gb_bf_mass_list(i) = gb_mass_list(i)
		enddo
		
		!min/max of x,y,z
		gbxmin = gb_xyz_list(1,1); gbxmax = gb_xyz_list(1,1);
		gbymin = gb_xyz_list(2,1); gbymax = gb_xyz_list(2,1);
		gbzmin = gb_xyz_list(3,1); gbzmax = gb_xyz_list(3,1);
		do i = 2, gb_num_xyz_mass
			x= gb_xyz_list(1,i); y= gb_xyz_list(2,i); z= gb_xyz_list(3,i)
			gbxmin = min(gbxmin,x); gbxmax = max(gbxmax,x)
			gbymin = min(gbymin,y); gbymax = max(gbymax,y)
			gbzmin = min(gbzmin,z); gbzmax = max(gbzmax,z)
		enddo
		!min/max of r
		call find_min_max(gb_r_list, num_halo, gbrmin, gbrmax)
		gbtotvol = vol_fun(gbrmin, gbrmax)
		
		if(print_info) then
			write(*,'(2x,A,f12.5,2x,f12.5)')  '  gbxmin / gbxmax = ', gbxmin, gbxmax
			write(*,'(2x,A,f12.5,2x,f12.5)')  '  gbymin / gbymax = ', gbymin, gbymax
			write(*,'(2x,A,f12.5,2x,f12.5)')  '  gbzmin / gbzmax = ', gbzmin, gbzmax
		endif
	end subroutine init_xyz_r_gb_mass_list
	
  !------------------------------------------
  ! read in the xyz_mass file
  !------------------------------------------
  	subroutine read_in_xyz_mass(xyz_mass_file_name)
  		character(len=char_len) :: xyz_mass_file_name
		integer :: i
		call count_line_number(xyz_mass_file_name, gb_num_xyz_mass)
		print*, '   num of halo = ', gb_num_xyz_mass, '...'
	!allocating arrays...
		if(allocated(gb_xyz_list)) &
			deallocate(gb_xyz_list)
		if(allocated(gb_mass_list)) &			
			deallocate(gb_mass_list)
		allocate(gb_xyz_list(3,gb_num_xyz_mass))
		allocate(gb_mass_list(gb_num_xyz_mass))
		open(unit=1,file=xyz_mass_file_name)
		do i = 1, gb_num_xyz_mass
			read(1,*) gb_xyz_list(1:3,i), gb_mass_list(i)
		enddo
	end subroutine read_in_xyz_mass


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
	subroutine do_cell_init(RSD, AP, rl_num_in_x, gv_print_info, do_not_init_xyz_mass, do_not_init_cellmat)
		integer :: RSD, AP, num_in_x
		real(dl) :: rl_num_in_x
		logical, optional :: gv_print_info, do_not_init_xyz_mass, do_not_init_cellmat
		logical :: print_info
		integer :: i, j, k, ix, iy, iz
		real(dl) :: x,y,z
		integer :: now_halo_num
		!maximal length of array to save index of halos in one cell
		integer, parameter :: max_incell_halo_num = 20000
		integer, parameter :: re_alo_num = 100

		real(dl), allocatable :: tmp_list(:)
		
		if(present(gv_print_info)) then
			print_info = gv_print_info
			else
			print_info = .false.
		endif

		if(present(do_not_init_xyz_mass)) then
			if(do_not_init_xyz_mass) then
				continue
				else
				call init_xyz_r_gb_mass_list(RSD, AP, print_info)
			endif
		else
			call init_xyz_r_gb_mass_list(RSD, AP, print_info)
		endif
		
		! For grid, use the same value as gbxmin/gbymin/gbzmin
! BEGIN TESTING

		if(fixgridrange) then
			gbgridxmin = fixgridxmin; gbgridxmax = fixgridxmax
			gbgridymin = fixgridymin; gbgridymax = fixgridymax
			gbgridzmin = fixgridzmin; gbgridzmax = fixgridzmax
		else			
			if(use_realrange .or. dotsbe) then
				gbgridxmin = gbrealxmin; gbgridxmax = gbrealxmax;
				gbgridymin = gbrealymin; gbgridymax = gbrealymax;
				gbgridzmin = gbrealzmin; gbgridzmax = gbrealzmax;
			else
				gbgridxmin = gbxmin; gbgridxmax = gbxmax
				gbgridymin = gbymin; gbgridymax = gbymax
				gbgridzmin = gbzmin; gbgridzmax = gbzmax
			endif
		endif

! END TESTING
		
	!basic info for cell division	
!		if(use_sd_unit_len) then
!			unit_len = sd_unit_len
!		else
		unit_len = (gbgridxmax - gbgridxmin) / rl_num_in_x
!		endif

		gb_n_cellx = int( (gbgridxmax-gbgridxmin)/unit_len + 0.5)
		gb_n_celly = int( (gbgridymax-gbgridymin)/unit_len + 0.5)
		gb_n_cellz = int( (gbgridzmax-gbgridzmin)/unit_len + 0.5)

		gbdeltax = unit_len
		gbdeltay = unit_len
		gbdeltaz = unit_len

		gbgridxmax = gbgridxmin + gbdeltax*dble(gb_n_cellx)
		gbgridymax = gbgridymin + gbdeltay*dble(gb_n_celly)
		gbgridzmax = gbgridzmin + gbdeltaz*dble(gb_n_cellz)

		if(print_info) then
			write(*,'(2x,A,f12.5,2x,f12.5)')  '  gbgridxmin / gbgridxmax = ', gbgridxmin, gbgridxmax
			write(*,'(2x,A,f12.5,2x,f12.5)')  '  gbgridymin / gbgridymax = ', gbgridymin, gbgridymax
			write(*,'(2x,A,f12.5,2x,f12.5)')  '  gbgridzmin / gbgridzmax = ', gbgridzmin, gbgridzmax
		endif
		
		gb_cell_vol = gbdeltax * gbdeltay * gbdeltaz
		if(print_info) then
			print*, 'rl_num_in_x =', rl_num_in_x, '...'
			print*,  '   gb_n_cellx, gbdeltax = ', gb_n_cellx, gbdeltax
			print*,  '   gb_n_celly, gbdeltay = ', gb_n_celly, gbdeltay
			print*,  '   gb_n_cellz, gbdeltaz = ', gb_n_cellz, gbdeltaz
			print*,  '   gb_cell_vol = ', gb_cell_vol
		endif

		if(present(do_not_init_cellmat)) then
			if(do_not_init_cellmat .eq. .true.) return
		endif
		
	!allocating gb_cell_mat 
		if(allocated(gb_cell_mat)) &
			deallocate(gb_cell_mat)
		allocate(gb_cell_mat(gb_n_cellx,gb_n_celly,gb_n_cellz))	
		allocate(tmp_list(max_incell_halo_num))
				
		do i = 1, gb_num_xyz_mass
			x= gb_xyz_list(1,i); y= gb_xyz_list(2,i); z= gb_xyz_list(3,i)
			if(x.ge.gbgridxmax .or. y.ge.gbgridymax .or. z.ge.gbgridzmax &
			   .or. x.le.gbgridxmin .or. y.le.gbgridymin .or. z.le.gbgridzmin) then
				cycle
			endif
			ix=max(min(int((x-gbgridxmin)/gbdeltax)+1,gb_n_cellx),1)
			iy=max(min(int((y-gbgridymin)/gbdeltay)+1,gb_n_celly),1)
			iz=max(min(int((z-gbgridzmin)/gbdeltaz)+1,gb_n_cellz),1)
			
			now_halo_num = gb_cell_mat(ix,iy,iz)%halo_num
			
			if(now_halo_num .eq. max_incell_halo_num) then
				print *, 'ERROR! incell halo num overflows!'
!				print *, 'ix,iy,iz,num = ', ix,iy,iz,num
!				print *, 'x,y,z = ', 
				stop
			endif
			
			if(mod(now_halo_num,re_alo_num) .eq. 0) then
				if(now_halo_num .eq. 0) then
					allocate(gb_cell_mat(ix,iy,iz)%list(re_alo_num))
				else
!					print*, 'ix,iy,iz,now_halo_num = ', ix,iy,iz,now_halo_num
					tmp_list(1:now_halo_num) = gb_cell_mat(ix,iy,iz)%list(1:now_halo_num)
					deallocate(gb_cell_mat(ix,iy,iz)%list)
					allocate(gb_cell_mat(ix,iy,iz)%list(now_halo_num+re_alo_num))
					gb_cell_mat(ix,iy,iz)%list(1:now_halo_num) = tmp_list(1:now_halo_num)
				endif
			endif
			
			now_halo_num = now_halo_num + 1
			gb_cell_mat(ix,iy,iz)%halo_num = now_halo_num
			gb_cell_mat(ix,iy,iz)%list(now_halo_num) = i
		enddo
		deallocate(tmp_list)
		if(print_info) then
			print *, '  Cell initialization done.'
		endif
	end subroutine do_cell_init	
	
	
  !------------------------------------------
  ! select out halos in the cells nearby the
  !  given position. make sure the number is 
  !  far larger than the required
  !------------------------------------------
  	subroutine nb_select(x,y,z,num,given_ra_ratio,selected_list)
  		real(dl), intent(in) :: x,y,z
  		real(dl), optional, intent(in) :: given_ra_ratio
  		integer, intent(in) :: num
  		integer, allocatable, intent(out) :: selected_list(:)
  		integer :: i,j,k,l1,l2,ix,iy,iz,requirednum
  		integer :: di, i1,i2, j1,j2, k1,k2, now_num
  		real(dl) :: dev, ra_ratio, ratio
  		ix = max(min(int((x-gbgridxmin)/gbdeltax)+1,gb_n_cellx),1)
  		iy = max(min(int((y-gbgridymin)/gbdeltay)+1,gb_n_celly),1)
  		iz = max(min(int((z-gbgridzmin)/gbdeltaz)+1,gb_n_cellz),1)
  		dev = max(abs(abs(gbgridxmin + (ix-1)*gbdeltax - x)/gbdeltax - 0.5), &
  			abs(abs(gbgridymin + (iy-1)*gbdeltay - y)/gbdeltay - 0.5), &
  			abs(abs(gbgridzmin + (iz-1)*gbdeltaz - z)/gbdeltaz - 0.5) )
  		if(present(given_ra_ratio)) then
  			ra_ratio = given_ra_ratio
  		else
  			ra_ratio = dft_ra_ratio
  		endif
  		ratio = ra_ratio * (1.0 + dev * 1.5)
  		requirednum = max(int(num*ratio), 40) ! at least 40 halos!!!
!  		print *, 'ix, iy, iz, dev, ratio, requirednum = ', ix, iy, iz, dev, ratio, requirednum
		di = 0
		do while(1.eq.1)
			call ilist(ix,di,1,gb_n_cellx,i1,i2)
			call ilist(iy,di,1,gb_n_celly,j1,j2)
			call ilist(iz,di,1,gb_n_cellz,k1,k2)
			now_num = 0
			do i = i1, i2
			do j = j1, j2
			do k = k1, k2
				now_num = now_num + gb_cell_mat(i,j,k)%halo_num
			enddo
			enddo
			enddo
			if(now_num .ge. requirednum) exit
			di = di + 1
		enddo
		allocate(selected_list(now_num))
		l1 = 1
		do i = i1, i2
		do j = j1, j2
		do k = k1, k2
			do l2 = 1, gb_cell_mat(i,j,k)%halo_num
				selected_list(l1) = gb_cell_mat(i,j,k)%list(l2)
				l1 = l1 + 1
			enddo
		enddo
		enddo
		enddo
	end subroutine nb_select
		

  !------------------------------------------
  ! select out halos in the cells nearby the
  !  given position. make sure the number is 
  !  far larger than the required
  !------------------------------------------
  	subroutine nb_select2(x,y,z,num,selected_list)
  		real(dl), intent(in) :: x,y,z
  		integer, intent(in) :: num
 ! 		real(dl), optional, intent(in) :: given_ra_ratio  		
  		integer, allocatable, intent(out) :: selected_list(:)
  		integer :: i,j,k,l1,l2,ix,iy,iz,requirednum
  		integer :: di, i1,i2, j1,j2, k1,k2, now_num
  		real(dl) :: dev, ra_ratio, ratio
  		
  		ix = max(min(int((x-gbgridxmin)/gbdeltax+0.01)+1,gb_n_cellx),1)
  		iy = max(min(int((y-gbgridymin)/gbdeltay+0.01)+1,gb_n_celly),1)
  		iz = max(min(int((z-gbgridzmin)/gbdeltaz+0.01)+1,gb_n_cellz),1)
  		dev = max(abs(abs(gbgridxmin + (ix-1)*gbdeltax - x)/gbdeltax - 0.5), &
  			abs(abs(gbgridymin + (iy-1)*gbdeltay - y)/gbdeltay - 0.5), &
  			abs(abs(gbgridzmin + (iz-1)*gbdeltaz - z)/gbdeltaz - 0.5) )
		
!		print *, 'ix, iy, iz, dev = ', ix, iy, iz, dev ! lxd
!		stop

  		requirednum = num
!  		print *, 'ix, iy, iz, dev, ratio, requirednum = ', ix, iy, iz, dev, ratio, requirednum
		di = 0
		do while(1.eq.1)
			call ilist(ix,di,1,gb_n_cellx,i1,i2)
			call ilist(iy,di,1,gb_n_celly,j1,j2)
			call ilist(iz,di,1,gb_n_cellz,k1,k2)
			now_num = 0
			do i = i1, i2
			do j = j1, j2
			do k = k1, k2
				now_num = now_num + gb_cell_mat(i,j,k)%halo_num
			enddo
			enddo
			enddo
			if(now_num .ge. requirednum) exit
			di = di + 1
		enddo
		di = ceiling((di+0.5)*sqrt(3.0)-0.5)
		call ilist(ix,di,1,gb_n_cellx,i1,i2)
		call ilist(iy,di,1,gb_n_celly,j1,j2)
		call ilist(iz,di,1,gb_n_cellz,k1,k2)
		now_num = 0
		do i = i1, i2
		do j = j1, j2
		do k = k1, k2
			now_num = now_num + gb_cell_mat(i,j,k)%halo_num
		enddo
		enddo
		enddo
		allocate(selected_list(now_num))
		l1 = 1
		do i = i1, i2
		do j = j1, j2
		do k = k1, k2
			do l2 = 1, gb_cell_mat(i,j,k)%halo_num
				selected_list(l1) = gb_cell_mat(i,j,k)%list(l2)
				l1 = l1 + 1
			enddo
		enddo
		enddo
		enddo
	end subroutine nb_select2

  !------------------------------------------
  ! select out cells which are surrounded to 
  !  a center, with distance di
  ! not efficient, abondoned!
  !------------------------------------------
	subroutine nbcell_list(ix,iy,iz,di,narray,istart,iend,nbclist)
		! dummy
		integer, intent(in) :: ix, iy, iz, di, narray,istart
		integer, intent(out) :: iend,nbclist(3,narray)
		! local
		integer :: n, i, j, k, nowi 

		if(di .eq. 0) then
			iend = istart
			nbclist(1,istart) = ix
			nbclist(2,istart) = iy
			nbclist(3,istart) = iz
			nowi = istart+1
		else
			iend = istart + (2*di*4)*(2*di+1) + 2*(2*di-1)**2 - 1
			nowi = istart			
			do k = iz-di, iz+di
				i = ix+di
				do j = iy-di,iy+di
					nbclist(1,nowi) = i
					nbclist(2,nowi) = j
					nbclist(3,nowi) = k
					nowi = nowi + 1
				enddo
				i = ix-di
				do j = iy-di,iy+di
					nbclist(1,nowi) = i
					nbclist(2,nowi) = j
					nbclist(3,nowi) = k
					nowi = nowi + 1
				enddo
				j = iy-di
				do i = ix-di+1,ix+di-1
					nbclist(1,nowi) = i
					nbclist(2,nowi) = j
					nbclist(3,nowi) = k
					nowi = nowi + 1
				enddo
				j = iy+di
				do i = ix-di+1,ix+di-1
					nbclist(1,nowi) = i
					nbclist(2,nowi) = j
					nbclist(3,nowi) = k
					nowi = nowi + 1
				enddo
			enddo
			k = iz-di
			do i = ix-di+1, ix+di-1
				do j = iy-di+1, iy+di-1
					nbclist(1,nowi) = i
					nbclist(2,nowi) = j
					nbclist(3,nowi) = k
					nowi = nowi + 1
				enddo
			enddo
			k = iz+di
			do i = ix-di+1, ix+di-1
				do j = iy-di+1, iy+di-1
					nbclist(1,nowi) = i
					nbclist(2,nowi) = j
					nbclist(1,nowi) = k
					nowi = nowi + 1
				enddo
			enddo
		endif
		
!		print *, 'di, istart     = ', di, istart
!		print *, 'iend, nowi-1   = ', iend, nowi-1
!		print *, '  List of coords:'
!		do i = istart, iend
!			print *, '     nowi: ', i, '; coord: ', nbclist(1:3,i)
!		enddo
!		print *
	end subroutine nbcell_list


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
  ! estimating rho and gradient rho based on
  !  cubic spline kernel; fixed radius
  !------------------------------------------
  	subroutine nb_fixmd_list(x,y,z,num,rho,drhodx,drhody,drhodz)
		! DUMMY ARGUMENTS
  		real(dl), intent(in) :: x,y,z
  		integer, intent(out) :: num
  		real(dl), intent(out) :: rho,drhodx,drhody,drhodz
  		! LOCAL VARIABLES
  		integer, parameter :: max_num = 10000 !maximal # of neary halos 
  		real(dl) :: distance_array(max_num), xyz_mass_array(4,max_num)
  		integer :: ix,iy,iz, di, i,j,k,l, nownum,nowindex
  		real(dl) :: r0(3), h, nowr, mass, dweight

  		ix = max(min(int((x-gbgridxmin)/gbdeltax)+1,gb_n_cellx),1)
  		iy = max(min(int((y-gbgridymin)/gbdeltay)+1,gb_n_celly),1)
  		iz = max(min(int((z-gbgridzmin)/gbdeltaz)+1,gb_n_cellz),1)

		di = floor(gb_fixmd / gbdeltax + 0.5_dl)
		
!		if(.true.) then
!			print *, gb_fixmd / gbdeltax, di
!			print_di = .false.
!		endif
		
		num = 0
		r0(1)=x; r0(2)=y; r0(3)=z;
		do i = ix-di, ix+di
		do j = iy-di, iy+di
		do k = iz-di, iz+di
			do l = 1, gb_cell_mat(i,j,k)%halo_num
				nowindex = gb_cell_mat(i,j,k)%list(l)
				nowr = distance(gb_xyz_list(1:3,nowindex),r0)
				if(nowr < gb_fixmd) then
					num = num+1
					distance_array(num) = nowr
					xyz_mass_array(1:3, num) = gb_xyz_list(1:3,nowindex)
					xyz_mass_array(4,num) = gb_mass_list(nowindex)
				endif
			enddo
		enddo
		enddo
		enddo

		if(num > max_num) then
			print *, 'ERROR (nb_fixmd_list): # of halos overflow: ', num, max_num
			stop
		endif
		
		if(num .eq. 0) return
			
  		h = gb_fixmd / 2.0
		rho = 0; drhodx=0;drhody=0;drhodz=0;
		do i = 1, num
			nowr = distance_array(i)
			mass = xyz_mass_array(4,i)
!			print *, i, mass !TESTING			
			rho = rho + mass*w_kernel(nowr, h)
			dweight = der_w_kernel(nowr,h)
			drhodx = drhodx + mass*(x-xyz_mass_array(1,i)) / nowr * dweight
			drhody = drhody + mass*(y-xyz_mass_array(2,i)) / nowr * dweight
			drhodz = drhodz + mass*(z-xyz_mass_array(3,i)) / nowr * dweight
		enddo
!		stop !TESTING
	end subroutine nb_fixmd_list


  !------------------------------------------
  ! estimating rho based on
  !  cubic spline kernel; fixed radius
  ! Same as above, but only calc rho
  !------------------------------------------
  	subroutine nb_fixmd_rho(x,y,z,num,rho)
		! DUMMY ARGUMENTS
  		real(dl), intent(in) :: x,y,z
  		integer, intent(out) :: num
  		real(dl), intent(out) :: rho
  		! LOCAL VARIABLES
  		integer, parameter :: max_num = 10000 !maximal # of neary halos 
  		real(dl) :: distance_array(max_num), xyz_mass_array(4,max_num)
  		integer :: ix,iy,iz, di, i,j,k,l, nownum,nowindex
  		real(dl) :: r0(3), h, nowr, mass

  		ix = max(min(int((x-gbgridxmin)/gbdeltax)+1,gb_n_cellx),1)
  		iy = max(min(int((y-gbgridymin)/gbdeltay)+1,gb_n_celly),1)
  		iz = max(min(int((z-gbgridzmin)/gbdeltaz)+1,gb_n_cellz),1)

		di = floor(gb_fixmd / gbdeltax + 0.5_dl)
		
		num = 0
		r0(1)=x; r0(2)=y; r0(3)=z;
		do i = ix-di, ix+di
		do j = iy-di, iy+di
		do k = iz-di, iz+di
			do l = 1, gb_cell_mat(i,j,k)%halo_num
				nowindex = gb_cell_mat(i,j,k)%list(l)
				nowr = distance(gb_xyz_list(1:3,nowindex),r0)
				if(nowr < gb_fixmd) then
					num = num+1
					distance_array(num) = nowr
					xyz_mass_array(1:3, num) = gb_xyz_list(1:3,nowindex)
					xyz_mass_array(4,num) = gb_mass_list(nowindex)
				endif
			enddo
		enddo
		enddo
		enddo

		if(num > max_num) then
			print *, 'ERROR (nb_fixmd_list): # of halos overflow: ', num, max_num
			stop
		endif
		
		if(num .eq. 0) return
			
  		h = gb_fixmd / 2.0
		rho = 0
		do i = 1, num
			nowr = distance_array(i)
			mass = xyz_mass_array(4,i)
!			print *, i, mass !TESTING			
			rho = rho + mass*w_kernel(nowr, h)
		enddo
!		stop !TESTING
	end subroutine nb_fixmd_rho


  !------------------------------------------
  ! estimating rho and gradient rho based on
  !  cubic spline kernel
  !------------------------------------------
  	subroutine nb_list0(x,y,z,num,rho,drhodx,drhody,drhodz,max_dist)
		! DUMMY ARGUMENTS
  		real(dl), intent(in) :: x,y,z
  		integer, intent(in) :: num
  		real(dl), intent(out) :: rho,drhodx,drhody,drhodz,max_dist
  		! LOCAL VARIABLES
  		integer, allocatable :: index_array(:)
  		real(dl), allocatable :: distance_array(:), xyz_mass_array(:,:), tmp(:), smda(:)
  		integer, allocatable :: tmpindex(:), smlablist(:)
  		integer :: i,j,n,nowindex
  		real(dl) :: h, r0(3), r, mass, dweight

  		call nb_select(x,y,z,num,selected_list=index_array)
  		n=size(index_array)

  		allocate(tmpindex(n), distance_array(n))
		r0(1)=x; r0(2)=y; r0(3)=z;
		do i = 1, n
			distance_array(i) = distance(gb_xyz_list(1:3,index_array(i)),r0)
			tmpindex(i) = index_array(i)
		enddo

		allocate(smlablist(num),smda(num))
		call ltlablist2(distance_array,n,num,smda,smlablist)

		allocate(xyz_mass_array(4,num))
		
		do i = 1, num
			nowindex = tmpindex(smlablist(i))
			xyz_mass_array(1:3,i) = gb_xyz_list(1:3,nowindex)
			xyz_mass_array(4,i) = gb_mass_list(nowindex)
!		  		print *, i, real(xyz_mass_array(1:5,i))
 		enddo

  		max_dist = maxval(smda)
  		h = max_dist / 2.0
		rho = 0; drhodx=0;drhody=0;drhodz=0;
		do i = 1, num
			r = distance_array(smlablist(i))
			mass = xyz_mass_array(4,i)
!			print *, i, mass !TESTING			
			rho = rho + mass*w_kernel(r, h)
			dweight = der_w_kernel(r,h)
			drhodx = drhodx + mass*(x-xyz_mass_array(1,i)) / r * dweight
			drhody = drhody + mass*(y-xyz_mass_array(2,i)) / r * dweight
			drhodz = drhodz + mass*(z-xyz_mass_array(3,i)) / r * dweight
		enddo
!		stop !TESTING
	end subroutine nb_list0

  !------------------------------------------
  ! estimating rho and gradient rho based on
  !  cubic spline kernel
  !------------------------------------------
  	subroutine nb_seg_list(x,y,z,num,rho,drhodx,drhody,drhodz,max_dist)
  		! Dummy
  		real(dl), intent(in) :: x,y,z
  		integer, intent(in) :: num
  		real(dl), intent(out) :: rho,drhodx,drhody,drhodz,max_dist
  		! Local
  		integer, allocatable :: index_array(:)
  		real(dl), allocatable :: distance_array(:), xyz_mass_array(:,:)
  		integer :: i,nowindex,n
  		real(dl) :: h, r0(3), r, mass, dweight
  		call nb_select(x,y,z,num,selected_list = index_array)
  		n  = size(index_array)
  		allocate(xyz_mass_array(4,n), distance_array(n))
  		r0(1)=x; r0(2)=y; r0(3)=z;
  		do i = 1, n
  			distance_array(i) = distance(gb_xyz_list(1:3,index_array(i)),r0)
  		enddo
  		call Qsort2(distance_array,index_array,n)
  		max_dist = distance_array(num)
		if(distance_array(num-5) < seg_cut_dist) then
			max_dist = 1.0e10
			return
		endif 
!		do i = 1, num
!			if(distance_array(i)>seg_cut_dist) exit
!		enddo
!		if(num-i+1 .lt. 5) then
!			max_dist = 1.0e10
!			return
!		endif
		h = max_dist / 2.0
		rho = 0; drhodx=0;drhody=0;drhodz=0;
!		do i = int(num*seg_cut_ratio+1.5), num
		do i = 1, num
			r = distance_array(i)
			if(r < seg_cut_dist) cycle
			nowindex = index_array(i)
			mass = gb_mass_list(nowindex)
			rho = rho + mass*w_kernel(r, h)
			dweight = der_w_kernel(r,h)
			drhodx = drhodx + mass*(x-gb_xyz_list(1,nowindex)) / r * dweight
			drhody = drhody + mass*(y-gb_xyz_list(2,nowindex)) / r * dweight
			drhodz = drhodz + mass*(z-gb_xyz_list(3,nowindex)) / r * dweight
		enddo
	end subroutine nb_seg_list
		
end module ap_smooth	
	

