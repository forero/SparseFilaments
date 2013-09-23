


!####################################
!This module does smooth
!####################################
module ap_smooth
use ap_settings_init

	implicit none
	
! settings of the programe	
	real(dl) :: dft_ra_ratio = 4.0d0
	
!	real(dl) :: sep_distance = 17.065
!	real(dl) :: sd_unit_len = sep_distance * (dble(smooth_num)**(1.0/3.0)) ! unit length used to dividing cells, ...
!	logical :: use_sd_unit_len = .false.
	
	
! nuisance settings	
!	character(len=char_len) :: rho_gradient_file_name
	integer :: gb_num_xyz_mass
	real(dl), allocatable :: gb_xyz_list(:,:), gb_mass_list(:), gb_bf_mass_list(:), gb_r_list(:)
	real(dl) :: gbxmin, gbxmax, gbymin, gbymax, gbzmin, gbzmax
	real(dl) :: gbrmin, gbrmax, gbtotvol
	real(dl) :: unit_len
	integer :: gb_n_cellx, gb_n_celly, gb_n_cellz
	real(dl) :: gbdeltax, gbdeltay, gbdeltaz
	real(dl) :: gb_cell_vol

	type :: halo_list
		integer :: halo_num = 0
		integer, allocatable :: list(:)
	end type
	
	! array containing which halo in which cell
	type(halo_list), allocatable :: gb_cell_mat(:,:,:)
	
	type :: pixelinfo
		real(dl), allocatable :: xyzrlist(:,:)
		integer, allocatable :: indexlist(:)
		real(dl) :: maxdist
	end type	
	
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
			        write(*,'(A,i4,i4,A)'), '  Initializing grids of cells with RSD, AP = ', RSD, AP, '  ...'
		endif

		call smooth_clean_up()

		allocate(gb_xyz_list(3,num_halo), gb_mass_list(num_halo), gb_bf_mass_list(num_halo), gb_r_list(num_halo))
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
	subroutine do_cell_initialize(RSD, AP, num_in_x, gv_print_info, do_not_init_xyz_mass)
		integer :: RSD, AP, num_in_x
		logical, optional :: gv_print_info, do_not_init_xyz_mass
		logical :: print_info
		integer :: i, j, k, ix, iy, iz
		real(dl) :: x,y,z
		integer :: now_halo_num
		!maximal length of array to save index of halos in one cell
		integer, parameter :: max_incell_halo_num = 10000
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
			write(*,'(2x,A,f10.5,2x,f10.5)')  '  gbxmin / gbxmax = ', gbxmin, gbxmax
			write(*,'(2x,A,f10.5,2x,f10.5)')  '  gbymin / gbymax = ', gbymin, gbymax
			write(*,'(2x,A,f10.5,2x,f10.5)')  '  gbzmin / gbzmax = ', gbzmin, gbzmax
		endif
		
	!basic info for cell division	
!		if(use_sd_unit_len) then
!			unit_len = sd_unit_len
!		else
			unit_len = (gbxmax - gbxmin) / dble(num_in_x)
!		endif

		gb_n_cellx = int( (gbxmax-gbxmin)/unit_len + 0.5)
		gb_n_celly = int( (gbymax-gbymin)/unit_len + 0.5)
		gb_n_cellz = int( (gbzmax-gbzmin)/unit_len + 0.5)

		gbdeltax = (gbxmax-gbxmin) / dble(gb_n_cellx)
		gbdeltay = (gbymax-gbymin) / dble(gb_n_celly)
		gbdeltaz = (gbzmax-gbzmin) / dble(gb_n_cellz)

		gb_cell_vol = gbdeltax * gbdeltay * gbdeltaz

		if(print_info) then
			print*,  '   gb_n_cellx, gbdeltax = ', gb_n_cellx, gbdeltax
			print*,  '   gb_n_celly, gbdeltay = ', gb_n_celly, gbdeltay
			print*,  '   gb_n_cellz, gbdeltaz = ', gb_n_cellz, gbdeltaz
			print*,  '   gb_cell_vol = ', gb_cell_vol
		endif

		
	!allocating gb_cell_mat 
		if(allocated(gb_cell_mat)) &
			deallocate(gb_cell_mat)
		allocate(gb_cell_mat(gb_n_cellx,gb_n_celly,gb_n_cellz))	
		allocate(tmp_list(max_incell_halo_num))
				
		do i = 1, gb_num_xyz_mass
			x= gb_xyz_list(1,i); y= gb_xyz_list(2,i); z= gb_xyz_list(3,i)
			ix=max(min(int((x-gbxmin)/gbdeltax)+1,gb_n_cellx),1)
			iy=max(min(int((y-gbymin)/gbdeltay)+1,gb_n_celly),1)
			iz=max(min(int((z-gbzmin)/gbdeltaz)+1,gb_n_cellz),1)
			
			now_halo_num = gb_cell_mat(ix,iy,iz)%halo_num
			
			if(now_halo_num .eq. max_incell_halo_num) then
				print *, 'ERROR! incell halo num overflows'
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
	end subroutine
	
  !------------------------------------------
  ! initialize the Calculation
  !------------------------------------------
	subroutine do_cell_init2(RSD, AP, rl_num_in_x, gv_print_info, do_not_init_xyz_mass)
		integer :: RSD, AP, num_in_x
		real(dl) :: rl_num_in_x
		logical, optional :: gv_print_info, do_not_init_xyz_mass
		logical :: print_info
		integer :: i, j, k, ix, iy, iz
		real(dl) :: x,y,z
		integer :: now_halo_num
		!maximal length of array to save index of halos in one cell
		integer, parameter :: max_incell_halo_num = 10000
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
			write(*,'(2x,A,f10.5,2x,f10.5)')  '  gbxmin / gbxmax = ', gbxmin, gbxmax
			write(*,'(2x,A,f10.5,2x,f10.5)')  '  gbymin / gbymax = ', gbymin, gbymax
			write(*,'(2x,A,f10.5,2x,f10.5)')  '  gbzmin / gbzmax = ', gbzmin, gbzmax
		endif
		
	!basic info for cell division	
!		if(use_sd_unit_len) then
!			unit_len = sd_unit_len
!		else
		unit_len = (gbxmax - gbxmin) / rl_num_in_x
!		endif

		gb_n_cellx = int( (gbxmax-gbxmin)/unit_len + 0.5)
		gb_n_celly = int( (gbymax-gbymin)/unit_len + 0.5)
		gb_n_cellz = int( (gbzmax-gbzmin)/unit_len + 0.5)

		gbdeltax = unit_len
		gbdeltay = unit_len
		gbdeltaz = unit_len

		gb_cell_vol = gbdeltax * gbdeltay * gbdeltaz
		if(print_info) then
			print*, 'rl_num_in_x =', rl_num_in_x, '...'
			print*,  '   gb_n_cellx, gbdeltax = ', gb_n_cellx, gbdeltax
			print*,  '   gb_n_celly, gbdeltay = ', gb_n_celly, gbdeltay
			print*,  '   gb_n_cellz, gbdeltaz = ', gb_n_cellz, gbdeltaz
			print*,  '   gb_cell_vol = ', gb_cell_vol
		endif

		
	!allocating gb_cell_mat 
		if(allocated(gb_cell_mat)) &
			deallocate(gb_cell_mat)
		allocate(gb_cell_mat(gb_n_cellx,gb_n_celly,gb_n_cellz))	
		allocate(tmp_list(max_incell_halo_num))
				
		do i = 1, gb_num_xyz_mass
			x= gb_xyz_list(1,i); y= gb_xyz_list(2,i); z= gb_xyz_list(3,i)
			ix=max(min(int((x-gbxmin)/gbdeltax)+1,gb_n_cellx),1)
			iy=max(min(int((y-gbymin)/gbdeltay)+1,gb_n_celly),1)
			iz=max(min(int((z-gbzmin)/gbdeltaz)+1,gb_n_cellz),1)
			
			now_halo_num = gb_cell_mat(ix,iy,iz)%halo_num
			
			if(now_halo_num .eq. max_incell_halo_num) then
				print *, 'ERROR! incell halo num overflows'
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
	end subroutine do_cell_init2
	
	
	
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
  		ix = max(min(int((x-gbxmin)/gbdeltax)+1,gb_n_cellx),1)
  		iy = max(min(int((y-gbymin)/gbdeltay)+1,gb_n_celly),1)
  		iz = max(min(int((z-gbzmin)/gbdeltaz)+1,gb_n_cellz),1)
  		dev = max(abs(abs(gbxmin + (ix-1)*gbdeltax - x)/gbdeltax - 0.5), &
  			abs(abs(gbymin + (iy-1)*gbdeltay - y)/gbdeltay - 0.5), &
  			abs(abs(gbzmin + (iz-1)*gbdeltaz - z)/gbdeltaz - 0.5) )
  		if(present(given_ra_ratio)) then
  			ra_ratio = given_ra_ratio
  		else
  			ra_ratio = dft_ra_ratio
  		endif
  		ratio = ra_ratio * (1.0 + dev * 1.5)
  		requirednum = int(num*ratio)
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
  ! central position of the cell
  !------------------------------------------	
  	subroutine cell_pos(ix,iy,iz,x,y,z)
		integer ix, iy, iz
  		real(dl) :: x,y,z
  		x = gbxmin + (ix-0.5)*gbdeltax
  		y = gbymin + (iy-0.5)*gbdeltay
  		z = gbzmin + (iz-0.5)*gbdeltaz
	end subroutine cell_pos


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
			rho = rho + mass*w_kernel(r, h)
			dweight = der_w_kernel(r,h)
			drhodx = drhodx + mass*(x-xyz_mass_array(1,i)) / r * dweight
			drhody = drhody + mass*(y-xyz_mass_array(2,i)) / r * dweight
			drhodz = drhodz + mass*(z-xyz_mass_array(3,i)) / r * dweight
		enddo
	end subroutine nb_list0

  !------------------------------------------
  ! estimating rho and gradient rho based on
  !  cubic spline kernel
  !------------------------------------------
  	subroutine nb_listinput(x,y,z,num,rho,drhodx,drhody,drhodz,pixel)
		! DUMMY ARGUMENTS
  		real(dl), intent(in) :: x,y,z
  		integer, intent(in) :: num
  		real(dl), intent(out) :: rho,drhodx,drhody,drhodz
  		type(pixelinfo), intent(in) :: pixel
  		! LOCAL VARIABLES
  		integer, allocatable :: index_array(:)
  		real(dl), allocatable :: distance_array(:), xyz_mass_array(:,:), tmp(:), smda(:)
  		integer, allocatable :: tmpindex(:), smlablist(:)
  		integer :: i,j,n,nowindex
  		real(dl) :: h, r0(3), r, mass, dweight
  		
		h = pixel%maxdist / 2.0
		rho = 0; drhodx=0;drhody=0;drhodz=0;
		do i = 1, num
			r = pixel%xyzrlist(4,i)
			mass = gb_mass_list(pixel%indexlist(i))
			rho = rho + mass*w_kernel(r, h)
			dweight = der_w_kernel(r,h)
			drhodx = drhodx + mass*(x-pixel%xyzrlist(1,i)) / r * dweight
			drhody = drhody + mass*(y-pixel%xyzrlist(2,i)) / r * dweight
			drhodz = drhodz + mass*(z-pixel%xyzrlist(3,i)) / r * dweight
		enddo	
	end subroutine nb_listinput	


  !------------------------------------------
  ! estimating rho and gradient rho based on
  !  cubic spline kernel
  !------------------------------------------
  	subroutine nb_listoutput(x,y,z,num,rho,drhodx,drhody,drhodz,max_dist,pixel)
		! DUMMY ARGUMENTS
  		real(dl), intent(in) :: x,y,z
  		integer, intent(in) :: num
  		real(dl), intent(out) :: rho,drhodx,drhody,drhodz,max_dist
  		type(pixelinfo), intent(out) :: pixel
  		! LOCAL VARIABLES
  		integer, allocatable :: index_array(:)
  		real(dl), allocatable :: distance_array(:), xyz_mass_array(:,:), tmp(:), smda(:)
  		integer, allocatable :: tmpindex(:), smlablist(:)
  		integer :: i,j,n,nowindex
  		real(dl) :: h, r0(3), r, mass, dweight

		if(.not.allocated(pixel%indexlist)) allocate(pixel%indexlist(num))
		if(.not.allocated(pixel%xyzrlist)) allocate(pixel%xyzrlist(4,num))

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
			!!! Save results to pixel
			pixel%indexlist(i) = nowindex
 		enddo

  		max_dist = maxval(smda)
  		pixel%maxdist = max_dist
  		h = max_dist / 2.0
		rho = 0; drhodx=0;drhody=0;drhodz=0;
		do i = 1, num
			r = distance_array(smlablist(i))
			mass = xyz_mass_array(4,i)
			rho = rho + mass*w_kernel(r, h)
			dweight = der_w_kernel(r,h)
			drhodx = drhodx + mass*(x-xyz_mass_array(1,i)) / r * dweight
			drhody = drhody + mass*(y-xyz_mass_array(2,i)) / r * dweight
			drhodz = drhodz + mass*(z-xyz_mass_array(3,i)) / r * dweight
			!!! Save results to pixel
			pixel%xyzrlist(1:3,i) = xyz_mass_array(1:3,i)
			pixel%xyzrlist(4,i) = r
		enddo
	end subroutine nb_listoutput


  !------------------------------------------
  ! estimating rho and gradient rho based on
  !  cubic spline kernel
  !------------------------------------------
  	subroutine nb_list(x,y,z,num,rho,drhodx,drhody,drhodz,max_dist,pixel,acoutput)
		! DUMMY ARGUMENTS
  		real(dl), intent(in) :: x,y,z
  		integer, intent(in) :: num
  		real(dl), intent(out) :: rho,drhodx,drhody,drhodz,max_dist
  		type(pixelinfo), optional :: pixel
  		logical, optional :: acoutput
  		! LOCAL VARIABLES
  		integer, allocatable :: index_array(:)
  		real(dl), allocatable :: distance_array(:), xyz_mass_array(:,:), tmp(:), smda(:)
  		integer, allocatable :: tmpindex(:), smlablist(:)
  		integer :: i,j,n,nowindex
  		real(dl) :: h, r0(3), r, mass, dweight
  		logical :: havepixel = .false.
!  		logical, optional :: hh !lxd

		! check length of the given index, xyz_distance arrays  		
  		if(present(pixel) .or. present(acoutput)) then
  			if(.not. present(pixel) .or. .not. present(acoutput)) then
  				print *, 'ERROR! pixel, actionisoutput shall be given together!'
  				stop
  			endif
  			havepixel = .true.
  		endif
  		
  		if(havepixel) then
  			if(.not.acoutput) then !use as iput!
!  				if(size(pixel%xyzrarray,2) .ne. num .or. size(pixel%index_array).ne.num) then
 ! 					print *, 'ERROR! Len of xyz_distance_array must be ', num
  !					stop
  !				endif
  		
  				max_dist = pixel%maxdist !gv_xyzrarray(4,num)
  				h = max_dist / 2.0
				rho = 0; drhodx=0;drhody=0;drhodz=0;
				do i = 1, num
					r = pixel%xyzrlist(4,i)
					mass = gb_mass_list(pixel%indexlist(i))
					rho = rho + mass*w_kernel(r, h)
					dweight = der_w_kernel(r,h)
					drhodx = drhodx + mass*(x-pixel%xyzrlist(1,i)) / r * dweight
					drhody = drhody + mass*(y-pixel%xyzrlist(2,i)) / r * dweight
					drhodz = drhodz + mass*(z-pixel%xyzrlist(3,i)) / r * dweight
!					print *, 'acinput: i, xyzr, index = ', pixel%xyzrlist(1:4,i), pixel%indexlist(i)
!					print *, 'acinput: i, r, mass, rho, drhodx, drhody, drhodz = ', i, r, mass, rho, drhodx, drhody, drhodz
!				if(present(hh).and.hh) then
!					print *, '   i,r,mass = ', i, r, mass
!				endif
				enddo	
!				print *, 'acinput: ', rho, drhodx, drhody, drhodz
  				return
  			else	
  				if(.not.allocated(pixel%indexlist)) allocate(pixel%indexlist(num))
  				if(.not.allocated(pixel%xyzrlist)) allocate(pixel%xyzrlist(4,num))
!  				if(size(gv_index_array).ne.num.or.size(gv_xyzrarray,1).ne.4.or.size(gv_xyzrarray,2).ne.num) then
! 					print*, 'ERROR (nb_list)! Check the size of xyzrarray, indexarray!'
!  					stop
! 				endif
	  		endif	
  		endif
		
  		call nb_select(x,y,z,num,selected_list=index_array)
  		n=size(index_array)

		! Use bubble sort rather than quick sort (only find out first num smallest elements)
		!  But quick sort may be faster when num is large
		
  		allocate(tmpindex(n), distance_array(n))
		r0(1)=x; r0(2)=y; r0(3)=z;
		do i = 1, n
			distance_array(i) = distance(gb_xyz_list(1:3,index_array(i)),r0)
			tmpindex(i) = index_array(i)
		enddo

		allocate(smlablist(num),smda(num))
!		call ltlablist1(A=distance_array,nA=n,numsm=num,numtolin=0, &
!				smlablist=smlablist,splist=smda)!,Aminout,Amaxout,markin)
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
			rho = rho + mass*w_kernel(r, h)
			dweight = der_w_kernel(r,h)
			drhodx = drhodx + mass*(x-xyz_mass_array(1,i)) / r * dweight
			drhody = drhody + mass*(y-xyz_mass_array(2,i)) / r * dweight
			drhodz = drhodz + mass*(z-xyz_mass_array(3,i)) / r * dweight
!			print *, 'acoutput: i, r, mass, rho, drhodx, drhody, drhodz = ', i, r, mass, rho, drhodx, drhody, drhodz
!			if(present(hh).and.hh) then
!				print *, '   i,r,mass = ', i, r, mass
!			endif
		enddo
!		print *, 'acoutput: ', rho, drhodx, drhody, drhodz
		if(havepixel) then
			if(acoutput) then
				do i = 1, num
					pixel%xyzrlist(1:3,i) = xyz_mass_array(1:3,i)
					pixel%xyzrlist(4,i) = distance_array(smlablist(i))
					pixel%indexlist(i) = tmpindex(smlablist(i))
					!print *, 'acoutput: i, xyzr, index = ', i, pixel%xyzrlist(1:4,i), pixel%indexlist(i)
				enddo
				pixel%maxdist = max_dist
			endif
		endif
!		if(present(hh).and.hh) then!lxd
!			print *, 'rho, drhodx, drhody, drhodz = ', rho, drhodx, drhody, drhodz
!		endif
	end subroutine nb_list
	
!	subroutine grid_rho_drho_list(RSD, AP, num, num_in_x, gv_print_info,v_check_boundary, gv_cb_adjust_ratio, &
!		pos_list, distance_list, rho_list,  drho_list, boundary_rmin, boundary_rmax, max_dist_list)
!		integer :: RSD, AP, num, num_in_x
!		logical,  optional :: gv_print_info, gv_check_boundary
!		real(dl), optional :: gv_cb_adjust_ratio
!		logical  :: print_info, check_boundary 
!		real(dl) :: cb_adjust_ratio, boundary_rmin, boundary_rmax
!		real(dl), allocatable :: pos_list(:,:), distance_list(:), rho_list(:),  drho_list(:,:), max_dist_list(:)
!		integer, allocatable :: index_list(:)
!		integer :: ix, iy, iz, n1, n2

!		if(present(gv_print_info)) then
!			print_info = gv_print_info
!			else
!			print_info = .false.
!		endif
		
!		if(present(gv_check_boundary)) then
!			check_boundary = gv_check_boundary
!			else
!			check_boundary = .false.
!		endif

!		if(present(gv_cb_adjust_ratio)) then
!			cb_adjust_ratio = gv_cb_adjust_ratio
!			else
!			cb_adjust_ratio = 1.0_dl
!		endif

!		call do_cell_initialize(RSD, AP, num_in_x, print_info)
!	end subroutine grid_rho_drho_list
		
end module ap_smooth	
	



