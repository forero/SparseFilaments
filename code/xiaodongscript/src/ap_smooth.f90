


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
	integer :: num_xyz_mass
	real(dl), allocatable :: xyz_list(:,:), mass_list(:), bf_mass_list(:), r_list(:)
	real(dl) :: xmin, xmax, ymin, ymax, zmin, zmax
	real(dl) :: gbrmin, gbrmax, gbtotvol
	real(dl) :: unit_len
	integer :: n_cellx, n_celly, n_cellz
	real(dl) :: deltax, deltay, deltaz
	real(dl) :: cell_vol

	type :: halo_list
		integer :: halo_num = 0
		integer, allocatable :: list(:)
	end type
	
	! array containing which halo in which cell
	type(halo_list), allocatable :: cell_mat(:,:,:)
	
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
		if(allocated(xyz_list)) &
			deallocate(xyz_list)
		if(allocated(mass_list)) &
			deallocate(mass_list)
		if(allocated(bf_mass_list)) &
			deallocate(bf_mass_list)
		if(allocated(cell_mat)) &
			deallocate(cell_mat)
		if(allocated(r_list)) &
			deallocate(r_list)			
	end subroutine

  !------------------------------------------
  ! initialize the xyz_mass array
  !------------------------------------------
	subroutine init_xyz_r_mass_list(RSD, AP, gv_print_info)
		integer :: RSD, AP, i
		real(dl) :: x, y, z, ratio, r
		logical, optional :: gv_print_info
		logical :: print_info

		if(present(gv_print_info)) then
			if(gv_print_info) &
			        write(*,'(A,i4,i4,A)'), '  Initializing grids of cells with RSD, AP = ', RSD, AP, '  ...'
		endif

		call smooth_clean_up()

		allocate(xyz_list(3,num_halo), mass_list(num_halo), bf_mass_list(num_halo), r_list(num_halo))
		num_xyz_mass = num_halo
		do i = 1, num_halo
			x = halo_info(i)%x
			y = halo_info(i)%y
			z = halo_info(i)%z
			r = halo_info(i)%r
			if    (RSD .eq. 0 .and. AP .eq. 0) then
				r_list(i) = r 
			elseif(RSD .eq. 1 .and. AP .eq. 0) then
				r_list(i) = halo_info(i)%r_RSD
			elseif(RSD .eq. 0 .and. AP .gt. 0) then
				r_list(i) = halo_info(i)%r_AP(AP)
			elseif(RSD .eq. 1 .and. AP .gt. 0) then
				r_list(i) = halo_info(i)%r_AP_RSD(AP)
			else
				print *, 'ERROR! No such RSD, AP : ', RSD, AP
				stop
			endif
			ratio = r_list(i) / r 
			xyz_list(1,i) = x * ratio
			xyz_list(2,i) = y * ratio
			xyz_list(3,i) = z * ratio
			mass_list(i) = halo_info(i)%mass
			bf_mass_list(i) = mass_list(i)
		enddo
	end subroutine init_xyz_r_mass_list
	
  !------------------------------------------
  ! read in the xyz_mass file
  !------------------------------------------
  	subroutine read_in_xyz_mass(xyz_mass_file_name)
  		character(len=char_len) :: xyz_mass_file_name
		integer :: i
		call count_line_number(xyz_mass_file_name, num_xyz_mass)
		print*, '   num of halo = ', num_xyz_mass, '...'
	!allocating arrays...
		if(allocated(xyz_list)) &
			deallocate(xyz_list)
		if(allocated(mass_list)) &			
			deallocate(mass_list)
		allocate(xyz_list(3,num_xyz_mass))
		allocate(mass_list(num_xyz_mass))
		open(unit=1,file=xyz_mass_file_name)
		do i = 1, num_xyz_mass
			read(1,*) xyz_list(1:3,i), mass_list(i)
		enddo
	end subroutine read_in_xyz_mass


  !------------------------------------------
  ! get the list of x,y,z at centers of cells
  !------------------------------------------
	subroutine get_cell_pos_list(pos_list)
		real(dl), allocatable :: pos_list(:,:)
		real(dl) :: x,y,z
		integer :: ix, iy, iz, i
		allocate(pos_list(3,n_cellx*n_celly*n_cellz))
		i = 1
		do ix = 1, n_cellx
		do iy = 1, n_celly
		do iz = 1, n_cellz
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
				call init_xyz_r_mass_list(RSD, AP, print_info)
			endif
		else
			call init_xyz_r_mass_list(RSD, AP, print_info)
		endif

		!min/max of x,y,z
		xmin = xyz_list(1,1); xmax = xyz_list(1,1);
		ymin = xyz_list(2,1); ymax = xyz_list(2,1);
		zmin = xyz_list(3,1); zmax = xyz_list(3,1);
		do i = 2, num_xyz_mass
			x= xyz_list(1,i); y= xyz_list(2,i); z= xyz_list(3,i)
			xmin = min(xmin,x); xmax = max(xmax,x)
			ymin = min(ymin,y); ymax = max(ymax,y)
			zmin = min(zmin,z); zmax = max(zmax,z)
		enddo
		!min/max of r
		call find_min_max(r_list, num_halo, gbrmin, gbrmax)
		gbtotvol = vol_fun(gbrmin, gbrmax)
		
		if(print_info) then
			write(*,'(2x,A,f10.5,2x,f10.5)')  '  xmin / xmax = ', xmin, xmax
			write(*,'(2x,A,f10.5,2x,f10.5)')  '  ymin / ymax = ', ymin, ymax
			write(*,'(2x,A,f10.5,2x,f10.5)')  '  zmin / zmax = ', zmin, zmax
		endif
		
	!basic info for cell division	
!		if(use_sd_unit_len) then
!			unit_len = sd_unit_len
!		else
			unit_len = (xmax - xmin) / dble(num_in_x)
!		endif

		n_cellx = int( (xmax-xmin)/unit_len + 0.5)
		n_celly = int( (ymax-ymin)/unit_len + 0.5)
		n_cellz = int( (zmax-zmin)/unit_len + 0.5)

		deltax = (xmax-xmin) / dble(n_cellx)
		deltay = (ymax-ymin) / dble(n_celly)
		deltaz = (zmax-zmin) / dble(n_cellz)

		cell_vol = deltax * deltay * deltaz

		if(print_info) then
			print*,  '   n_cellx, deltax = ', n_cellx, deltax
			print*,  '   n_celly, deltay = ', n_celly, deltay
			print*,  '   n_cellz, deltaz = ', n_cellz, deltaz
			print*,  '   cell_vol = ', cell_vol
		endif

		
	!allocating cell_mat 
		if(allocated(cell_mat)) &
			deallocate(cell_mat)
		allocate(cell_mat(n_cellx,n_celly,n_cellz))	
		allocate(tmp_list(max_incell_halo_num))
				
		do i = 1, num_xyz_mass
			x= xyz_list(1,i); y= xyz_list(2,i); z= xyz_list(3,i)
			ix=max(min(int((x-xmin)/deltax)+1,n_cellx),1)
			iy=max(min(int((y-ymin)/deltay)+1,n_celly),1)
			iz=max(min(int((z-zmin)/deltaz)+1,n_cellz),1)
			
			now_halo_num = cell_mat(ix,iy,iz)%halo_num
			
			if(now_halo_num .eq. max_incell_halo_num) then
				print *, 'ERROR! incell halo num overflows'
				stop
			endif
			
			if(mod(now_halo_num,re_alo_num) .eq. 0) then
				if(now_halo_num .eq. 0) then
					allocate(cell_mat(ix,iy,iz)%list(re_alo_num))
				else
!					print*, 'ix,iy,iz,now_halo_num = ', ix,iy,iz,now_halo_num
					tmp_list(1:now_halo_num) = cell_mat(ix,iy,iz)%list(1:now_halo_num)
					deallocate(cell_mat(ix,iy,iz)%list)
					allocate(cell_mat(ix,iy,iz)%list(now_halo_num+re_alo_num))
					cell_mat(ix,iy,iz)%list(1:now_halo_num) = tmp_list(1:now_halo_num)
				endif
			endif
			
			now_halo_num = now_halo_num + 1
			cell_mat(ix,iy,iz)%halo_num = now_halo_num
			cell_mat(ix,iy,iz)%list(now_halo_num) = i
		enddo
		deallocate(tmp_list)
		if(print_info) then
			print *, '  Cell initialization done.'
		endif
	end subroutine
	
	
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
  		ix = max(min(int((x-xmin)/deltax)+1,n_cellx),1)
  		iy = max(min(int((y-ymin)/deltay)+1,n_celly),1)
  		iz = max(min(int((z-zmin)/deltaz)+1,n_cellz),1)
  		dev = max(abs(abs(xmin + (ix-1)*deltax - x)/deltax - 0.5), &
  			abs(abs(ymin + (iy-1)*deltay - y)/deltay - 0.5), &
  			abs(abs(zmin + (iz-1)*deltaz - z)/deltaz - 0.5) )
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
			call ilist(ix,di,1,n_cellx,i1,i2)
			call ilist(iy,di,1,n_celly,j1,j2)
			call ilist(iz,di,1,n_cellz,k1,k2)
			now_num = 0
			do i = i1, i2
			do j = j1, j2
			do k = k1, k2
				now_num = now_num + cell_mat(i,j,k)%halo_num
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
			do l2 = 1, cell_mat(i,j,k)%halo_num
				selected_list(l1) = cell_mat(i,j,k)%list(l2)
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
  		x = xmin + (ix-0.5)*deltax
  		y = ymin + (iy-0.5)*deltay
  		z = zmin + (iz-0.5)*deltaz
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
			distance_array(i) = distance(xyz_list(1:3,index_array(i)),r0)
			tmpindex(i) = index_array(i)
		enddo

		allocate(smlablist(num),smda(num))
		call ltlablist2(distance_array,n,num,smda,smlablist)

		allocate(xyz_mass_array(4,num))
		
		do i = 1, num
			nowindex = tmpindex(smlablist(i))
			xyz_mass_array(1:3,i) = xyz_list(1:3,nowindex)
			xyz_mass_array(4,i) = mass_list(nowindex)
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
			mass = mass_list(pixel%indexlist(i))
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
			distance_array(i) = distance(xyz_list(1:3,index_array(i)),r0)
			tmpindex(i) = index_array(i)
		enddo

		allocate(smlablist(num),smda(num))
		call ltlablist2(distance_array,n,num,smda,smlablist)

		allocate(xyz_mass_array(4,num))
		
		do i = 1, num
			nowindex = tmpindex(smlablist(i))
			xyz_mass_array(1:3,i) = xyz_list(1:3,nowindex)
			xyz_mass_array(4,i) = mass_list(nowindex)
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
					mass = mass_list(pixel%indexlist(i))
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
			distance_array(i) = distance(xyz_list(1:3,index_array(i)),r0)
			tmpindex(i) = index_array(i)
		enddo

		allocate(smlablist(num),smda(num))
!		call ltlablist1(A=distance_array,nA=n,numsm=num,numtolin=0, &
!				smlablist=smlablist,splist=smda)!,Aminout,Amaxout,markin)
		call ltlablist2(distance_array,n,num,smda,smlablist)

		allocate(xyz_mass_array(4,num))
		
		do i = 1, num
			nowindex = tmpindex(smlablist(i))
			xyz_mass_array(1:3,i) = xyz_list(1:3,nowindex)
			xyz_mass_array(4,i) = mass_list(nowindex)
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
	



