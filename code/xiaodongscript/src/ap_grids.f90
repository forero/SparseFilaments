


!####################################
!This module does smooth
!####################################
module ap_grids
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
		real(dl) :: x, y, z
		real(dl) :: mass = 0.0_dl, rho = 0.0_dl ! if cs%use_num_density = .true., then mass is number, rho is number density
		logical :: has_mass = .false.
		logical :: has_diff(30) = .false.
		real(dl) :: diff(3,30)
	end type

	type(pixel), allocatable :: gbpixels(:,:,:)

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
		real(dl) :: x,y,z,r, rminsq,rmaxsq, px,py,pz, mass, ratio, time1, time2, time3

		if(print_info) then
		        write(*,'(A,i4,i4,A)'), '  Initializing grids of pixels with RSD, AP = ', RSD, AP, '  ...'
		endif
		
		call cpu_time(time1)
		if(do_xyzmassinit) then
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
		
		unit_len = (gbxmax - gbxmin) / dble(rl_num_in_x)

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
		
		! thresholds of x/y/z: pixels lies out of this range has no mass
		rminsq = (gbrmin + unit_len * (sqrt(3.0)/2.0))**2.0
		rmaxsq = (gbrmax - unit_len * (sqrt(3.0)/2.0))**2.0
		
		call grid_clean_up()
		allocate(gbpixels(gb_n_cellx, gb_n_celly, gb_n_cellz))
		
		if(use_num_density) gb_mass_list = 1.0_dl
		
		call cpu_time(time2)
		do i = 1, gb_num_xyz_mass
			if(mod(i,100000) .eq. 0) then
				test_print=.true.
			else
				test_print=.false.
			endif
			x= gb_xyz_list(1,i); y= gb_xyz_list(2,i); z= gb_xyz_list(3,i)
			r=gb_r_list(i)
			mass = gb_mass_list(i)
!			if(x .le. xmin .or. x .ge. xmax .or. y.le.ymin .or. y.ge.ymax .or.&
!			   z.le.zmin .or. z.ge.zmax) then

			ix=max(min(int((x-gbxmin)/gbdeltax)+1,gb_n_cellx),1)
			iy=max(min(int((y-gbymin)/gbdeltay)+1,gb_n_celly),1)
			iz=max(min(int((z-gbzmin)/gbdeltaz)+1,gb_n_cellz),1)

			ix2=int((x-gbxmin-gbdeltax/2.0)/gbdeltax+2)
			ix1=ix2-1
			iy2=int((y-gbymin-gbdeltay/2.0)/gbdeltay+2)
			iy1=iy2-1
			iz2=int((z-gbzmin-gbdeltaz/2.0)/gbdeltaz+2)
			iz1=iz2-1
			
			if(test_print) then
				print*
				print *, 'Halo index: ', i
				print *, '   x,y,z,r,mass=', real(x),real(y),real(z),real(r),real(mass)
				print *, '   ix,iy,iz=', ix,iy,iz
				print *, '   ix2,iy2,iz2=',ix2,iy2,iz2
			endif

			! if any of ix1/iy1/iz1 smaller than 1, assign this halo to this single pixel
			if(ix1.le.0 .or. iy1.le.0 .or. iz1.le.0 .or. &
				ix2.gt.gb_n_cellx .or. iy2.gt.gb_n_celly .or. iz2.gt.gb_n_cellz) then
				gbpixels(ix,iy,iz)%mass = gbpixels(ix,iy,iz)%mass + mass
				if(test_print)&
					print *, 'Cycle because ix/iy/iz too small/too large...'
				cycle
			endif
			
			! if the cell has lowest r touches the boundary, assign this halo to this single pixel
			call cell_pos(ix1,iy1,iz1,px,py,pz)
			if(px*px + py*py + pz*pz .le. rminsq) then
				gbpixels(ix,iy,iz)%mass = gbpixels(ix,iy,iz)%mass + mass
				if(test_print)&
					print *, 'Cycle because ix1/iy1/iz1 touches boundary...'
				cycle
			endif

			! if the cell has largest r touches the boundary, assign this halo to this single pixel
			call cell_pos(ix2,iy2,iz2,px,py,pz)
			if(px*px + py*py + pz*pz .ge. rmaxsq) then
				gbpixels(ix,iy,iz)%mass = gbpixels(ix,iy,iz)%mass + mass
				if(test_print)&
					print *, 'Cycle because ix2/iy2/iz2 touches boundary...'
				cycle
			endif

			! o.w., do CIC 
			do ipx=ix1,ix2
			do ipy=iy1,iy2
			do ipz=iz1,iz2
				call cell_pos(ipx,ipy,ipz,px,py,pz)
				ratio = (gbdeltax-abs(x-px))/gbdeltax * (gbdeltay-abs(y-py))/gbdeltay * (gbdeltaz-(z-pz))/gbdeltaz
				gbpixels(ipx,ipy,ipz)%mass = gbpixels(ipx,ipy,ipz)%mass + mass*ratio
				if(test_print) then
					print *, 'ipx,ipy,ipz=',ipx,ipy,ipz
					print *, 'px,py,pz=',real(px),real(py),real(pz)
					print *, 'ratio=',ratio
				endif
			enddo
			enddo
			enddo
		enddo
		call cpu_time(time3)
		
		print *, 'Time used in init_xyz_mass: ', time2-time1
		print *, 'Time used in init_pixels  : ', time3-time2
! 	We shall check this CIC subroutine carefully...		
	end subroutine init_pixels
		
end module ap_grids
	



