


!####################################
!This module does smooth
!####################################
module ap_grad_fields
use ap_smooth
use ap_mu_tools
	implicit none
	
	logical, public :: gb_dropval, gb_dropdval
	real(dl), public :: gb_dropvalratio(2), gb_dropdvalratio(2)
	
contains	

  !------------------------------------------
  ! estimating val/dval by using given list
  !------------------------------------------
	subroutine get_val_dval_list(num, print_info, check_boundary, boundary_rmin, boundary_rmax, cb_adjust_ratio, use_num_density, &
		pos_list, distance_list, val_list,  dval_list,  max_dist_list, indexlist, xyzrlist)
		logical  :: print_info, check_boundary, use_num_density, hh
		real(dl) :: cb_adjust_ratio, boundary_rmin, boundary_rmax, rmin, rmax,x,y,z,r,val,dvalx,dvaly,dvalz,max_dist
		real(dl), allocatable, intent(inout) :: pos_list(:,:)
		real(dl), allocatable :: distance_list(:), val_list(:),  dval_list(:,:), max_dist_list(:), tmp(:,:), tmpxyzrseg(:,:)
		integer, intent(in), optional :: indexlist(:,:)
		integer, allocatable :: tmpindexseg(:)
		real(dl), intent(in), optional :: xyzrlist(:,:,:)
		integer :: num, i, j, ix, iy, iz, n, n1, n2

		if(present(indexlist) .or. present(xyzrlist)) then
			if(.not.present(indexlist) .or. .not.present(xyzrlist)) then
				print *, '  ERROR (get_val_dval_list)! indexlist and xyzrlist must appear together!'
				stop
			endif
			allocate(tmpindexseg(num),tmpxyzrseg(4,num))
		endif

		if(size(pos_list,1).ne.3) then
			print *, '  ERROR (get_val_dval_list)! dim not 3!'
			stop
		endif
		n=size(pos_list,2)
		if(print_info) then
			print *, 'Calculating val/dval list ...'
		endif
		if(check_boundary.and.print_info) then
			write(*,'(2x,A,f10.5,2x,f10.5,A)') ' Chenck boundary using rmin, rmax = ', boundary_rmin, boundary_rmax, '  ...'
			allocate(tmp(9,n))
		else
			if(print_info)	print *, '  Does not check boundary...'
			allocate(distance_list(n),val_list(n),dval_list(3,n),max_dist_list(n))
		endif
		n2=0
		if(use_num_density) then
			mass_list = 1.0_dl
		endif
		do i = 1, n
			x=pos_list(1,i); y=pos_list(2,i); z=pos_list(3,i)
			r = sqrt(x*x+y*y+z*z)
			if(check_boundary) then
				if(r<boundary_rmin.or.r>boundary_rmax) cycle
			endif
			if(present(indexlist)) then
				call nb_list(x,y,z,num,val,dvalx,dvaly,dvalz,max_dist,tmpindexseg,tmpxyzrseg,actionisoutput=.false.)
			else
				call nb_list(x,y,z,num,val,dvalx,dvaly,dvalz,max_dist)
			endif
!			hh = .false.!lxd
!			if(i<2) hh = .true.
!			call nb_list(x,y,z,num,val,dvalx,dvaly,dvalz,max_dist,hh)
!			if(i<2) then
!				print *, 'i,r,val,dvalx,dvaly,dvalz,max_dist = ', i, r, val, dvalx, dvaly, dvalz, max_dist, mass_list(i)
!				print *
!			endif!!!lxd
			!whether check boundary
			if(check_boundary) then
				if(has_boundary_effect(x,y,z,max_dist,boundary_rmin,boundary_rmax,cb_adjust_ratio)) then
					cycle
				else
					n2=n2+1	
					tmp(1:9,n2) = (/x,y,z,r,val,dvalx,dvaly,dvalz,max_dist/)
				endif
			else
				distance_list(i) = r
				val_list(i) = val
				dval_list(1:3,i) = (/dvalx,dvaly,dvalz/)
				max_dist_list(i) = max_dist
			endif
				
		enddo
		! shall save data from tmp to arrays if check_boudnary
		if(check_boundary) then	
			deallocate(pos_list)
			allocate(pos_list(3,n2),distance_list(n2),val_list(n2),dval_list(3,n2),max_dist_list(n2))
			do i = 1, n2
				pos_list(1:3,i) = tmp(1:3,i)
				distance_list(i) = tmp(4,i)
				val_list(i) = tmp(5,i)
				dval_list(1:3,i) = tmp(6:8,i)
				max_dist_list(i) = tmp(9,i)
			enddo
			if(print_info) then
				print *, '    tot within boundary / tot without boundary effect = ', n, n2
				print *, '    ratio of skipped = ', (n-n2)/(n+0.0)
			endif
		endif
	end subroutine get_val_dval_list
		
  !------------------------------------------
  ! estimating rho/drho at grid points
  !------------------------------------------
	subroutine grid_rho_drho_list(RSD, AP, num, num_in_x, print_info, check_boundary, cb_adjust_ratio, use_num_density, &
		pos_list, distance_list, rho_list,  drho_list, boundary_rmin, boundary_rmax, max_dist_list, &
		indexlist, xyzrlist, mark_mat)
		integer :: RSD, AP, num, num_in_x
		logical  :: print_info, check_boundary, use_num_density
		real(dl) :: cb_adjust_ratio, boundary_rmin, boundary_rmax, rmin, rmax,x,y,z,r,rho,drhox,drhoy,drhoz,max_dist
		real(dl), allocatable :: pos_list(:,:), distance_list(:), rho_list(:),  drho_list(:,:), max_dist_list(:),&
			tmp(:,:), tmpxyzrlist(:,:,:), tmpxyzrseg(:,:)
		integer, allocatable :: tmpindexlist(:,:), tmpindexseg(:)
		integer, allocatable, optional, intent(out) :: indexlist(:,:)
		real(dl), allocatable, optional, intent(out) :: xyzrlist(:,:,:)
		integer :: i,j, ix, iy, iz, n, n1, n2
		integer, allocatable, optional :: mark_mat(:,:,:)
!
		if(present(indexlist) .or. present(xyzrlist)) then
			if(.not.present(indexlist) .or. .not.present(xyzrlist)) then
				print*, 'ERROR (grid_rho_drho_list)! indexlist and xyzrlist must appear together!'
				stop
			endif
		endif
		
		call do_cell_initialize(RSD, AP, num_in_x, print_info)
		
		if(use_num_density) mass_list = 1.0_dl
		
		rmin = minval(r_list)
		rmax = maxval(r_list)
		if(check_boundary) then
			boundary_rmin = rmin 
			boundary_rmax = rmax
			if(print_info) &
				write(*,'(2x,A,f10.5,2x,f10.5,A)') ' Chenck boundary using rmin, rmax = ', boundary_rmin, boundary_rmax, ' ...'
		endif
		if(print_info) &
			print *, '  Estimating rho_list and its gradient...'

		n = n_cellx*n_celly*n_cellz
		
		allocate(tmp(9,n))
		
		if(present(indexlist)) then 
			allocate(tmpindexlist(num,n),tmpxyzrlist(4,num,n),tmpxyzrseg(4,num),tmpindexseg(num))
		endif

		if(present(mark_mat)) then
			allocate(mark_mat(n_cellx,n_celly,n_cellz))
			mark_mat = 0
		endif
		
		n1=0; n2=0;
		do ix = 1, n_cellx
		do iy = 1, n_celly
		do iz = 1, n_cellz
			i=i+1;
			call cell_pos(ix,iy,iz,x,y,z)
			r = sqrt(x*x+y*y+z*z)
			if(check_boundary) then
				if(r<boundary_rmin.or.r>boundary_rmax) cycle
			endif
			n1=n1+1
			if(present(indexlist)) then
				call nb_list(x,y,z,num,rho,drhox,drhoy,drhoz,max_dist,tmpindexseg,tmpxyzrseg, actionisoutput=.true.)
			else
				call nb_list(x,y,z,num,rho,drhox,drhoy,drhoz,max_dist)
				!call nb_list(x,y,z,num,rho,drhox,drhoy,drhoz,max_dist,tmpindexseg,tmpxyzrseg, actionisoutput=.true.)
			endif
			if(check_boundary) then
				if(has_boundary_effect(x,y,z,max_dist,boundary_rmin,boundary_rmax,cb_adjust_ratio)) cycle
			endif	
			n2=n2+1	
			if(present(mark_mat)) mark_mat(ix,iy,iz) = n2
			tmp(1:9,n2) = (/x,y,z,r,rho,drhox,drhoy,drhoz,max_dist/)
			if(present(indexlist)) then 
				tmpindexlist(1:num,n2) =  tmpindexseg(1:num)
				tmpxyzrlist(1,1:num,n2) = tmpxyzrseg(1,1:num)
				tmpxyzrlist(2,1:num,n2) = tmpxyzrseg(2,1:num)
				tmpxyzrlist(3,1:num,n2) = tmpxyzrseg(3,1:num)
				tmpxyzrlist(4,1:num,n2) = tmpxyzrseg(4,1:num)
			endif
		enddo
		enddo
		enddo
		allocate(pos_list(3,n2),distance_list(n2),rho_list(n2),drho_list(3,n2),max_dist_list(n2))
		if(present(indexlist)) allocate(indexlist(num,n2),xyzrlist(4,1:num,n2))
		do i = 1, n2
			pos_list(1:3,i) = tmp(1:3,i)
			distance_list(i) = tmp(4,i)
			rho_list(i) = tmp(5,i)
			drho_list(1:3,i) = tmp(6:8,i)
			max_dist_list(i) = tmp(9,i)
			if(present(indexlist)) then
				indexlist(1:num,n2) = tmpindexlist(1:num,n2)
				do j = 1, 4
					xyzrlist(j,1:num,n2) = tmpxyzrlist(j,1:num,n2)
				enddo
			endif
		enddo
		
		if(present(indexlist))	deallocate(tmpindexlist,tmpxyzrlist)
		
		if(print_info) then
			print *, '  tot within boundary / tot without boundary effect = ', n1, n2
			print *, '  ratio of skipped = ', (n1-n2)/(n1+0.0)
		endif
	end subroutine grid_rho_drho_list
	
!get_rho_drho_list(num,print_info, check_boundary, boundary_rmin, boundary_rmax, cb_adjust_ratio, &
!		pos_list, distance_list, rho_list,  drho_list,  max_dist_list)

  !------------------------------------------
  ! estimating gradient of fields of 
  !  rho/delta/normed delta at grid points
  !------------------------------------------
	subroutine grid_rho_delta_list(RSD, AP, num, num_in_x, print_info, check_boundary, cb_adjust_ratio, remov_dist_ratio, &
		use_num_density, nbins_rhoav, use_intpl_rho, &
		!pos_list, distance_list, rho_list,  drho_list, boundary_rmin, boundary_rmax, max_dist_list, &
		!delta_list, ddelta_list, normed_delta_list, dnormed_delta_list, &
		drho_mu_data, ddelta_mu_data, dnormed_delta_mu_data)
		integer :: RSD, AP, num, num_in_x, nbins_rhoav
		logical  :: print_info, check_boundary, use_num_density, use_intpl_rho
		real(dl) :: cb_adjust_ratio, remov_dist_ratio
		real(dl) :: boundary_rmin, boundary_rmax, rmin, rmax,remov_dist
		real(dl) :: r,delta_r,rho_av,dvar_sqrt,rho1,rho2,tmpx,tmpswitch(10)
		real(dl), allocatable :: pos_list(:,:), distance_list(:), rho_list(:), drho_list(:,:), max_dist_list(:), &
			rho_av_list(:), rho_er_list(:), binned_r_list(:), delta_av_list(:), delta_er_list(:), delta_var_list(:), tmp(:,:), reflist(:)
		real(dl), allocatable :: delta_list(:), ddelta_list(:,:), normed_delta_list(:), dnormed_delta_list(:,:)
		real(dl), allocatable :: drho_mu_data(:), ddelta_mu_data(:), dnormed_delta_mu_data(:)
		integer :: i, j, i_sm, binned_i, n, i1, i2, i3
		integer, allocatable :: indexlist(:,:)
		real(dl), allocatable :: xyzrlist(:,:,:)
		integer :: resmooth_time=0, renorm_time=0, i_drop, i_quan
		logical :: fast_mode = .false. !fast mode still failed...
		logical :: drop_val(3), drop_dval(3)
		real(dl) :: dropvalratio(2,3), dropdvalratio(2,3)

		drop_val(1)  = gb_dropval
		drop_dval(1) = gb_dropdval
		dropvalratio(1:2,1) =  gb_dropvalratio(1:2)
		dropdvalratio(1:2,1) =  gb_dropdvalratio(1:2)   
		

		drop_val(2) = .false.
		drop_dval(2) = .false. 
		dropvalratio(1,2) = 0.0;    dropvalratio(2,2) = 0.0
		dropvalratio(1,3) = 0.0;    dropvalratio(2,3) = 0.0
		drop_val(3) = .false.  
		drop_dval(3) = .false.  
		dropdvalratio(1,2) = 0.0;    dropdvalratio(2,2) = 0.0
		dropdvalratio(1,3) = 0.0;    dropdvalratio(2,3) = 0.0		

		if(print_info) &
			print *, 'Estimating rho, delta, normed_delta...'
		

		call grid_rho_drho_list(RSD, AP, num, num_in_x, print_info, check_boundary, cb_adjust_ratio, use_num_density, &
			pos_list, distance_list, rho_list,  drho_list, boundary_rmin, boundary_rmax, max_dist_list)		

		call drop_pixels(rho_list, distance_list, pos_list, drho_list, &
			drop_val(1), dropvalratio(1:2,1), drop_dval(1), dropdvalratio(1:2,1))

		! mu data		
		call get_mu_from_gradient_list(pos_list, drho_list, drho_mu_data)

!		if(fast_mode) then
!			call grid_rho_drho_list(RSD, AP, num, num_in_x, print_info, check_boundary, cb_adjust_ratio, use_num_density, &
!				pos_list, distance_list, rho_list,  drho_list, boundary_rmin, boundary_rmax, max_dist_list, indexlist, xyzrlist)
!		else
!			call grid_rho_drho_list(RSD, AP, num, num_in_x, print_info, check_boundary, cb_adjust_ratio, use_num_density, &
!				pos_list, distance_list, rho_list,  drho_list, boundary_rmin, boundary_rmax, max_dist_list)
!		endif
		
 		!################################################3
  		!  get delta data (mass_list)
		!################################################3  		
		! get mean value of rho
		n = size(distance_list)
		rmin = minval(r_list)
		rmax = maxval(r_list)
		remov_dist = est_sm_sphe_r(vol_fun(rmin,rmax), num_halo, num)
		boundary_rmin = rmin + remov_dist*remov_dist_ratio
		boundary_rmax = rmax - remov_dist*remov_dist_ratio
		delta_r = (boundary_rmax - boundary_rmin) / (nbins_rhoav + 0.0)
		if(print_info) then
		        print *, '  adopting boundary_remove...'
		        print *, '  remov_dist, remov_dist_ratio = ', real(remov_dist), real(remov_dist_ratio)
        		print *, '  recommend # of particles = ', (vol_fun(rmin, rmax)/(4.0*const_pi/3.0*remov_dist**3.0))
		        print  *, '  after remov, rmin/rmax = ', boundary_rmin, boundary_rmax
		endif
  		
  		do i_sm = 1, 1+resmooth_time
  			call binned_quan(rho_list, distance_list, boundary_rmin, boundary_rmax, nbins_rhoav, &
		  		rho_av_list, rho_er_list, binned_r_list)
	  		do i = 1, size(r_list)
  				r = r_list(i)
  				binned_i = max(min(int((r-boundary_rmin)/delta_r)+1,nbins_rhoav),1)
  				if(use_intpl_rho) then
  					call ilist(binned_i, 1, 1, nbins_rhoav,i1,i3)
  					i2 = i1 + 1
  					rho_av = intpl_vl(r,binned_r_list(i1),rho_av_list(i1),binned_r_list(i2),rho_av_list(i2),binned_r_list(i3),rho_av_list(i3))
  					mass_list(i) = mass_list(i) / rho_av
  				else
  					mass_list(i) = mass_list(i) / rho_av_list(binned_i)
  				endif
  			enddo
  		enddo
  		
  		! delta_list
  		deallocate(distance_list,max_dist_list)
  		
  		if(fast_mode) then
	  		call get_val_dval_list(num, print_info, .false., boundary_rmin, boundary_rmax, cb_adjust_ratio, .false., &
				pos_list, distance_list, delta_list,  ddelta_list,  max_dist_list, indexlist, xyzrlist)
		else
			call get_val_dval_list(num, print_info, .false., boundary_rmin, boundary_rmax, cb_adjust_ratio, .false., &
				pos_list, distance_list, delta_list,  ddelta_list,  max_dist_list)
		endif



		call get_mu_from_gradient_list(pos_list, ddelta_list, ddelta_mu_data)
		!################################################3
		!  get normed_delta data (mass_list)
		!################################################3
  		do i_sm = 1, 1+renorm_time
			call binned_quan(delta_list, distance_list, boundary_rmin, boundary_rmax, nbins_rhoav, &
				delta_av_list, delta_er_list, binned_r_list, delta_var_list)

  			do i = 1, size(r_list)
  				r = r_list(i)
  				binned_i = max(min(int((r-boundary_rmin)/delta_r)+1,nbins_rhoav),1)
  				if(use_intpl_rho) then
  					call ilist(binned_i, 1, 1, nbins_rhoav,i1,i3)
  					i2 = i1 + 1
  					dvar_sqrt = intpl_vl(r,binned_r_list(i1),sqrt(delta_var_list(i1)),binned_r_list(i2),sqrt(delta_var_list(i2)),binned_r_list(i3),sqrt(delta_var_list(i3)))
  					mass_list(i) = mass_list(i) / dvar_sqrt
  				else
  					mass_list(i) = mass_list(i) / sqrt(delta_var_list(binned_i))
  				endif
  			enddo
  		enddo
  		
		!get normed_delta_list
		deallocate(distance_list,max_dist_list)
		

  		if(fast_mode) then
	  		call get_val_dval_list(num, print_info, .false., boundary_rmin, boundary_rmax, cb_adjust_ratio, .false., &
				pos_list, distance_list, normed_delta_list,  dnormed_delta_list,  max_dist_list, indexlist, xyzrlist)
		else
			call get_val_dval_list(num, print_info, .false., boundary_rmin, boundary_rmax, cb_adjust_ratio, .false., &
				pos_list, distance_list, normed_delta_list,  dnormed_delta_list,  max_dist_list)
		endif
		
		
		call get_mu_from_gradient_list(pos_list, dnormed_delta_list, dnormed_delta_mu_data)
	end subroutine grid_rho_delta_list

  !------------------------------------------
  ! Estimating gradient of fields of 
  !  rho/delta/normed delta at grid points
  !  (*finite diff method to estimate drho*)
  !  Have been checked to be not helpful...
  !------------------------------------------
	subroutine grid_rho_delta_list_fd(RSD, AP, num, num_in_x, print_info, check_boundary, cb_adjust_ratio, remov_dist_ratio, &
		use_num_density, nbins_rhoav, use_intpl_rho, &
		!pos_list, distance_list, rho_list,  drho_list, boundary_rmin, boundary_rmax, max_dist_list, &
		!delta_list, ddelta_list, normed_delta_list, dnormed_delta_list, &
		drho_mu_data, ddelta_mu_data, dnormed_delta_mu_data)
		integer :: RSD, AP, num, num_in_x, nbins_rhoav
		logical  :: print_info, check_boundary, use_num_density, use_intpl_rho
		real(dl) :: cb_adjust_ratio, remov_dist_ratio
		real(dl) :: boundary_rmin, boundary_rmax, rmin, rmax,remov_dist
		real(dl) :: r,delta_r,rho_av,dvar_sqrt,rho1,rho2,tmpx,tmpswitch(10)
		real(dl), allocatable :: pos_list(:,:), distance_list(:), rho_list(:), drho_list(:,:), max_dist_list(:), &
			rho_av_list(:), rho_er_list(:), binned_r_list(:), delta_av_list(:), delta_er_list(:), delta_var_list(:), tmp(:,:), reflist(:)
		real(dl), allocatable :: delta_list(:), ddelta_list(:,:), normed_delta_list(:), dnormed_delta_list(:,:)
		real(dl), allocatable :: drho_mu_data(:), ddelta_mu_data(:), dnormed_delta_mu_data(:)
		integer :: i, j, i_sm, binned_i, n, i1, i2, i3, ix,iy,iz, nx,ny,nz, numhasnb, jxleft,jxright,jyleft,jyright,jzleft,jzright
		integer, allocatable :: indexlist(:,:), mark_mat(:,:,:)
		real(dl), allocatable :: xyzrlist(:,:,:)
		integer :: resmooth_time=0, renorm_time=0, i_drop, i_quan
		logical :: fast_mode = .false. !fast mode still failed...
		logical :: drop_val(3), drop_dval(3)
		real(dl) :: dropvalratio(2,3), dropdvalratio(2,3)
		real(dl) :: x, y, z, dx, dy, dz, drhodx, drhody, drhodz

		drop_val(1)  = gb_dropval
		drop_dval(1) = gb_dropdval
		dropvalratio(1:2,1) =  gb_dropvalratio(1:2)
		dropdvalratio(1:2,1) =  gb_dropdvalratio(1:2)   
		

		drop_val(2) = .false.
		drop_dval(2) = .false. 
		dropvalratio(1,2) = 0.0;    dropvalratio(2,2) = 0.0
		dropvalratio(1,3) = 0.0;    dropvalratio(2,3) = 0.0
		drop_val(3) = .false.  
		drop_dval(3) = .false.  
		dropdvalratio(1,2) = 0.0;    dropdvalratio(2,2) = 0.0
		dropdvalratio(1,3) = 0.0;    dropdvalratio(2,3) = 0.0		

		if(print_info) &
			print *, 'Estimating rho, delta, normed_delta...'
		

		call grid_rho_drho_list(RSD, AP, num, num_in_x, print_info, check_boundary, cb_adjust_ratio, use_num_density, &
			pos_list, distance_list, rho_list,  drho_list, boundary_rmin, boundary_rmax, max_dist_list, mark_mat=mark_mat)		

		nx = size(mark_mat,1)
		ny = size(mark_mat,2)
		nz = size(mark_mat,3)

		numhasnb = 0
		allocate(tmp(8,size(rho_list)))
		
		do ix=2,nx-1
		do iy=2,ny-1
		do iz=2,nz-1
			j = mark_mat(ix,iy,iz)
			if(j.eq.0) cycle
			jxleft  =mark_mat(ix-1,iy,iz)
			jxright =mark_mat(ix+1,iy,iz)
			jyleft  =mark_mat(ix,iy-1,iz)
			jyright =mark_mat(ix,iy+1,iz)
			jzleft  =mark_mat(ix,iy,iz-1)
			jzright =mark_mat(ix,iy,iz+1)
			if(jxleft.eq.0 .or. jxright.eq.0 .or. &
			   jyleft.eq.0 .or. jyright.eq.0 .or. &
			   jzleft.eq.0 .or. jzright.eq.0) cycle
			numhasnb = numhasnb+1
			if(numhasnb .eq. 1) then
				dx = pos_list(1,jxright) - pos_list(1,j)
				dy = pos_list(2,jyright) - pos_list(2,j)
				dz = pos_list(3,jzright) - pos_list(3,j)
				dx=dx*2; dy=dy*2; dz=dz*2
			endif
			tmp(1,numhasnb) = rho_list(j)
			tmp(2,numhasnb) = distance_list(j)
			tmp(3:5,numhasnb) = pos_list(1:3,j)
			tmp(6,numhasnb) = (rho_list(jxright) - rho_list(jxleft)) / dx
			tmp(7,numhasnb) = (rho_list(jyright) - rho_list(jyleft)) / dy
			tmp(8,numhasnb) = (rho_list(jzright) - rho_list(jzleft)) / dz
		enddo
		enddo
		enddo

		print *, 'size(rho_list),numhasnb,ratio=', size(rho_list),numhasnb,numhasnb/(size(rho_list)+0.0)
		
		deallocate(rho_list, distance_list, pos_list, drho_list)
		allocate(rho_list(numhasnb), distance_list(numhasnb), pos_list(3,numhasnb), drho_list(3,numhasnb))
		do i = 1, numhasnb
			distance_list(i) = tmp(2,i)
			rho_list(i) = tmp(1,i)
			pos_list(1:3,i) = tmp(3:5,i)
			drho_list(1:3,i) = tmp(6:8,i)
		enddo


		call drop_pixels(rho_list, distance_list, pos_list, drho_list, &
			drop_val(1), dropvalratio(1:2,1), drop_dval(1), dropdvalratio(1:2,1))

		! mu data		
		call get_mu_from_gradient_list(pos_list, drho_list, drho_mu_data)

!		if(fast_mode) then
!			call grid_rho_drho_list(RSD, AP, num, num_in_x, print_info, check_boundary, cb_adjust_ratio, use_num_density, &
!				pos_list, distance_list, rho_list,  drho_list, boundary_rmin, boundary_rmax, max_dist_list, indexlist, xyzrlist)
!		else
!			call grid_rho_drho_list(RSD, AP, num, num_in_x, print_info, check_boundary, cb_adjust_ratio, use_num_density, &
!				pos_list, distance_list, rho_list,  drho_list, boundary_rmin, boundary_rmax, max_dist_list)
!		endif
		
 		!################################################3
  		!  get delta data (mass_list)
		!################################################3  		
		! get mean value of rho
		n = size(distance_list)
		rmin = minval(r_list)
		rmax = maxval(r_list)
		remov_dist = est_sm_sphe_r(vol_fun(rmin,rmax), num_halo, num)
		boundary_rmin = rmin + remov_dist*remov_dist_ratio
		boundary_rmax = rmax - remov_dist*remov_dist_ratio
		delta_r = (boundary_rmax - boundary_rmin) / (nbins_rhoav + 0.0)
		if(print_info) then
		        print *, '  adopting boundary_remove...'
		        print *, '  remov_dist, remov_dist_ratio = ', real(remov_dist), real(remov_dist_ratio)
        		print *, '  recommend # of particles = ', (vol_fun(rmin, rmax)/(4.0*const_pi/3.0*remov_dist**3.0))
		        print  *, '  after remov, rmin/rmax = ', boundary_rmin, boundary_rmax
		endif
  		
  		do i_sm = 1, 1+resmooth_time
  			call binned_quan(rho_list, distance_list, boundary_rmin, boundary_rmax, nbins_rhoav, &
		  		rho_av_list, rho_er_list, binned_r_list)
	  		do i = 1, size(r_list)
  				r = r_list(i)
  				binned_i = max(min(int((r-boundary_rmin)/delta_r)+1,nbins_rhoav),1)
  				if(use_intpl_rho) then
  					call ilist(binned_i, 1, 1, nbins_rhoav,i1,i3)
  					i2 = i1 + 1
  					rho_av = intpl_vl(r,binned_r_list(i1),rho_av_list(i1),binned_r_list(i2),rho_av_list(i2),binned_r_list(i3),rho_av_list(i3))
  					mass_list(i) = mass_list(i) / rho_av
  				else
  					mass_list(i) = mass_list(i) / rho_av_list(binned_i)
  				endif
  			enddo
  		enddo
  		
  		! delta_list
  		deallocate(distance_list,max_dist_list)
  		
  		if(fast_mode) then
	  		call get_val_dval_list(num, print_info, .false., boundary_rmin, boundary_rmax, cb_adjust_ratio, .false., &
				pos_list, distance_list, delta_list,  ddelta_list,  max_dist_list, indexlist, xyzrlist)
		else
			call get_val_dval_list(num, print_info, .false., boundary_rmin, boundary_rmax, cb_adjust_ratio, .false., &
				pos_list, distance_list, delta_list,  ddelta_list,  max_dist_list)
		endif



		call get_mu_from_gradient_list(pos_list, ddelta_list, ddelta_mu_data)
		!################################################3
		!  get normed_delta data (mass_list)
		!################################################3
  		do i_sm = 1, 1+renorm_time
			call binned_quan(delta_list, distance_list, boundary_rmin, boundary_rmax, nbins_rhoav, &
				delta_av_list, delta_er_list, binned_r_list, delta_var_list)

  			do i = 1, size(r_list)
  				r = r_list(i)
  				binned_i = max(min(int((r-boundary_rmin)/delta_r)+1,nbins_rhoav),1)
  				if(use_intpl_rho) then
  					call ilist(binned_i, 1, 1, nbins_rhoav,i1,i3)
  					i2 = i1 + 1
  					dvar_sqrt = intpl_vl(r,binned_r_list(i1),sqrt(delta_var_list(i1)),binned_r_list(i2),sqrt(delta_var_list(i2)),binned_r_list(i3),sqrt(delta_var_list(i3)))
  					mass_list(i) = mass_list(i) / dvar_sqrt
  				else
  					mass_list(i) = mass_list(i) / sqrt(delta_var_list(binned_i))
  				endif
  			enddo
  		enddo
  		
		!get normed_delta_list
		deallocate(distance_list,max_dist_list)
		

  		if(fast_mode) then
	  		call get_val_dval_list(num, print_info, .false., boundary_rmin, boundary_rmax, cb_adjust_ratio, .false., &
				pos_list, distance_list, normed_delta_list,  dnormed_delta_list,  max_dist_list, indexlist, xyzrlist)
		else
			call get_val_dval_list(num, print_info, .false., boundary_rmin, boundary_rmax, cb_adjust_ratio, .false., &
				pos_list, distance_list, normed_delta_list,  dnormed_delta_list,  max_dist_list)
		endif
		
		
		call get_mu_from_gradient_list(pos_list, dnormed_delta_list, dnormed_delta_mu_data)
	end subroutine grid_rho_delta_list_fd



  !------------------------------------------
  ! dropping pixels in lists
  !------------------------------------------
  	subroutine drop_pixels(val_list, distance_list, pos_list, dval_list, dropval, dropvalratio, dropdval, dropdvalratio)
  		real(dl), allocatable :: val_list(:), distance_list(:), pos_list(:,:), dval_list(:,:), reflist(:), tmp(:,:)
  		real(dl) :: tmpx, dropvalratio(2), dropdvalratio(2)
  		integer, allocatable :: indexarray(:)
  		integer :: i, j, i_drop, n, i1, i2, tmpi
  		logical :: dropval, dropdval

		n = size(val_list)
		if(size(distance_list).ne.n .or. size(pos_list,2).ne.n .or. size(dval_list,2).ne. n) then
			print *, 'ERROR (drop_pixels)! Check length of arrays: ', n,size(distance_list),size(pos_list,2),size(dval_list,2)
			stop
		endif
		do i_drop = 1, 2
			if(i_drop.eq.1 .and. .not.dropval) cycle
			if(i_drop.eq.2 .and. .not.dropdval) cycle
			
			if(i_drop .eq. 1) then
				i1 = n*dropvalratio(1); i2 = n*(1-dropvalratio(2));
			else
				i1 = n*dropdvalratio(1); i2 = n*(1-dropdvalratio(2));
			endif
			
			i1 = max(min(i1,n),1); i2 = max(min(i2,n),1);
			if(i2 .le. i1) then
				print *, 'ERROR (drop_pixels)! Dropping index: i1, i2 = ', i1, i2
				stop
			endif
			
			! in some cases no nessicity to do sort ...
			if(i1<=1 .and. i2 >= n) then
!				print *, 'No sort!!!: i_drop=', i_drop
				cycle
			endif
			
			allocate(reflist(n),indexarray(n),tmp(8,n))
			do i = 1, n
				tmp(1,i)=val_list(i); tmp(2,i)=distance_list(i)
				tmp(3:5,i)=pos_list(1:3,i); tmp(6:8,i)=dval_list(1:3,i)	
				indexarray(i) = i
				if(i_drop .eq. 1) then
					reflist(i) = val_list(i)
				elseif(i_drop .eq. 2) then
					reflist(i) = dval_list(1,i)**2.0+dval_list(2,i)**2.0+dval_list(3,i)**2.0
				endif
			enddo
			
			! Do it by using quick sort !!!
			call Qsort2(reflist, indexarray, n)
			
			!do i = 1, n
			!	do j = n, i+1, -1
			!		if(reflist(j-1) > reflist(j)) then
			!			tmpx=reflist(j-1); reflist(j-1)=reflist(j); reflist(j)=tmpx
  			!			tmpi=indexarray(j-1); indexarray(j-1)=indexarray(j); indexarray(j)=tmpi
			!		endif
			!	enddo
			!enddo
			
			! drop some high/low density pixels ...
			n = i2-i1+1
			deallocate(val_list,dval_list,pos_list,distance_list)
			allocate(val_list(n),dval_list(3,n),pos_list(3,n),distance_list(n))
			j = 1
			do i = i1, i2
				val_list(j) = tmp(1,indexarray(i))
				distance_list(j) = tmp(2,indexarray(i))
				pos_list(1:3,j) = tmp(3:5,indexarray(i))
				dval_list(1:3,j) = tmp(6:8,indexarray(i))
				j = j + 1
			enddo
			deallocate(tmp,reflist,indexarray)
		enddo
	end subroutine drop_pixels


  !------------------------------------------
  ! est recommended # of points
  !------------------------------------------
  	real(dl) function est_num_in_x(smnum)
  		integer :: smnum, i
  		real(dl) :: rmin, rmax, x1,x2,y1,y2,z1,z2
  		rmin = halo_info(1)%r; rmax = rmin
  		x1 = halo_info(1)%x; x2 = x1
  		y1 = halo_info(1)%y; y2 = y1
  		z1 = halo_info(1)%z; z2 = z1
  		do i = 2, num_halo
  			rmin = min(rmin, halo_info(i)%r)
  			rmax = max(rmax, halo_info(i)%r)
  			x1 = min(x1, halo_info(i)%x)
  			x2 = max(x2, halo_info(i)%x)
  			y1 = min(y1, halo_info(i)%y)
  			y2 = max(y2, halo_info(i)%y)
  			z1 = min(z1, halo_info(i)%z)
  			z2 = max(z2, halo_info(i)%z)
  		enddo
		
!		print *, x1,x2,y1,y2,z1,z2, (dble(num_halo)/dble(smnum)), (x2-x1)*(y2-y1)*(z2-z2), vol_fun(rmin,rmax)
		
		est_num_in_x = ( (dble(num_halo)/dble(smnum)) * (x2-x1)*(y2-y1)*(z2-z1) / vol_fun(rmin,rmax) )**(1.0_dl/3.0_dl)
  	
  	end function est_num_in_x

end module ap_grad_fields	
	



