
! Please check dft_ra_ratio in ap_smooth.f90.
!
! Quick find is successful...
!  Now, we can make quick rho_delta...
!   For each pixel, we save its xyzrlist & indexlist, for further usage ...
!   Adjust # of cells to make nb_list more efficient...

!####################################
!This module does smooth
!####################################
module ap_grad_fields
use ap_smooth
use ap_mu_tools
	implicit none
	
	logical, public :: gb_dropmethod(2) ! 1 means ratio; 2 means variance
	logical, public :: gb_dropval(2), gb_dropdval(2)
	real(dl), public :: gb_dropvalratio(2,2), gb_dropdvalratio(2,2)
	logical, public :: gradfieldprintinfo = .true.
	
contains	

  !------------------------------------------
  ! estimating val/dval by using given list
  !------------------------------------------
	subroutine get_val_dval_list(smnum, print_info, pos_list, val_list, dval_list, pilist)
		! Dummy
		integer :: smnum
		logical  :: print_info
		real(dl), allocatable, intent(in) :: pos_list(:,:)
		real(dl), allocatable, intent(out) :: val_list(:), dval_list(:,:)
		type(pixelinfo), optional, intent(in) :: pilist(:)
		! LOCAL
		real(dl) :: x,y,z,r,val,dvalx,dvaly,dvalz,max_dist
		integer :: i, j, ix, iy, iz, n
		real(dl), allocatable :: tmp(:,:), tmpxyzrseg(:,:)
		integer, allocatable :: tmpindexseg(:)

		if(size(pos_list,1).ne.3) then
			print *, '  ERROR (get_val_dval_list)! dim not 3!'
			stop
		endif
		
		n=size(pos_list,2)
		
		if(print_info) then
			print *, 'Calculating val/dval list ...'
		endif

		allocate(val_list(n),dval_list(3,n))

		do i = 1, n
			x=pos_list(1,i); y=pos_list(2,i); z=pos_list(3,i)
			r = sqrt(x*x+y*y+z*z)
			if(present(pilist)) then
				call nb_listinput(x,y,z,smnum,val,dvalx,dvaly,dvalz,pilist(i))
			else
				call nb_list(x,y,z,smnum,val,dvalx,dvaly,dvalz,max_dist)
			endif
			val_list(i) = val
			dval_list(1:3,i) = (/dvalx,dvaly,dvalz/)
		enddo
	end subroutine get_val_dval_list


!----------------------------------------------------
!----------------------------------------------------
!  Further tries (multiple calculation of chisq)
!...

  !------------------------------------------
  ! estimating rho/drho at grid points
  !------------------------------------------
	subroutine grid_rho_drho_list(RSD, AP, cs, changenuminx, pos_list, rho_list, drho_list, boundary_rmin, boundary_rmax, pilist)
		!Dummy
		type(chisq_settings) :: cs
		integer, intent(in) :: RSD, AP 
		real(dl), intent(in) :: changenuminx
		real(dl) :: boundary_rmin, boundary_rmax
		real(dl), allocatable :: max_dist_list(:)
		type(pixelinfo), allocatable, optional :: pilist(:)
		real(dl), allocatable, intent(out) :: pos_list(:,:),  rho_list(:),  drho_list(:,:)

		! Local
		real(dl) :: x,y,z,r,rho,drhox,drhoy,drhoz,max_dist
		real(dl), allocatable :: tmp(:,:), tmpxyzrlist(:,:,:), tmpxyzrseg(:,:)
		integer :: i,j, ix, iy, iz, n, n1, n2
		
		call do_cell_init2(RSD, AP, real(cs%num_in_x)*(1.0_dl+changenuminx), cs%print_info)
		
		if(cs%use_num_density) then
			gb_mass_list = 1.0_dl
		else
			gb_mass_list = gb_bf_mass_list
		endif
		
		if(cs%check_boundary) then
			boundary_rmin = gbrmin 
			boundary_rmax = gbrmax
			if(cs%print_info) &
				write(*,'(2x,A,f10.5,2x,f10.5,A)') ' Chenck boundary using rmin, rmax = ', boundary_rmin, boundary_rmax, ' ...'
		endif
		if(cs%print_info) &
			print *, '  Estimating rho_list and its gradient...'

		n = gb_n_cellx*gb_n_celly*gb_n_cellz
		
		allocate(tmp(7,n))
		
		if(present(pilist)) then
			allocate(pilist(n))
		endif
		
		! n1 counts how many pixels in shell
		! n2 counts how many pixels in shell and has no boundary effect
		n1=0; n2=0;
		do ix = 1, gb_n_cellx
		do iy = 1, gb_n_celly
		do iz = 1, gb_n_cellz
			i=i+1;
			call cell_pos(ix,iy,iz,x,y,z)
			r = sqrt(x*x+y*y+z*z)
			if(cs%check_boundary) then
				if(r<boundary_rmin.or.r>boundary_rmax) cycle
			endif
			n1=n1+1
			if(present(pilist)) then
				call nb_listoutput(x,y,z,cs%smnum,rho,drhox,drhoy,drhoz,max_dist, &
					pilist(n2+1))
			else
				call nb_list0(x,y,z,cs%smnum,rho,drhox,drhoy,drhoz,max_dist)
			endif
			if(cs%check_boundary) then
				if(has_boundary_effect(x,y,z,max_dist,boundary_rmin,boundary_rmax,cs%cb_adjust_ratio)) cycle
			endif	
			n2=n2+1	
			tmp(1:7,n2) = (/x,y,z,rho,drhox,drhoy,drhoz/)
		enddo
		enddo
		enddo
		
!		if(gbtp) then
!			print *, '  # of pixels inside shell:            ', n1
!			print *, '  # of pixels without boundary effect: ', n2
!		endif
		
		allocate(pos_list(3,n2),rho_list(n2),drho_list(3,n2))
		do i = 1, n2
			pos_list(1:3,i) = tmp(1:3,i)
			rho_list(i) = tmp(4,i)
			drho_list(1:3,i) = tmp(5:7,i)
		enddo
		
		if(cs%print_info) then
			print *, '  tot within boundary / tot without boundary effect = ', n1, n2
			print *, '  ratio of skipped = ', (n1-n2)/(n1+0.0)
		endif
	end subroutine grid_rho_drho_list
	

  !------------------------------------------
  ! estimating gradient of fields of 
  !  rho/delta/normed delta at grid points
  !------------------------------------------
	subroutine grid_rho_delta_list(RSD, AP, cs, changenuminx, drho_mu_data, ddelta_mu_data, dnormed_delta_mu_data)
		! Dummy
		integer, intent(in) :: RSD, AP
		type(chisq_settings) :: cs
		real(dl), intent(in) :: changenuminx
		real(dl), allocatable, intent(out) :: drho_mu_data(:), ddelta_mu_data(:), dnormed_delta_mu_data(:)
		! Local
		real(dl) :: boundary_rmin, boundary_rmax, remov_dist
		real(dl) :: r,delta_r,rho_av,dvar_sqrt,rho1,rho2,tmpx,tmpswitch(10)
		real(dl), allocatable :: pos_list(:,:), distance_list(:), rho_list(:), drho_list(:,:), absdrholist(:),  &
			rho_av_list(:), rho_er_list(:), binned_r_list(:)
		real(dl), allocatable :: delta_av_list(:), delta_er_list(:), delta_var_list(:)
		real(dl), allocatable :: delta_list(:), ddelta_list(:,:), normed_delta_list(:), dnormed_delta_list(:,:)
		type(pixelinfo), allocatable :: pilist(:)
		real(dl), allocatable :: tmp1d(:), tmp(:,:), reflist(:), rho_list2(:), drho_list2(:,:), distance_list2(:)
		integer :: i, j, k, i_sm, binned_i, n, n1, n2, n3, i1, i2, i3
		integer, allocatable :: markdrop(:)
		integer :: resmooth_time=0, renorm_time=0
		logical :: fast_mode = .true., print_time = .false.
		real(dl) :: time0, time1, time2, time3, time4, time5, time6, time7

		call cpu_time(time0)
		
		if(cs%print_info) print *, 'Estimating rho, delta, normed_delta: changenuminx = ', real(changenuminx), '...'

		if(fast_mode) then
			call grid_rho_drho_list(RSD, AP, cs, changenuminx, pos_list, rho_list, drho_list, &
				boundary_rmin, boundary_rmax, pilist)		
		else
			call grid_rho_drho_list(RSD, AP, cs, changenuminx, pos_list, rho_list, drho_list, &
				boundary_rmin, boundary_rmax)
		endif
		call cpu_time(time1)		

		if(.not.gb_dropmethod(1) .and. .not.gb_dropmethod(2)) goto 99 ! no drop
		
		n1 = size(rho_list)
		allocate(absdrholist(n1),markdrop(n1),reflist(n1))
		do i = 1, n1
			absdrholist(i)=sqrt(drho_list(1,i)**2.0+drho_list(2,i)**2.0+drho_list(3,i)**2.0)
			markdrop(i) = 0
		enddo
		if(gb_dropmethod(1)) then
			if(gb_dropval(1)) then
				reflist = rho_list
				call mark_drop_pixels(reflist, markdrop, n1, int(gb_dropvalratio(1,1)*n1)-1, int(n1*(1-gb_dropvalratio(2,1)))+1)
			endif
			if(gb_dropdval(1)) then
				reflist = absdrholist
				call mark_drop_pixels(reflist, markdrop, n1, int(gb_dropdvalratio(1,1)*n1)-1, int(n1*(1-gb_dropdvalratio(2,1)))+1)
			endif
			if(gradfieldprintinfo) then
				n2 = 0
				do i = 1, n1
					if(markdrop(i).eq.0) n2 = n2+1
				enddo
				write(*,'(1x,A,i7,f8.5)'), '  left_# / drop_ratio  in drop 1:      ', n2, (n1-n2)/(n1+0.0)
				gradfieldprintinfo = .false.
			endif
		endif
		
		if(gb_dropmethod(2)) then
			call sigma_drop_pixels(rho_list, pos_list,  absdrholist, markdrop, n1, &
				gb_dropval(2), gb_dropvalratio(1:2,2), gb_dropdval(2), gb_dropdvalratio(1:2,2))
		endif

		call drop_pixels(rho_list, pos_list, drho_list, markdrop)
		n3 = size(rho_list)
		if(fast_mode) then
			j = 0
			do i = 1, n1
				if(markdrop(i) .ne. 0) cycle
				j = j + 1
				if(i.eq.j) cycle
				pilist(j)%maxdist = pilist(i)%maxdist
				do k = 1, cs%smnum
					pilist(j)%xyzrlist(1:4,k) = pilist(i)%xyzrlist(1:4,k)
					pilist(j)%indexlist(k) = pilist(i)%indexlist(k)
				enddo
			enddo
			do i = j+1, size(rho_list)
				deallocate(pilist(i)%xyzrlist)
				deallocate(pilist(i)%indexlist)
			enddo
		endif
				
!		if(gbtp) then
!			write(*,'(1x,A,i7,f8.5)'), '  left_# / drop_ratio (tot drop):      ', n3, (n1-n3)/(n1+0.0)
!		endif
		deallocate(reflist, absdrholist, markdrop)

		call cpu_time(time2)		
99		continue

		! mu data		
		call get_mu_from_gradient_list(pos_list, drho_list, drho_mu_data)

 		!################################################3
  		!  get delta data (mass_list)
		!################################################3  		
		! get mean value of rho
		n = size(rho_list)
		
		allocate(distance_list(n))
		do i = 1, n
			distance_list(i) = sqrt(pos_list(1,i)**2.0+pos_list(2,i)**2.0+pos_list(3,i)**2.0)
		enddo
		
!		call find_min_max(r_list,size(r_list),rmin,rmax)
		remov_dist = est_sm_sphe_r(vol_fun(gbrmin,gbrmax), num_halo, cs%smnum)
		boundary_rmin = gbrmin + remov_dist*cs%remov_dist_ratio
		boundary_rmax = gbrmax - remov_dist*cs%remov_dist_ratio
		delta_r = (boundary_rmax - boundary_rmin) / (cs%nbins_rhoav + 0.0)
		if(cs%print_info) then
		        print *, '  adopting boundary_remove...'
		        print *, '  remov_dist, remov_dist_ratio = ', real(remov_dist), real(cs%remov_dist_ratio)
        		print *, '  recommend # of particles = ', (vol_fun(gbrmin, gbrmax)/(4.0*const_pi/3.0*remov_dist**3.0))
		        print  *, '  after remov, rmin/rmax = ', boundary_rmin, boundary_rmax
		endif
  		
  		do i_sm = 1, 1+resmooth_time
  			call binned_quan(rho_list, distance_list, boundary_rmin, boundary_rmax, cs%nbins_rhoav, &
		  		rho_av_list, rho_er_list, binned_r_list)
	  		do i = 1, size(gb_r_list)
  				r = gb_r_list(i)
  				binned_i = max(min(int((r-boundary_rmin)/delta_r)+1,cs%nbins_rhoav),1)
  				if(cs%use_intpl_rho) then
  					call ilist(binned_i, 1, 1, cs%nbins_rhoav,i1,i3)
  					i2 = i1 + 1
  					rho_av = intpl_vl(r,binned_r_list(i1),rho_av_list(i1),binned_r_list(i2),rho_av_list(i2),binned_r_list(i3),rho_av_list(i3))
  					gb_mass_list(i) = gb_mass_list(i) / rho_av
  				else
  					gb_mass_list(i) = gb_mass_list(i) / rho_av_list(binned_i)
  				endif
  			enddo
  		enddo
  		
		call cpu_time(time3)
  		
  		if(fast_mode) then
			call get_val_dval_list(cs%smnum, cs%print_info, pos_list, delta_list, ddelta_list, pilist)
		else
			call get_val_dval_list(cs%smnum, cs%print_info, pos_list, delta_list, ddelta_list)
		endif

		call get_mu_from_gradient_list(pos_list, ddelta_list, ddelta_mu_data)
		call cpu_time(time4)
		
		!################################################3
		!  get normed_delta data (mass_list)
		!################################################3
  		do i_sm = 1, 1+renorm_time
			call binned_quan(delta_list, distance_list, boundary_rmin, boundary_rmax, cs%nbins_rhoav, &
				delta_av_list, delta_er_list, binned_r_list, delta_var_list)

  			do i = 1, size(gb_r_list)
  				r = gb_r_list(i)
  				binned_i = max(min(int((r-boundary_rmin)/delta_r)+1,cs%nbins_rhoav),1)
  				if(cs%use_intpl_rho) then
  					call ilist(binned_i, 1, 1, cs%nbins_rhoav,i1,i3)
  					i2 = i1 + 1
  					dvar_sqrt = intpl_vl(r,binned_r_list(i1),sqrt(delta_var_list(i1)),binned_r_list(i2),&
  						sqrt(delta_var_list(i2)),binned_r_list(i3),sqrt(delta_var_list(i3)))
  					gb_mass_list(i) = gb_mass_list(i) / dvar_sqrt
  				else
  					gb_mass_list(i) = gb_mass_list(i) / sqrt(delta_var_list(binned_i))
  				endif
  			enddo
  		enddo
  		
  		if(fast_mode) then
			call get_val_dval_list(cs%smnum, cs%print_info, pos_list, normed_delta_list, dnormed_delta_list,  pilist)
		else
			call get_val_dval_list(cs%smnum, cs%print_info, pos_list, normed_delta_list, dnormed_delta_list)
		endif
		
		call get_mu_from_gradient_list(pos_list, dnormed_delta_list, dnormed_delta_mu_data)
		call cpu_time(time5)
		deallocate(pos_list,rho_list,drho_list,delta_list,ddelta_list,normed_delta_list,dnormed_delta_list,distance_list)
		call cpu_time(time6)
		if(fast_mode) then
			do i = 1, n3
				deallocate(pilist(i)%indexlist, pilist(i)%xyzrlist)
			enddo
			deallocate(pilist)
		endif
		call cpu_time(time7)

		if(print_time) then
			print *, 'Time used in grid_rho_drho_list: ', time1-time0
			print *, 'Time used in drop: ', time2-time1
			print *, 'Time used in smooth rho: ', time3-time2
			print *, 'Time used in delta list: ', time4-time3
			print *, 'Time used in smooth delta and ndelta list: ', time5-time4
			print *, 'Time used in deallocate many lists: ', time6-time5
			print *, 'Time used in deallocate pilist: ', time7 - time6
		endif
	end subroutine grid_rho_delta_list
	
  !------------------------------------------
  ! estimating gradient of fields of 
  !  rho/delta/normed delta at grid points
  !------------------------------------------
	subroutine gd_mldprho_chi2s(RSD, AP, cs, changenuminx, chisqlist)
		! Dummy
		integer, intent(in) :: RSD, AP
		type(chisq_settings) :: cs
		real(dl), intent(in) :: changenuminx
		real(dl), intent(out) :: chisqlist(:)
		! Local
		real(dl) :: boundary_rmin, boundary_rmax
		integer :: i, j, k, idrop
		integer :: n, n1, n2
		real(dl), allocatable :: drho_mu_data(:)
		real(dl), allocatable :: pos_list0(:,:), rho_list0(:), drho_list0(:,:)
		real(dl), allocatable :: pos_list(:,:), rho_list(:), drho_list(:,:), absdrholist(:), reflist(:)
		integer, allocatable :: markdrop(:)
		logical :: print_time = .false.
		real(dl) :: tmp, time0, time1, time2, time3, time4, time5, time6, time7

		call cpu_time(time0)
		
		if(cs%print_info) print *, 'Estimating rho, delta, normed_delta: changenuminx = ', real(changenuminx), '...'

!		print *, '(gd_mldprho_chi2s) Test A...'
		call grid_rho_drho_list(RSD, AP, cs, changenuminx, pos_list0, rho_list0, drho_list0, &
				boundary_rmin, boundary_rmax)
!		print *, '(gd_mldprho_chi2s) Test B...'

		call cpu_time(time1)		
		
		if(size(chisqlist).lt.cs%numdrop) then
			print *, 'ERROR (gd_mldprho_chi2s)! size of chisqlist less than ', cs%numdrop
			stop
		endif
		
		n1 = size(rho_list0)
!		print *, '(gd_mldprho_chi2s) Test C...'
		allocate(absdrholist(n1),markdrop(n1),reflist(n1))
		do i = 1, n1
			absdrholist(i)=sqrt(drho_list0(1,i)**2.0+drho_list0(2,i)**2.0+drho_list0(3,i)**2.0)
		enddo

!		print *, '(gd_mldprho_chi2s) Test D...'
		do idrop = 1, cs%numdrop
			! cycle if no drop
			if(.not.cs%dropval(idrop) .and. .not.cs%dropdval(idrop)) then
				call get_mu_from_gradient_list(pos_list0, drho_list0, drho_mu_data)
				chisqlist(idrop) = chisq_of_mu_data2(drho_mu_data)
				deallocate(drho_mu_data)
			else
				markdrop = 0
				if(cs%dropval(idrop)) then
					reflist = rho_list0
					call mark_drop_pixels(reflist, markdrop, n1,&
						int(cs%lowdropvalratio(idrop)*n1)-1, int(n1*( 1.0_dl-cs%highdropvalratio(idrop) )) + 1)
				endif
				if(cs%dropdval(idrop)) then
					reflist = absdrholist
					call mark_drop_pixels(reflist, markdrop, n1,&
						int(cs%lowdropdvalratio(idrop)*n1)-1, int((1-cs%highdropdvalratio(idrop))*n1)+1)

				endif

				if(gradfieldprintinfo) then
					n2 = 0
					do i = 1, n1
						if(markdrop(i).eq.0) n2 = n2+1
					enddo
					if(cs%dropval(idrop)) &
						write(*,'(1x,A,f6.2,1x,f6.2)') '  Drop rho ratio: low/high  =     ', &
							cs%lowdropvalratio(idrop), cs%highdropvalratio(idrop)
					if(cs%dropdval(idrop)) &
						write(*,'(1x,A,f6.2,1x,f6.2)') '  Drop drho ratio: low/high =     ', &
							cs%lowdropdvalratio(idrop), cs%highdropdvalratio(idrop)
					write(*,'(1x,A,i7,f8.5)'), '  left_# / drop_ratio  in drop:        ', n2, (n1-n2)/(n1+0.0)
				endif

				call drop_pixels2(pos_list0, rho_list0, drho_list0, markdrop, pos_list, rho_list, drho_list)
				call get_mu_from_gradient_list(pos_list, drho_list, drho_mu_data)
				chisqlist(idrop) = chisq_of_mu_data2(drho_mu_data)

				deallocate(pos_list,rho_list,drho_list,drho_mu_data)
!				print *, 'gd_mldprho_chi2s J'
			endif
		enddo
				
		call cpu_time(time2)		

		deallocate(absdrholist,markdrop,reflist)

		if(print_time .or. gradfieldprintinfo) then
			print *, '  Time used in grid_rho_drho_list: ', real(time1-time0)
			print *, '  Time used in drop: ', real(time2-time1)
		endif
		gradfieldprintinfo = .false.
!		print *, '(gd_mldprho_chi2s) Test E...'		
	end subroutine gd_mldprho_chi2s
	

  !------------------------------------------
  ! Mark the pixels that shall be dropped
  !------------------------------------------
	subroutine mark_drop_pixels(reflist, markdrop, n, i1, i2)
		!dummy
		real(dl) :: reflist(n)
		integer :: markdrop(n), n, i1, i2
		!local
		integer :: i, indexarray(n)
		if(i1.le.0 .and. i2.ge.n) then
			print *, 'Warning (mark_drop_pixels)! No need to drop: ', i1, i2
			return
		endif
		do i = 1, n
			indexarray(i) = i
		enddo
		call Qsort2(reflist, indexarray, n)
		if(i1 .ge. 1) then
			do i = 1, i1
				markdrop(indexarray(i)) = 1
			enddo
		endif
		if(i2 .le. n) then
			do i = i2, n
				markdrop(indexarray(i)) = 1
			enddo
		endif
	end subroutine mark_drop_pixels
	

  !------------------------------------------
  ! dropping pixels in lists
  !------------------------------------------
  	subroutine drop_pixels(val_list, pos_list, dval_list, markdrop)
  		real(dl), allocatable :: val_list(:), distance_list(:), pos_list(:,:), dval_list(:,:), reflist(:), tmp(:,:)
  		integer :: markdrop(:)
  		integer :: i, j,  n, m
		n = size(val_list)
		if(size(pos_list,2).ne.n .or. size(dval_list,2).ne. n.or.size(markdrop).ne.n) then
			print *, 'ERROR (drop_pixels)! Check length of arrays: ', n,size(pos_list,2),size(dval_list,2),size(markdrop)
			stop
		endif
		m = 0
		do i = 1, n
			if(markdrop(i) .eq. 0) m = m+1
		enddo
		if(m.eq.n) return
		allocate(tmp(7,n))
		do i = 1, n
			tmp(1,i)=val_list(i); 
			tmp(2:4,i)=pos_list(1:3,i); tmp(5:7,i)=dval_list(1:3,i)	
		enddo
		deallocate(val_list,dval_list,pos_list)
		allocate(val_list(m),dval_list(3,m),pos_list(3,m))
		j = 0
		do i = 1,n
			if(markdrop(i) .eq. 0) then
				j = j+1
				val_list(j) = tmp(1,i)
				pos_list(1:3,j) = tmp(2:4,i)
				dval_list(1:3,j) = tmp(5:7,i)
			endif
		enddo
		deallocate(tmp)
	end subroutine drop_pixels

  !------------------------------------------
  ! dropping pixels in lists
  !------------------------------------------
  	subroutine drop_pixels2(pos_list, val_list, dval_list, markdrop, pos_list2, val_list2, dval_list2)
		! Dummy
  		real(dl), intent(in) :: pos_list(:,:), val_list(:), dval_list(:,:)
  		integer, intent(in) :: markdrop(:)
  		real(dl), intent(out), allocatable :: pos_list2(:,:), val_list2(:), dval_list2(:,:)
		! local
  		integer :: i, j,  n, m
		n = size(val_list)
		if( size(pos_list,2).ne.n .or. size(dval_list,2).ne. n.or.size(markdrop).ne.n) then
			print *, 'ERROR (drop_pixels2)! Check length of arrays: ', n,size(pos_list,2),size(dval_list,2),size(markdrop)
			stop
		endif
		m = 0
		do i = 1, n
			if(markdrop(i) .eq. 0) m = m+1
		enddo
		allocate(pos_list2(3,m),val_list2(m),dval_list2(3,m))
		j = 0
		do i = 1,n
			if(markdrop(i) .eq. 0) then
				j = j+1
				pos_list2(1:3,j) = pos_list(1:3,i)	
				val_list2(j) = val_list(i)
				dval_list2(1:3,j) = dval_list(1:3,i)	
			endif
		enddo
	end subroutine drop_pixels2


  !------------------------------------------
  ! dropping pixels in lists
  !------------------------------------------
  	subroutine sigma_drop_pixels(val_list, pos_list, absdrholist, markdrop, n, dropval, dropvalsigma, dropdval, dropdvalsigma )
		!DUMMY
  		real(dl), intent(in) :: val_list(n), pos_list(3,n),  absdrholist(n)
  		integer, intent(in) :: n
  		integer, intent(inout) :: markdrop(n)
  		real(dl) :: dropvalsigma(2), dropdvalsigma(2)
  		logical :: dropval, dropdval
  		! LOCAL
  		real(dl) :: x, valmean, valvar, dvalmean, dvalvar, var1, var2, dvar1, dvar2
  		integer :: i


		if(.not. dropval .and. .not. dropdval) return

		if(dropval) then
			call get_mean_var(val_list, valmean, valvar)
			var1 = valmean - sqrt(valvar)*dropvalsigma(1)
			var2 = valmean + sqrt(valvar)*dropvalsigma(2)
			do i = 1, n
				x = val_list(i)
				if(x<var1 .or. x>var2) markdrop(i) = 1
			enddo
!			if(gbtp) &
!				print*, '   sigma_drop_pixels: Dropping var according to (mean, var, var1, var2):', &
!					real(valmean), real(valvar), real(var1), real(var2), '...'
		endif
		
		if(dropdval) then
			call get_mean_var(absdrholist, dvalmean, dvalvar)
			dvar1 = dvalmean - sqrt(dvalvar)*dropdvalsigma(1)
			dvar2 = dvalmean + sqrt(dvalvar)*dropdvalsigma(2)
			do i = 1, n
				x = absdrholist(i)
				if(x<dvar1 .or. x>dvar2) markdrop(i) = 1
			enddo
!			if(gbtp) &
!				print*, '   sigma_drop_pixels: Dropping dvar according to (mean, var, var1, var2):', &
!					real(dvalmean), real(dvalvar),  real(dvar1), real(dvar2), '...'			
		endif

	end subroutine sigma_drop_pixels


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
	



