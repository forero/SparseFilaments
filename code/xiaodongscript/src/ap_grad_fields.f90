
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

!!! IMPORTANT SETTINGS	
! using fixed radius smoothing kernel
	logical, public :: gb_use_fixmd = .true.
	real(dl), public :: gb_fixmd = 40.0_dl

!!! IMPORTANT SETTINGS
! using cuts in smooth kernel (ignore nearby halos)
	logical, public :: gb_do_seg_cut = .false.
	real(dl), public :: gb_seg_cut_dist = 5.0_dl

!! IMPORTANT SETTINGS	
	integer, parameter :: nbinchisq = 4
	integer, parameter :: nbinchisqlist(nbinchisq) = (/2,3,4,5/)
	
	logical, public :: gb_dropmethod(2) = .false. ! 1 means ratio; 2 means variance
	logical, public :: gb_dropval(2), gb_dropdval(2)
	real(dl), public :: gb_dropvalratio(2,2), gb_dropdvalratio(2,2)
	logical, public :: gradfieldprintinfo = .true.
	
	
!############################## Commenting testing codes
	!BEGIN TESTING
	character(len=char_len) :: tsbestr
	!END TESTING
	
	
contains	

  !------------------------------------------
  ! estimating val/dval by using given list
  !------------------------------------------
	subroutine get_val_dval_list(smnum, print_info, pos_list, val_list, dval_list, n, save_in_cellmat)
		! Dummy
		integer, intent(in) :: smnum
		logical  :: print_info
		integer, intent(in) :: n
		real(dl), intent(in) :: pos_list(3,n)
		real(dl), intent(out) :: val_list(n), dval_list(3,n)
		logical, optional, intent(in) :: save_in_cellmat
		! LOCAL
		real(dl) :: x,y,z,r,val,dvalx,dvaly,dvalz,max_dist
		integer :: i, j, ix, iy, iz, num
		logical :: touchbdflag,vcorect_er_flag
		
		if(print_info) then
			write(*,'(A,i9)'), '   (get_val_dval_list begins) Calculating val/dval based on pos_list. Total #: ', n
			if(present(save_in_cellmat)) then
				write(*,'(26x,A)') 'Results will be saved into gb_cell_mat.'
			else
				write(*,'(26x,A)') 'Results will NOT be saved into gb_cell_mat.'
			endif
			write(*,'(26x,A,L2,L2)') 'use_fixmd / do_seg_cut = ', gb_use_fixmd, gb_do_seg_cut
			if(gb_use_fixmd) then
				write(*,'(26x,A,f10.3)') 'Fixed smoothing length gb_fixmd = ', real(gb_fixmd)
				if(gb_keepzerorho) then
					write(*,'(26x,A,i3,A)') 'Minimal in-sphere #:', &
						gb_min_smnum, '. Small-# pixels will be saved with rho=0.'
				else
					write(*,'(26x,A,i3,A)') 'Minimal in-sphere #:', &
						gb_min_smnum, '. Small-# pixels will be skipped with rhodrho_has_er=T'
				endif
			else
				write(*,'(26x,A,i4)') 'Use # of nearby halos = ', smnum
			endif
			if(gb_do_seg_cut)  then
				write(*,'(26x,A,f10.3)') 'Applying minimal-r cut: ', real(gb_seg_cut_dist)
			endif
		endif
		
		do i = 1, n
			x=pos_list(1,i); y=pos_list(2,i); z=pos_list(3,i)
			if(gb_use_fixmd) then
				call nb_fixmd_list(x,y,z,num,val,dvalx,dvaly,dvalz,&
					gb_fixmd,gb_do_seg_cut,gb_seg_cut_dist,touchbdflag,vcorect_er_flag)!,.true.) !Testing
			else
				call nb_list0(x,y,z,smnum,val,dvalx,dvaly,dvalz,max_dist)
			endif
			if(touchbdflag.or.vcorect_er_flag) then
				print *, 'ERROR (get_val_dval_list)!!! There shall be no error!: touchbdflag, vcorect_er_flag = ', &
					touchbdflag, vcorect_er_flag
			endif
			val_list(i) = val
			dval_list(1:3,i) = (/dvalx,dvaly,dvalz/)
			if(present(save_in_cellmat)) then
				if(save_in_cellmat) then
					call cell_index(ix,iy,iz,x,y,z)
					gb_cell_mat(ix,iy,iz)%rhodrhos(1) = val
					gb_cell_mat(ix,iy,iz)%rhodrhos(2) = dvalx
					gb_cell_mat(ix,iy,iz)%rhodrhos(3) = dvaly
					gb_cell_mat(ix,iy,iz)%rhodrhos(4) = dvalz
				endif
			endif		
		enddo	
		if(print_info) print *, '  (get_val_dval_list done)'
	end subroutine get_val_dval_list


  !------------------------------------------
  ! estimating rho/drho at grid points
  !------------------------------------------
	subroutine grid_rho_drho_list(smnum, printinfo, pos_list, rho_list, drho_list, ixiyizlist)
		!Dummy
		integer, intent(in) :: smnum 
		logical, intent(in) :: printinfo
		real(dl), allocatable, intent(out) :: pos_list(:,:),  rho_list(:),  drho_list(:,:)
		integer, allocatable, optional, intent(out) :: ixiyizlist(:,:)
		! Local
		real(dl) :: x,y,z,r,rho,drhox,drhoy,drhoz,max_dist
		integer :: i,j, ix, iy, iz, n, n1, n2, fixmdsmnum
		real(dl) :: numhasbd,numhasver,numsmallsmnum, numavg, numsqavg, countnum
		logical :: touchbdflag, vcorect_er_flag, checkbd = .true.
		real(dl) :: bdratio = 1.0_dl
		
		if(printinfo) then
			print *, '  (grid_rho_drho_list begins) Estimating rho_list and its gradient...'
			if(present(ixiyizlist)) then
				write(*,'(26x,A)') 'Output ix/iy/iz for each pixel...'
			else
				write(*,'(26x,A)') 'Do not output ix/iy/iz...'
			endif
			write(*,'(26x,A,L2,L2)') 'use_fixmd / gb_do_seg_cut = ', gb_use_fixmd, gb_do_seg_cut
			if(gb_use_fixmd) then
				write(*,'(26x,A,f10.3)') 'Fixed smoothing length gb_fixmd = ', real(gb_fixmd)
				if(gb_keepzerorho) then
					write(*,'(26x,A,i3,A)') 'Minimal in-sphere #:', &
						gb_min_smnum, '. Small-# pixels will be saved with rho=0.'
				else
					write(*,'(26x,A,i3,A)') 'Minimal in-sphere #:', &
						gb_min_smnum, '. Small-# pixels will be skipped with has_bd_effect=T'
				endif
			else
				write(*,'(26x,A,i4)') 'Use # of nearby halos = ', smnum
			endif
			if(gb_do_seg_cut)  then
				write(*,'(26x,A,f10.3)') 'Applying minimal-r cut: ', real(gb_seg_cut_dist)
			endif
		endif
		
		if(printinfo.and.checkbd) then
			write(*,'(26x,A,f10.5,2x,f10.5,A)') 'Check boundary using rmin, rmax = ', gbrmin,gbrmax, '.'
		endif

		n = gb_n_cellx*gb_n_celly*gb_n_cellz
		
		! n1 counts how many pixels in shell; n2 counts how many pixels in shell and has no boundary effect
		i=0; n1=0; n2=0;
		numavg = 0.0;
		numsqavg = 0.0;
		numsmallsmnum = 0.0;
		numhasbd = 0.0;
		numhasver = 0.0;
		do ix = 1, gb_n_cellx
		do iy = 1, gb_n_celly
		do iz = 1, gb_n_cellz
			i=i+1;
			call cell_pos(ix,iy,iz,x,y,z)
			r = sqrt(x*x+y*y+z*z)
			if(checkbd) then
				if(r<gbrmin.or.r>gbrmax) then
					gb_cell_mat(ix,iy,iz)%rhodrho_has_er = .true.
					cycle
				endif
			endif
			n1=n1+1
			! Estimating rho / drhos ...
			if(gb_use_fixmd) then
				if(checkbd) then
					if(has_boundary_effect(x,y,z,gb_fixmd,bdratio)) then
						gb_cell_mat(ix,iy,iz)%rhodrho_has_er = .true.
						numhasbd = numhasbd + 1.0
						cycle
					endif
				endif	
				call nb_fixmd_list(x,y,z,fixmdsmnum,rho,drhox,drhoy,drhoz,&
					gb_fixmd,gb_do_seg_cut,gb_seg_cut_dist,touchbdflag,vcorect_er_flag)
				! skip points which has error in v correction
				if(vcorect_er_flag) then
					gb_cell_mat(ix,iy,iz)%rhodrho_has_er = .true.
					numhasver = numhasver + 1.0
					cycle
				endif
				if(fixmdsmnum .le. gb_min_smnum) then
					numsmallsmnum = numsmallsmnum+1.0_dl
					if(.not.gb_keepzerorho) then
						gb_cell_mat(ix,iy,iz)%rhodrho_has_er = .true.
						cycle
					endif
				endif
				if(touchbdflag .eq. .true.) then
					print *, 'ERROR (grid_rho_drho_list)!!! We expect no touchbd affair:', x,y,z
					gb_cell_mat(ix,iy,iz)%rhodrho_has_er = .true.
					numhasbd = numhasbd + 1.0
					cycle
!					stop
				endif
				numavg = numavg + fixmdsmnum
				numsqavg = numsqavg + fixmdsmnum**2.0
			else
				call nb_list0(x,y,z,smnum,rho,drhox,drhoy,drhoz,max_dist)
				if(checkbd) then
					if(has_boundary_effect(x,y,z,max_dist,bdratio)) then
						gb_cell_mat(ix,iy,iz)%rhodrho_has_er = .true.
						cycle
					endif
				endif	
			endif
			
			gb_cell_mat(ix,iy,iz)%rhodrho_has_er = .false.
			n2=n2+1	
			gb_cell_mat(ix,iy,iz)%rhodrhos(1) = rho
			gb_cell_mat(ix,iy,iz)%rhodrhos(2) = drhox
			gb_cell_mat(ix,iy,iz)%rhodrhos(3) = drhoy
			gb_cell_mat(ix,iy,iz)%rhodrhos(4) = drhoz
		enddo
		enddo
		enddo

		if(printinfo) then
			if(gb_use_fixmd) then
				numavg = dble(numavg) / dble(n2)
				write(*,'(26x,i7,A,f5.2,A,i3,A)') int(numsmallsmnum+0.1), ' (',numsmallsmnum/dble(n1)*100.0,&
					'%) has in-sphere-# <= ',gb_min_smnum,'.'
				write(*,'(26x,i7,A,f5.2,A)') int(numhasbd+0.1), ' (',numhasbd/dble(n1)*100.0,&
					'%) skipped due to bd effect.'
				write(*,'(26x,i7,A,f5.2,A)') int(numhasver+0.1), ' (',numhasver/dble(n1)*100.0,&
					'%) skipped due to error in vel-correction.'
				write(*,'(26x,A,f11.4,A,f11.4)') 'Average # of in-sphere halos = ', &
					 real(numavg), '  +/-', real(sqrt(numsqavg/dble(n2)-numavg**2.0))
			endif
		endif

		allocate(pos_list(3,n2),rho_list(n2),drho_list(3,n2))
		i = 0
		do ix = 1, gb_n_cellx
		do iy = 1, gb_n_celly
		do iz = 1, gb_n_cellz
			if(.not.gb_cell_mat(ix,iy,iz)%rhodrho_has_er) then
				i =i+1
				call cell_pos(ix,iy,iz,x,y,z)
				pos_list(1,i)=x; pos_list(2,i)=y; pos_list(3,i)=z; 
				rho_list(i) = gb_cell_mat(ix,iy,iz)%rhodrhos(1)
				drho_list(1:3,i) = gb_cell_mat(ix,iy,iz)%rhodrhos(2:4)
			endif
		enddo
		enddo
		enddo

		if(present(ixiyizlist)) then 
			allocate(ixiyizlist(3,n2))
			do i = 1, n2
				call cell_index(ixiyizlist(1,i),ixiyizlist(2,i),ixiyizlist(3,i),&
					pos_list(1,i),pos_list(2,i),pos_list(3,i))
			enddo
		endif

		if(printinfo) then
			write(*,'(24x,A,1x,i9)')    'Total # of pixels: ', n2
			write(*,'(26x,A,f8.4)')     'ratio of boundary skipped = ', dble(n1-n2)/(n1+0.0)
			write(*,'(26x,A,f8.4)')     'ratio-of-used in cell     = ', dble(n2)/dble(gb_n_cellx*gb_n_celly*gb_n_cellz)
			print *, '  (grid_rho_drho_list done)'
		endif
	end subroutine grid_rho_drho_list

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
!			print *, 'Warning (mark_drop_pixels)! No need to drop: ', i1, i2
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
