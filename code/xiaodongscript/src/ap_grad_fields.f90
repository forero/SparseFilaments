
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
	
	logical, public :: gb_dropmethod(2) = .false. ! 1 means ratio; 2 means variance
	logical, public :: gb_dropval(2), gb_dropdval(2)
	real(dl), public :: gb_dropvalratio(2,2), gb_dropdvalratio(2,2)
	logical, public :: gradfieldprintinfo = .true.
	integer, parameter :: gb_numwei = 5
	
!##### information of a pixel
	type :: pixelinfo
		real :: rho, drhox, drhoy, drhoz
		logical :: has_bd_effect, has_drho
	end type
	
!############################## Commenting testing codes
	character(len=char_len) :: tsbestr, tsbefilestr
	!BEGIN TESTING
	integer :: testi, nbfileunit = 666
	real(dl) :: nowmu, muleft = 0.103_dl, muright = 0.1031_dl, testxyzarray(10000)
	!END TESTING
	
contains	

  !------------------------------------------
  ! estimating val/dval by using given list
  !------------------------------------------
	subroutine get_val_dval_list(smnum, print_info, pos_list, val_list, dval_list, save_in_cellmat)
		! Dummy
		integer :: smnum
		logical  :: print_info
		real(dl), allocatable, intent(in) :: pos_list(:,:)
		real(dl), allocatable, intent(out) :: val_list(:), dval_list(:,:)
		logical, optional, intent(in) :: save_in_cellmat
		! LOCAL
		real(dl) :: x,y,z,r,val,dvalx,dvaly,dvalz,max_dist
		integer :: i, j, ix, iy, iz, n

		if(size(pos_list,1).ne.3) then
			print *, '  ERROR (get_val_dval_list)! dim not 3!'
			stop
		endif
		
		n=size(pos_list,2)
		
		if(print_info) then
			print *, '  (get_val_dval_list) Calculating val/dval list ...'
		endif

		allocate(val_list(n),dval_list(3,n))

		do i = 1, n
			x=pos_list(1,i); y=pos_list(2,i); z=pos_list(3,i)
!			r = sqrt(x*x+y*y+z*z)
			
			if(use_fixmd) then
				call nb_fixmd_list(x,y,z,smnum,val,dvalx,dvaly,dvalz)
			elseif(use_seg) then
				call nb_seg_list(x,y,z,smnum,val,dvalx,dvaly,dvalz,max_dist)
			else
				call nb_list0(x,y,z,smnum,val,dvalx,dvaly,dvalz,max_dist)
			endif
			val_list(i) = val
			dval_list(1:3,i) = (/dvalx,dvaly,dvalz/)
			if(present(save_in_cellmat)) then
				if(save_in_cellmat) then
					call cell_index(ix,iy,iz,x,y,z)
					gb_cell_mat(ix,iy,iz)%rho = val
				endif
			endif
		enddo			
	end subroutine get_val_dval_list


  !------------------------------------------
  ! estimating list of rho/drho by using
  !  finite difference
  ! To be tested...
  !------------------------------------------
	subroutine fd_rhodrho_list(RSD, AP, cs, changenuminx, pos_list, rho_list, drho_list)
		!Dummy
		type(chisq_settings) :: cs
		integer, intent(in) :: RSD, AP
		real(dl), intent(in) :: changenuminx
		real(dl), allocatable, intent(out) :: pos_list(:,:),  rho_list(:),  drho_list(:,:)
		!Local
		type(pixelinfo), allocatable :: pixelgrid(:,:,:)
		integer :: i, ix,iy,iz, n1,n2, smnum, di
		real(dl) :: x,y,z, rho,drhox,drhoy,drhoz, r,max_dist, numavg,numsqavg
		
		call do_cell_init(RSD, AP, real(cs%num_in_x)*(1.0_dl+changenuminx), cs%print_info)
		
		if(cs%use_num_density) then
			gb_mass_list = 1.0_dl
		else
			gb_mass_list = gb_bf_mass_list
		endif
		
		if(cs%check_boundary .and. cs%print_info) &
				write(*,'(2x,A,f10.5,2x,f10.5,A)') &
					' Chenck boundary using rmin, rmax = ', gbrmin, gbrmax, ' ...'

		if(cs%print_info) &
			print *, '  Estimating rho_list and its gradient...'

		allocate(pixelgrid(gb_n_cellx,gb_n_celly,gb_n_cellz))

		i=0; n1=0; n2=0;		
		do ix = 1, gb_n_cellx
		do iy = 1, gb_n_celly
		do iz = 1, gb_n_cellz
			i=i+1;
			call cell_pos(ix,iy,iz,x,y,z)
			r = dsqrt(x*x+y*y+z*z)
			if(r<gbrmin.or.r>gbrmax) then
				pixelgrid(ix,iy,iz)%has_bd_effect = .true.
				cycle
			endif
			n1=n1+1
			! Estimate rho. Divided into two cases (whether use fixed radius of smoothing kernel)
			if(use_fixmd) then
				if(has_boundary_effect(x,y,z,gb_fixmd,cs%cb_adjust_ratio)) then
					pixelgrid(ix,iy,iz)%has_bd_effect = .true.
					cycle
				endif
				call nb_fixmd_rho(x,y,z,smnum,rho)
				if(smnum .eq. 0) then
					pixelgrid(ix,iy,iz)%has_bd_effect = .true. ! in fact not boundary effect but no in-cell halo..
					cycle
				endif
				numavg = numavg + smnum
				numsqavg = numsqavg + smnum**2.0
				max_dist = gb_fixmd
			else
				if(use_seg) then
					call nb_seg_list(x,y,z,cs%smnum,rho,drhox,drhoy,drhoz,max_dist)
				else
					call nb_list0(x,y,z,cs%smnum,rho,drhox,drhoy,drhoz,max_dist)
				endif
				if(has_boundary_effect(x,y,z,max_dist,cs%cb_adjust_ratio)) then
					pixelgrid(ix,iy,iz)%has_bd_effect = .true.
					cycle
				endif
			endif
			pixelgrid(ix,iy,iz)%has_bd_effect = .false.
			pixelgrid(ix,iy,iz)%rho = rho
			n2=n2+1	
		enddo
		enddo
		enddo
		
		if(use_fixmd) then
			numavg = numavg / dble(n2)
			print *, '# of nearby halos = ', real(numavg), '+/-', real(dsqrt(numsqavg/dble(n2)-numavg**2.0))
		endif
		
!	This define which neareast neighbor pixel to estimate derivative		
		di = 1
		do ix = 1+di, gb_n_cellx-di
		do iy = 1+di, gb_n_celly-di
		do iz = 1+di, gb_n_cellz-di
			if(pixelgrid(ix-di,iy,iz)%has_bd_effect .or. &
 			   pixelgrid(ix,iy-di,iz)%has_bd_effect .or. &
			   pixelgrid(ix,iy,iz-di)%has_bd_effect .or. &
			   pixelgrid(ix+di,iy,iz)%has_bd_effect .or. &
			   pixelgrid(ix,iy+di,iz)%has_bd_effect .or. &
			   pixelgrid(ix,iy,iz+di)%has_bd_effect) then
			   	pixelgrid(ix,iy,iz)%has_drho = .false.
			   	cycle
			endif
			pixelgrid(ix,iy,iz)%has_drho = .true.
			pixelgrid(ix,iy,iz)%drhox = pixelgrid(ix+di,iy,iz)%rho - pixelgrid(ix-di,iy,iz)%rho
			pixelgrid(ix,iy,iz)%drhoy = pixelgrid(ix,iy+di,iz)%rho - pixelgrid(ix,iy-di,iz)%rho
			pixelgrid(ix,iy,iz)%drhoz = pixelgrid(ix,iy,iz+di)%rho - pixelgrid(ix,iy,iz-di)%rho
		enddo
		enddo
		enddo
		
		
	end subroutine fd_rhodrho_list


!----------------------------------------------------
!----------------------------------------------------
!  Further tries (multiple calculation of chisq)
!...

  !------------------------------------------
  ! estimating rho/drho at grid points
  !------------------------------------------
	subroutine grid_rho_drho_list(RSD, AP, cs, changenuminx, pos_list, rho_list, drho_list, boundary_rmin, boundary_rmax, ixiyizlist)
		!Dummy
		type(chisq_settings) :: cs
		integer, intent(in) :: RSD, AP 
		real(dl), intent(in) :: changenuminx
		real(dl), allocatable, intent(out) :: pos_list(:,:),  rho_list(:),  drho_list(:,:)
		real(dl), intent(out) :: boundary_rmin, boundary_rmax
		integer, allocatable, optional, intent(out) :: ixiyizlist(:,:)
		! Local
		real(dl) :: x,y,z,r,rho,drhox,drhoy,drhoz,max_dist
		real(dl), allocatable :: tmp(:,:), tmpxyzrlist(:,:,:), tmpxyzrseg(:,:)
		integer :: i,j, ix, iy, iz, n, n1, n2, smnum
		real(dl) :: numavg, numsqavg, countnum
		
		call do_cell_init(RSD, AP, real(cs%num_in_x)*(1.0_dl+changenuminx), cs%print_info)
		
		if(cs%use_num_density) then
			gb_mass_list = 1.0_dl
		else
			gb_mass_list = gb_bf_mass_list
		endif
		
		if(cs%print_info) print *, '  (grid_rho_drho_list) Estimating rho_list and its gradient...'

		if(cs%check_boundary) then
			boundary_rmin = gbrmin 
			boundary_rmax = gbrmax
			if(cs%print_info) then
				write(*,'(2x,A,f10.5,2x,f10.5,A)') &
					' (grid_rho_drho_list) Chenck boundary using rmin, rmax = ', boundary_rmin, boundary_rmax, ' ...'
			endif
		endif

		n = gb_n_cellx*gb_n_celly*gb_n_cellz
		allocate(tmp(7,n))
		
		! n1 counts how many pixels in shell; n2 counts how many pixels in shell and has no boundary effect
		i=0; n1=0; n2=0;
		numavg = 0.0;
		numsqavg = 0.0;
		do ix = 1, gb_n_cellx
		do iy = 1, gb_n_celly
		do iz = 1, gb_n_cellz
			i=i+1;
			call cell_pos(ix,iy,iz,x,y,z)
			r = dsqrt(x*x+y*y+z*z)
			if(cs%check_boundary) then
				if(r<boundary_rmin.or.r>boundary_rmax) then
					gb_cell_mat(ix,iy,iz)%has_bd_effect = .true.
					cycle
				endif
			endif
			n1=n1+1
			if(use_fixmd) then
				if(cs%check_boundary) then
					if(has_boundary_effect(x,y,z,gb_fixmd,cs%cb_adjust_ratio)) then
						gb_cell_mat(ix,iy,iz)%has_bd_effect = .true.
						cycle
					endif
				endif	
				call nb_fixmd_list(x,y,z,smnum,rho,drhox,drhoy,drhoz)
				if(smnum .le. 5) then
!					print *, 'Skip pixel with smnum = ',smnum,': x,y,z = ', x,y,z
					gb_cell_mat(ix,iy,iz)%has_bd_effect = .true.
					cycle
				endif
				gb_cell_mat(ix,iy,iz)%has_bd_effect = .false.
				numavg = numavg + smnum
				numsqavg = numsqavg + smnum**2.0
				max_dist = gb_fixmd
			else
				if(use_seg) then
					call nb_seg_list(x,y,z,cs%smnum,rho,drhox,drhoy,drhoz,max_dist)
				else
					call nb_list0(x,y,z,cs%smnum,rho,drhox,drhoy,drhoz,max_dist)
				endif
				if(cs%check_boundary) then
					if(has_boundary_effect(x,y,z,max_dist,cs%cb_adjust_ratio)) then
						gb_cell_mat(ix,iy,iz)%has_bd_effect = .true.
						cycle
					endif
				endif	
			endif
			
			n2=n2+1	
			tmp(1:7,n2) = (/x,y,z,rho,drhox,drhoy,drhoz/)
			gb_cell_mat(ix,iy,iz)%rho = rho
			
		enddo
		enddo
		enddo
		
		numavg = numavg / dble(n2)
		print *, '# of nearby halos = ', real(numavg), '+/-', real(dsqrt(numsqavg/dble(n2)-numavg**2.0))
		
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
		
		if(present(ixiyizlist)) then 
			allocate(ixiyizlist(3,n2))
			do i = 1, n2
				call cell_index(ixiyizlist(1,i),ixiyizlist(2,i),ixiyizlist(3,i),&
					pos_list(1,i),pos_list(2,i),pos_list(3,i))
			enddo
		endif
		
		if(cs%print_info) then
			print *, '  tot within boundary / tot without boundary effect = ', n1, n2
			print *, '  ratio of skipped = ', (n1-n2)/(n1+0.0)
		endif
	end subroutine grid_rho_drho_list
	
  !------------------------------------------
  ! estimating gradient of fields of 
  !  rho/delta/normed delta at grid points
  !------------------------------------------
	subroutine gd_mldprho_chi2s(RSD, AP, cs, changenuminx, chisqlist, dfchisqlist, numchisq)
		! Dummy
		integer, intent(in) :: RSD, AP, numchisq
		type(chisq_settings) :: cs
		real(dl), intent(in) :: changenuminx
		real(dl), intent(out) :: chisqlist(numchisq), dfchisqlist(numchisq)
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
		!variables used to estimate mu at different r bins
		real(dl), allocatable :: r_list(:), muavlist(:), muerlist(:), ravlist(:), muvarlist(:)
		real(dl) :: rmin, rmax
		integer :: num_bins
		!tmp variables used for testing
		real(dl) :: xyz1,xyz2,r1,r2,tsbedropratio,nowx,nowy,nowz,nowrsq,rmaxsq,invrsq

		if(numchisq.ne.cs%numdrop) then
			print *, 'ERROR (gd_mldprho_chi2s)! size of chisqlist mismathces than ', cs%numdrop
			stop
		endif

		call cpu_time(time0)
		
		if(cs%print_info) print *, 'Estimating rho, delta, normed_delta: changenuminx = ', real(changenuminx), '...'

!		print *, '(gd_mldprho_chi2s) Test A...'
		call grid_rho_drho_list(RSD, AP, cs, changenuminx, pos_list0, rho_list0, drho_list0, &
				boundary_rmin, boundary_rmax)
!		print *, '(gd_mldprho_chi2s) Test B...'

		call cpu_time(time1)		
		
		n1 = size(rho_list0)
!		print *, '(gd_mldprho_chi2s) Test C...'
		allocate(absdrholist(n1),markdrop(n1),reflist(n1))
		markdrop = 0
		do i = 1, n1
			absdrholist(i)=sqrt(drho_list0(1,i)**2.0+drho_list0(2,i)**2.0+drho_list0(3,i)**2.0)
		enddo

		num_bins = 2
		allocate(muavlist(num_bins), muerlist(num_bins), ravlist(num_bins), muvarlist(num_bins))
		do idrop = 1, cs%numdrop
			! cycle if no drop
			if(.not.cs%dropval(idrop) .and. .not.cs%dropdval(idrop)) then
				call list_get_mu(pos_list0, drho_list0, drho_mu_data)
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
				call list_get_mu(pos_list, drho_list, drho_mu_data)
				chisqlist(idrop) = chisq_of_mu_data2(drho_mu_data)
				
				! chisq from differantial of mu
				allocate(r_list(size(rho_list)))
				do i = 1, size(r_list)
					r_list(i) = sqrt(pos_list(1,i)**2.0+pos_list(2,i)**2.0+pos_list(3,i)**2.0)
				enddo
				call find_min_max(r_list, size(r_list), rmin, rmax)
								
				do i = 1, size(drho_mu_data)
					drho_mu_data(i) = abs(drho_mu_data(i))
				enddo
			  	call binned_quan(drho_mu_data, r_list, rmin, rmax, size(drho_mu_data), num_bins, &	
  					muavlist, muerlist, ravlist, muvarlist)
				dfchisqlist(idrop) = abs(muavlist(2) - muavlist(1))**2.0 / ((muerlist(1)**2.0 + (muerlist(2)**2.0)))
				if(gradfieldprintinfo .and. .false.) then
	  				print *
					print *, 'num_bins: ', num_bins
					print *, 'rmin/rmax: ', rmin, rmax
					print *, 'muavlist: ', muavlist
					print *, 'muerlist: ', muerlist
					print *, 'ravlist: ', ravlist
					print *
					print *, 'chisq of diff mu: ', dfchisqlist(idrop)
					print *
				endif
				deallocate(r_list,pos_list,rho_list,drho_list,drho_mu_data)
!				print *, 'gd_mldprho_chi2s J'
!				print *, chisqlist
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
  ! estimating gradient of fields of 
  !  rho/delta/normed delta at grid points
  !------------------------------------------
	subroutine gd_mldprho_mlchi2s(RSD, AP, cs, changenuminx, rhochisqlist, rhodfchisqlist, &
		dltchisqlist, dltdfchisqlist, ndltchisqlist, ndltdfchisqlist, numchisq)
		! Dummy
		integer, intent(in) :: RSD, AP, numchisq
		type(chisq_settings) :: cs
		real(dl), intent(in) :: changenuminx
		real(dl), intent(out) :: rhochisqlist(numchisq), rhodfchisqlist(numchisq), dltchisqlist(numchisq), &
			dltdfchisqlist(numchisq), ndltchisqlist(numchisq), ndltdfchisqlist(numchisq)
		! Local
		real(dl) :: boundary_rmin, boundary_rmax
		integer :: i, j, k, idrop
		integer :: n1, n2
		!variables used to save the mu data
		real(dl), allocatable :: drho_mu_data(:), ddlt_mu_data(:), dndlt_mu_data(:)
		!vairables used to save info before drop
		real(dl), allocatable :: poslist0(:,:), rholist0(:), drholist0(:,:), &
			dltlist0(:), ddltlist0(:,:), ndltlist0(:), dndltlist0(:,:)
		!variables used to save info after drop
		real(dl), allocatable :: poslist(:,:), rholist(:), drholist(:,:), &
			dltlist(:), ddltlist(:,:), ndltlist(:), dndltlist(:,:)
		!variables used to drop pixels
		real(dl), allocatable :: refvallist(:), refdvallist(:), refvallist0(:), refdvallist0(:)
		integer, allocatable :: markdrop(:)
		!variables used to re-scale mass 
		real(dl), allocatable :: distancelist(:)
		integer :: nsm = 1 !do smoothing for only 1 time
		!variables used to get time consumed
		logical :: print_time = .false.
		real(dl) :: tmp, time0, time1, time2, time3, time4!, time5, time6, time7
		!variables used to estimate mu at different r bins
		real(dl), allocatable :: rlist(:), muavlist(:), muerlist(:), ravlist(:), muvarlist(:)
		real(dl) :: rmin, rmax
		integer :: num_bins
		!tmp variables used for testing
		real(dl) :: xyz1,xyz2,r1,r2,tsbedropratio,nowx,nowy,nowz,nowrsq,rmaxsq,invrsq

		if(numchisq.ne.cs%numdrop) then
			print *, 'ERROR (gd_mldprho_chi2s)! size of chisqlist mismathces than ', cs%numdrop
			stop
		endif
		
		! Get rho/drho
		call cpu_time(time0)
		if(cs%print_info) print *, 'Estimating rho, delta, normed_delta: changenuminx = ', real(changenuminx), '...'
		call grid_rho_drho_list(RSD, AP, cs, changenuminx, poslist0, rholist0, drholist0, &
			boundary_rmin, boundary_rmax)

		! Rescale mass list: normalized according to <rho>(r)
		n1 = size(rholist0)
		allocate(distancelist(n1))
		do i = 1, n1
			distancelist(i) = sqrt(poslist0(1,i)**2.0+poslist0(2,i)**2.0+poslist0(3,i)**2.0)
		enddo
		call avg_mlist(distancelist, rholist0, n1, cs, 1)

		! Get dlt/ddlt
		call cpu_time(time1)		
		call get_val_dval_list(cs%smnum, cs%print_info, poslist0, dltlist0, ddltlist0)

		! Rescal mass list: normalized according to <dlt dlt>(r)
		call sqvar_mlist(distancelist, dltlist0, n1, cs, 1)
	
		! Get ndlt/dndlt
		call cpu_time(time2)		
		call get_val_dval_list(cs%smnum, cs%print_info, poslist0, ndltlist0, dndltlist0)

		! Preparing for dropping
		call cpu_time(time3)
		allocate(refvallist(n1),refdvallist(n1),refvallist0(n1),refdvallist0(n1),markdrop(n1))	

		markdrop = 0
!		refvallist0 = rholist0
		refvallist0 = dltlist0
		do i = 1, n1
!			refdvallist0(i)=sqrt(drholist0(1,i)**2.0+drholist0(2,i)**2.0+drholist0(3,i)**2.0)
			refdvallist(i)=sqrt(ddltlist0(1,i)**2.0+ddltlist0(2,i)**2.0+ddltlist0(3,i)**2.0)
		enddo	
	
		num_bins = 2
		allocate(muavlist(num_bins), muerlist(num_bins), ravlist(num_bins), muvarlist(num_bins))
		do idrop = 1, cs%numdrop
			! cycle if no drop
			if(.not.cs%dropval(idrop) .and. .not.cs%dropdval(idrop)) then
				! chisq of rho, dlt, ndlt
				call list_get_mu(poslist0, drholist0, drho_mu_data)
				call list_get_mu(poslist0, ddltlist0, ddlt_mu_data)
				call list_get_mu(poslist0, dndltlist0, dndlt_mu_data)
				rhochisqlist(idrop) = chisq_of_mu_data2(drho_mu_data)
				dltchisqlist(idrop) = chisq_of_mu_data2(ddlt_mu_data)
				ndltchisqlist(idrop) = chisq_of_mu_data2(dndlt_mu_data)
				deallocate(drho_mu_data, ddlt_mu_data, dndlt_mu_data)
			else
				markdrop = 0
				refvallist = refvallist0
				refdvallist =refdvallist0
				if(cs%dropval(idrop)) then
					call mark_drop_pixels(refvallist, markdrop, n1,&
						int(cs%lowdropvalratio(idrop)*n1)-1, int(n1*( 1.0_dl-cs%highdropvalratio(idrop) )) + 1)
				endif
				if(cs%dropdval(idrop)) then
					call mark_drop_pixels(refdvallist, markdrop, n1,&
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
				
				! Applying dropping & calculating chisq
				call drop_pixels2(poslist0, rholist0, drholist0, markdrop, poslist, rholist, drholist)
				call drop_pixels2(poslist0, dltlist0, ddltlist0, markdrop, poslist, dltlist, ddltlist)
				call drop_pixels2(poslist0, ndltlist0, dndltlist0, markdrop, poslist, ndltlist, dndltlist)
				call list_get_mu(poslist, drholist, drho_mu_data)
				call list_get_mu(poslist, ddltlist, ddlt_mu_data)
				call list_get_mu(poslist, dndltlist, dndlt_mu_data)
				rhochisqlist(idrop) = chisq_of_mu_data2(drho_mu_data)
				dltchisqlist(idrop) = chisq_of_mu_data2(ddlt_mu_data)
				ndltchisqlist(idrop) = chisq_of_mu_data2(dndlt_mu_data)
				
				! chisq from differantial of mu
				! Preparing for calculating df chisqs
				allocate(rlist(size(rholist)))
				do i = 1, size(rlist)
					rlist(i) = sqrt(poslist(1,i)**2.0+poslist(2,i)**2.0+poslist(3,i)**2.0)
				enddo
				call find_min_max(rlist, size(rlist), rmin, rmax)
								
				do i = 1, size(drho_mu_data)
					drho_mu_data(i) = abs(drho_mu_data(i))
					ddlt_mu_data(i) = abs(ddlt_mu_data(i))
					dndlt_mu_data(i) = abs(dndlt_mu_data(i))
				enddo
			  	call binned_quan(drho_mu_data, rlist, rmin, rmax, size(rlist), num_bins, &
			  		muavlist, muerlist, ravlist, muvarlist)
			  	rhodfchisqlist(idrop) = abs(muavlist(2) - muavlist(1))**2.0 / ((muerlist(1)**2.0 + (muerlist(2)**2.0)))
			  	call binned_quan(ddlt_mu_data, rlist, rmin, rmax, size(rlist), num_bins, &
			  		muavlist, muerlist, ravlist, muvarlist)
			  	dltdfchisqlist(idrop) = abs(muavlist(2) - muavlist(1))**2.0 / ((muerlist(1)**2.0 + (muerlist(2)**2.0)))
			  	call binned_quan(dndlt_mu_data, rlist, rmin, rmax, size(rlist), num_bins, &
			  		muavlist, muerlist, ravlist, muvarlist) 
			  	ndltdfchisqlist(idrop) = abs(muavlist(2) - muavlist(1))**2.0 / ((muerlist(1)**2.0 + (muerlist(2)**2.0)))
				
				deallocate(rlist,poslist,rholist,drholist,drho_mu_data,&
					dltlist,ddltlist,ddlt_mu_data,ndltlist,dndltlist,dndlt_mu_data)
!				print *, 'gd_mldprho_chi2s J'
			endif
		enddo
		call cpu_time(time4)		
		deallocate(markdrop,refvallist,refdvallist)

		if(print_time .or. gradfieldprintinfo) then
			print *, '  Time used in grid_rho_drho_list: ', real(time1-time0)
			print *, '  Time used in dlt list: ', real(time2-time1)			
			print *, '  Time used in ndlt list: ', real(time3-time2)	
			print *, '  Time used in drop: ', real(time4-time3)
		endif
		gradfieldprintinfo = .false.
!		print *, '(gd_mldprho_chi2s) Test E...'	
	end subroutine gd_mldprho_mlchi2s
	
	
  !------------------------------------------
  ! Reset the mass_list: normalized by avg;
  !------------------------------------------
	subroutine avg_mlist(distance_list, quan_list, n, cs, nsm)
		!Dummy
		integer :: n, nsm
		type(chisq_settings) :: cs
		real(dl) :: distance_list(n), quan_list(n), remov_dist_ratio
		!Local
		real(dl) :: remov_dist, boundary_rmin, boundary_rmax, delta_r, r, quan_av
		real(dl), allocatable :: quan_av_list(:), quan_er_list(:), binned_r_list(:)
		integer :: i_sm, i, i1, i2, i3, binned_i
	
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
  		
  		allocate(quan_av_list(cs%nbins_rhoav), quan_er_list(cs%nbins_rhoav), binned_r_list(cs%nbins_rhoav))
  		do i_sm = 1, nsm
  			call binned_quan(quan_list, distance_list, boundary_rmin, boundary_rmax, n, cs%nbins_rhoav, &
		  		quan_av_list, quan_er_list, binned_r_list)
	  		do i = 1, size(gb_r_list)
  				r = gb_r_list(i)
  				binned_i = max(min(int((r-boundary_rmin)/delta_r)+1,cs%nbins_rhoav),1)
  				if(cs%use_intpl_rho) then
  					call ilist(binned_i, 1, 1, cs%nbins_rhoav,i1,i3)
  					i2 = i1 + 1
  					quan_av = intpl_vl(r,binned_r_list(i1),quan_av_list(i1),binned_r_list(i2),quan_av_list(i2),binned_r_list(i3),quan_av_list(i3))
  					gb_mass_list(i) = gb_mass_list(i) / quan_av
  				else
  					gb_mass_list(i) = gb_mass_list(i) / quan_av_list(binned_i)
  				endif
  			enddo
  		enddo
  	end subroutine avg_mlist
  	

  !------------------------------------------
  ! Reset the mass_list: normalized by avg;
  !------------------------------------------
	subroutine sqvar_mlist(distance_list, quan_list, n, cs, nsm)
		!Dummy
		integer :: n, nsm
		type(chisq_settings) :: cs
		real(dl) :: distance_list(n), quan_list(n), remov_dist_ratio
		!Local
		real(dl) :: remov_dist, boundary_rmin, boundary_rmax, delta_r, r, dvar_sqrt
		real(dl), allocatable :: quan_av_list(:), quan_er_list(:), quan_var_list(:), binned_r_list(:)
		integer :: i_sm, i, i1, i2, i3, binned_i
	
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
  		
  		allocate(quan_av_list(cs%nbins_rhoav), quan_er_list(cs%nbins_rhoav), binned_r_list(cs%nbins_rhoav), &
  			quan_var_list(cs%nbins_rhoav))
  		do i_sm = 1, nsm
			call binned_quan(quan_list, distance_list, boundary_rmin, boundary_rmax, n, cs%nbins_rhoav, &
				quan_av_list, quan_er_list, binned_r_list, quan_var_list)

  			do i = 1, size(gb_r_list)
  				r = gb_r_list(i)
  				binned_i = max(min(int((r-boundary_rmin)/delta_r)+1,cs%nbins_rhoav),1)
  				if(cs%use_intpl_rho) then
  					call ilist(binned_i, 1, 1, cs%nbins_rhoav,i1,i3)
  					i2 = i1 + 1
  					dvar_sqrt = intpl_vl(r,binned_r_list(i1),sqrt(quan_var_list(i1)),binned_r_list(i2),&
  						sqrt(quan_var_list(i2)),binned_r_list(i3),sqrt(quan_var_list(i3)))
  					gb_mass_list(i) = gb_mass_list(i) / dvar_sqrt
  				else
  					gb_mass_list(i) = gb_mass_list(i) / sqrt(quan_var_list(binned_i))
  				endif
  			enddo
  		enddo
  	end subroutine sqvar_mlist


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
  	

  !------------------------------------------
  ! estimating gradient of fields of 
  !  rho/delta/normed delta at grid points
  !------------------------------------------
	subroutine gd_mdwei_chi2s(RSD, AP, cs, changenuminx, chisqlist, numchisq)
		! Dummy
		integer, intent(in) :: RSD, AP
		type(chisq_settings) :: cs
		integer :: numchisq
		real(dl), intent(in) :: changenuminx
		real(dl), intent(out) :: chisqlist(numchisq)
		! Local
		real(dl) :: boundary_rmin, boundary_rmax, rhomin, rhomean
		integer :: i, n
		real(dl), allocatable :: pos_list(:,:), rho_list(:), drho_list(:,:), absdrholist(:), weight(:)
		real(dl), allocatable :: drho_mu_data(:)

!		print *, 'gdmdwei A-2'

		if(cs%print_info) print *, 'Estimating rho, delta, normed_delta: changenuminx = ', real(changenuminx), '...'
		call grid_rho_drho_list(RSD, AP, cs, changenuminx, pos_list, rho_list, drho_list, &
				boundary_rmin, boundary_rmax)
!		print *, 'gdmdwei A-1'
		call list_get_mu(pos_list, drho_list, drho_mu_data)
!		print *, 'gdmdwei A'		
		n = size(rho_list)
		allocate(absdrholist(n),weight(n))
		do i = 1, n
			absdrholist(i)=sqrt(drho_list(1,i)**2.0+drho_list(2,i)**2.0+drho_list(3,i)**2.0)
		enddo

!		print *, 'gdmdwei B'
		! unweighted chisq
		chisqlist(1) = chisq_of_mu_data2(drho_mu_data)
		
		! weighted by 1/rho
		do i = 1, n
			weight(i) = 1.0_dl / rho_list(i)
		enddo
		
		chisqlist(2) = weighted_chisq(drho_mu_data, weight, n)
		
		! weighted by 1/sqrt(rho)
		do i = 1, n
			weight(i) = 1.0_dl / sqrt(rho_list(i))
		enddo
		chisqlist(3) = weighted_chisq(drho_mu_data, weight, n)
				
		! weighted by 1/(log(rho/rhomin)+1)
		rhomin = minval(rho_list)
		do i = 1, n
			weight(i) = 1.0_dl / (log(rho_list(i)/rhomin)+1.0)
		enddo
		chisqlist(4) = weighted_chisq(drho_mu_data, weight, n)
		
		! switch from log(<rho>/rho)+1.0 to 1/(log(rho/<rho>)+1) from rho less to larger than <rho>
		call get_mean_var(rho_list, rhomean)
		do i = 1, n
			if(rho_list(i) .le. rhomean) then
				weight(i) = log(rhomean/rho_list(i)) + 1.0_dl
			else				
				weight(i) = 1.0_dl / (log(rho_list(i)/rhomean)+1.0)
			endif
		enddo
		chisqlist(5) = weighted_chisq(drho_mu_data, weight, n)
		deallocate(pos_list,rho_list,drho_list,absdrholist,weight,drho_mu_data)	
	end subroutine gd_mdwei_chi2s

  	
  !------------------------------------------
  ! estimating gradient of fields of 
  !  rho/delta/normed delta at grid points:
  ! chisq from evloving effect
  !------------------------------------------
	subroutine gd_mldprho_dschi2s(RSD, AP, cs, changenuminx, chisqlist)
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
		!variables used to estimate mu at different r bins
		real(dl), allocatable :: r_list(:), muavlist(:), muerlist(:), ravlist(:), muvarlist(:)
		real(dl) :: rmin, rmax
		integer :: num_bins
		
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

		num_bins = 2
		allocate(muavlist(num_bins), muerlist(num_bins), ravlist(num_bins), muvarlist(num_bins))
!		print *, '(gd_mldprho_chi2s) Test D...'
		do idrop = 1, cs%numdrop
			! cycle if no drop
			if(.not.cs%dropval(idrop) .and. .not.cs%dropdval(idrop)) then
				call list_get_mu(pos_list0, drho_list0, drho_mu_data)
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
				call list_get_mu(pos_list, drho_list, drho_mu_data)
				chisqlist(idrop) = chisq_of_mu_data2(drho_mu_data)
!				print *, 'chisq = ', chisqlist(idrop)
				
				allocate(r_list(size(rho_list)))
				do i = 1, size(r_list)
					r_list(i) = sqrt(pos_list(1,i)**2.0+pos_list(2,i)**2.0+pos_list(3,i)**2.0)
				enddo
				call find_min_max(r_list, size(r_list), rmin, rmax)
								
				do i = 1, size(drho_mu_data)
					drho_mu_data(i) = abs(drho_mu_data(i))
				enddo
			  	call binned_quan(drho_mu_data, r_list, rmin, rmax, size(drho_mu_data), num_bins, &	
  					muavlist, muerlist, ravlist, muvarlist)
				if(gradfieldprintinfo) then
	  				print *
					print *, 'num_bins: ', num_bins
					print *, 'rmin/rmax: ', rmin, rmax
					print *, 'muavlist: ', muavlist
					print *, 'muerlist: ', muerlist
					print *, 'ravlist: ', ravlist
					print *
					print *, 'chisq of diff mu: ', chisqlist(idrop)
					print *
				endif
				chisqlist(idrop) = abs(muavlist(2) - muavlist(1))**2.0 / ((muerlist(1)**2.0 + (muerlist(2)**2.0)))
				
				deallocate(r_list,pos_list,rho_list,drho_list,drho_mu_data)
!				print *, 'gd_mldprho_chi2s J'
			endif
		enddo
				
		call cpu_time(time2)		

		deallocate(absdrholist,markdrop,reflist)

		if(print_time .or. gradfieldprintinfo) then
			print *, '  Time used in grid_rho_drho_list: ', real(time1-time0)
			print *, '  Time used in drop: ', real(time2-time1)
			print *
		endif
		gradfieldprintinfo = .false.
!		print *, '(gd_mldprho_chi2s) Test E...'		
	end subroutine gd_mldprho_dschi2s 	

end module ap_grad_fields	
	



