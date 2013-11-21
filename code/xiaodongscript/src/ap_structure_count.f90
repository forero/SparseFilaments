
! This module counts # of structures (e.g., >1.7 sigma peaks) to probe expansion history

module ap_structure_count

use ap_grad_fields

	implicit none

	type :: connect_mat_element
		integer :: binnum, crflag = -1 ! flag for cr (connected region)
		logical :: marked = .false., checked = .false.
	end type
	type(connect_mat_element), allocatable, public :: connect_mat(:,:,:)

	
contains
	
  !------------------------------------------
  ! initialization of stucount
  !  I think it may be better to input
  !   size rather than # of pixel 
  !------------------------------------------
	subroutine strucount(RSD, AP, cs, pixelsize)
		! DUMMY
		integer, intent(in) :: RSD, AP
		type(chisq_settings) :: cs
		real(dl), intent(in) :: pixelsize
		! Local
		real(dl) :: changenuminx, rl_num_in_x, boundary_rmin,boundary_rmax, rmin,rmax,deltar
		real(dl), allocatable :: pos_list(:,:), rho_list(:), drho_list(:,:), distance_list(:)
		real(dl), allocatable :: rho_av_list(:), rho_er_list(:), binned_r_list(:)
		integer :: n, i, ix,iy,iz, binnum, nowcrflag
		integer, allocatable :: ixiyizlist(:,:), binnedpixelnum(:), binnedflagnum(:)
		logical :: alreadydone
		! Testing
		real(dl) :: testx
		integer :: testnum

		if(cs%print_info) then
		        write(*,'(A,i4,i4,A)'), '  (strucount_init) Initializing grids of pixels with RSD, AP = ', RSD, AP, '  ...'
		endif
		
		call init_xyz_r_gb_mass_list(RSD, AP, cs%print_info)

		! fixed grid range. x,y,z start from 0 !!!NOTICE
		fixgridrange = .true.
		fixgridxmin = 0.0_dl;	fixgridxmax = gbxmax
		fixgridymin = 0.0_dl;	fixgridymax = gbymax
		fixgridzmin = 0.0_dl;	fixgridzmax = gbzmax
		rl_num_in_x = (fixgridxmax - fixgridxmin) / pixelsize
		changenuminx = rl_num_in_x / dble(cs%num_in_x) - 1.0

		if(cs%print_info) then
			if(fixgridrange) &
				print *, '  (strucount_init) Use fix grid range: For x we use ', &
					real(fixgridxmin), real(fixgridxmax)
			print *, '  (strucount_init) Use rl_num_in_x = ', rl_num_in_x
		endif
		
		! calculating rho at each grid point
		! notice that whether using num density has been included in cs
		call grid_rho_drho_list(RSD, AP, cs, changenuminx, pos_list, rho_list, drho_list, boundary_rmin, boundary_rmax,&
			ixiyizlist)

		! seperate the sample into different bins; find out the middle value of rho in each bin
		n = size(rho_list)
		allocate(distance_list(n))
		do i = 1, n
			distance_list(i) = dsqrt(pos_list(1,i)**2.0+pos_list(2,i)**2.0+pos_list(3,i)**2.0)
		enddo
		
		! removing evolving in mass list (divided by *middle* value of rho in each bin)
		call medrho_mlist(distance_list, rho_list, n, cs, 1)

		! list of rhos after removoing evolving effect
		call get_val_dval_list(cs%smnum, cs%print_info, pos_list, rho_list, drho_list, save_in_cellmat = .true.)
		call medrho_mlist(distance_list, rho_list, n, cs, 1) ! For test: coefficient shall be close to 1 or 0
				
		! <rho> at each bin
		allocate(rho_av_list(cs%nbins_rhoav), rho_er_list(cs%nbins_rhoav), binned_r_list(cs%nbins_rhoav))
		call find_min_max(distance_list, n, rmin, rmax)
		call binned_quan(rho_list, distance_list, rmin, rmax, n, cs%nbins_rhoav, &
		  		rho_av_list, rho_er_list, binned_r_list)
		
		!... We shall count # of connected regions in each bin
		if(allocated(connect_mat)) deallocate(connect_mat)
		allocate(connect_mat(gb_n_cellx,gb_n_celly,gb_n_cellz),binnedpixelnum(cs%nbins_rhoav))
		binnedpixelnum = 0
		deltar = (rmax - rmin) / dble(cs%nbins_rhoav)
  		do i = 1, n
	  		binnum = max(min(int((distance_list(i)-rmin)/deltar)+1,cs%nbins_rhoav),1)
	  		binnedpixelnum(binnum) = binnedpixelnum(binnum) + 1
			ix=ixiyizlist(1,i);iy=ixiyizlist(2,i);iz=ixiyizlist(3,i);
			connect_mat(ix,iy,iz)%binnum = binnum
			! structures: rho > 1.7 sigma
	  		if(rho_list(i) .ge. rho_av_list(binnum) + 1.7_dl*rho_er_list(binnum)) then
	  			connect_mat(ix,iy,iz)%marked = .true.
	  		endif
  		enddo
  		
  		nowcrflag = 1
  		allocate(binnedflagnum(cs%nbins_rhoav))
  		binnedflagnum = 0
  		do ix=1,gb_n_cellx
  		do iy=1,gb_n_celly
  		do iz=1,gb_n_cellz
  			if(connect_mat(ix,iy,iz)%marked.eq..false.) then
  				connect_mat(ix,iy,iz)%checked = .true.; cycle
  			elseif(connect_mat(ix,iy,iz)%checked.eq..true.) then
  				if(connect_mat(ix,iy,iz)%crflag .eq. -1) then
  					print *, 'Strang! Checked marked shall not has crflag==-1! '
  					print *, ix,iy,iz, connect_mat(ix,iy,iz)
  				endif
  				cycle
  			else
  				binnum = connect_mat(ix,iy,iz)%binnum 
	  			call conecregion(ix,iy,iz,binnum,nowcrflag)
  				nowcrflag = nowcrflag + 1
  				binnedflagnum(binnum) = binnedflagnum(binnum) + 1
  			endif
		enddo
		enddo
		enddo
		
		print *, 'Total # of cr: ', sum(binnedflagnum)
		do i =1,cs%nbins_rhoav
			print *, '#, ratio of cr in bin ', i,':', binnedflagnum(i), binnedflagnum(i) / real(binnedpixelnum(i))
		enddo
		
		! Testing the result by drawing ...
		open(unit=89071,file='Test/CRinfo.txt')
  		do ix=1,gb_n_cellx
  		do iy=1,gb_n_celly
  		do iz=1,gb_n_cellz
  			if(connect_mat(ix,iy,iz)%binnum.eq.4 .and. connect_mat(ix,iy,iz)%marked.eq..true.) then
  				write(89071,'(3(i3),i5)') ix,iy,iz,connect_mat(ix,iy,iz)%crflag
  			endif
  			if(connect_mat(ix,iy,iz)%marked .and. (connect_mat(ix,iy,iz)%crflag.eq.-1)) then
!  				print *, 'ERROR!!!! Marked pixels shall has crflag!'
 ! 				print *, ix,iy,iz, connect_mat(ix,iy,iz)
!  				stop
  			elseif((.not.connect_mat(ix,iy,iz)%marked) .and. (connect_mat(ix,iy,iz)%crflag.ne.-1)) then
				print *, 'ERROR!!!! Unmarked pixels shall not have crflag!'
 ! 				stop
  			endif
  		enddo
  		enddo
  		enddo
  		close(89071)
			  			
!		call medrho_mlist(distance_list, rho_list, n, cs, 1)
	end subroutine strucount


  !------------------------------------------
  ! recursive subroutine which searches for
  !  a connected region
  ! All marked pixels which belongs to a same
  !  connected region will be assigned crflag
  ! It does two things:
  !  (1). If it's not satisfied or already correctly assigned, return with alreadydone = .true.
  !  (2). If it's satisfied and not assigned; then:
  !       (a). Assign crflag to it;
  !       (b). Recursivly do the same procedure for nearby pixels
  !------------------------------------------
	recursive subroutine conecregion(ix,iy,iz,checkbinnum,crflag)
		! Dummy
		integer, intent(in) :: ix,iy,iz,crflag,checkbinnum
		!local
		integer :: i,j,k
		! return if already checked, not marked, or not in same bin
		if(connect_mat(ix,iy,iz)%checked .or. (.not.connect_mat(ix,iy,iz)%marked) &
 		   .or. connect_mat(ix,iy,iz)%binnum.ne.checkbinnum) then
			connect_mat(ix,iy,iz)%checked = .true.
			return
		endif

		! assign crflag to the pixel
!		write(*,'(A,4(i3),1x,i5)') 'Assigning crflag: ix,iy,iz,binnum,crflag = ',ix,iy,iz,nowbinnum,crflag
		if(connect_mat(ix,iy,iz)%crflag .ne. -1) then
			print *, '(conecregion) Strange! Shall be -1!', ix,iy,iz, connect_mat(ix,iy,iz)%crflag
		endif
		connect_mat(ix,iy,iz)%crflag = crflag
		
		! recursively check nearby pixels
		do i = max(1,ix-1), min(ix+1,gb_n_cellx)
		do j = max(1,iy-1), min(iy+1,gb_n_celly)
		do k = max(1,iz-1), min(iz+1,gb_n_cellz)
			call conecregion(i,j,k,checkbinnum,crflag)
		enddo
		enddo
		enddo
		return
	end subroutine conecregion
  !------------------------------------------
  ! Reset the mass_list: normalized by avg;
  !------------------------------------------
	subroutine medrho_mlist(distance_list, quan_list, num_quan, cs, nsm)
		!Dummy
		integer :: num_quan, nsm
		type(chisq_settings) :: cs
		real(dl) :: distance_list(num_quan), quan_list(num_quan), remov_dist_ratio
		!Local
		real(dl), allocatable :: quan_av_list(:), quan_er_list(:), binned_r_list(:), quan_var_list(:), quan_medval_list(:)
		integer, parameter :: polyorder = 3
		real(dl) :: r, rmin, rmax, polyrho, polycoef(polyorder+1)
		integer :: i_sm, i
		!Variables used in Test
!		character(len=char_len) :: outputfile
!		real(dl) :: x
		
		if(cs%print_info) then
			print *, '  (medrho_mlist) estimating medium rho for ', cs%nbins_rhoav, ' bins, distance range ', real(gbrmin), ' to ', real(gbrmax), '...'
		endif
		
		allocate(quan_av_list(cs%nbins_rhoav), quan_er_list(cs%nbins_rhoav), binned_r_list(cs%nbins_rhoav), &
			quan_var_list(cs%nbins_rhoav), quan_medval_list(cs%nbins_rhoav))
  		do i_sm = 1, nsm
  			call find_min_max(distance_list, num_quan, rmin, rmax)
  			call binned_quan(quan_list, distance_list, rmin, rmax, num_quan, cs%nbins_rhoav, &
		  		quan_av_list, quan_er_list, binned_r_list, quan_var_list, quan_medval_list)
!		  	print *, 'binned_r_list:   ', binned_r_list
!		  	print *, 'quan_medval_list:', quan_medval_list
		  	call poly_fit(binned_r_list,quan_medval_list,polycoef,cs%nbins_rhoav,polyorder)
		  	print *, 'polycoef: ', polycoef
!  			call random_seed()
!  			open(unit=872,file='Test/pixelrlist.txt')
!			open(unit=873,file='Test/pixelrho.txt')
!  			open(unit=874,file='Test/polyrlist.txt')
!			open(unit=875,file='Test/polyrho.txt')
	  		do i = 1, num_halo
  				r = gb_r_list(i)
  				polyrho = poly(r,polycoef,polyorder)
  				gb_mass_list(i) = gb_mass_list(i) / polyrho
!  				call random_number(x)
!  				if(x .le. 100.0 / dble(num_halo)) then
!					write(874,*) r
!  					write(875,*) polyrho
!  				endif
  			enddo
!  			do i = 1, num_quan
!  				call random_number(x)
!  				if(x .le. 3000.0 / dble(num_quan)) then
!  					write(872,*) distance_list(i)
!  					write(873,*) quan_list(i)
!  				endif
!  			enddo
!  			close(872);close(873);close(874);close(875);
  			! output for testing

  			
  		enddo
  	end subroutine medrho_mlist

end module ap_structure_count
