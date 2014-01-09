! Check 1: with the same cri why # of CRs increase? (it shall decrease if we donot split CRs)
! Check 2: why chisqmpA changes a little? (there is a one pixel difference)

! This module counts # of structures (e.g., >1.7 sigma peaks) to probe expansion history

module ap_structure_count

use ap_grad_fields

	implicit none

	! Settings
	
	logical, public :: strucount_checksaddle = .false. ! whether checking saddle points
	real(dl), public, parameter :: cutsigma = 1.7_dl, crifrac = 0.9554345372317531_dl !
	logical, allocatable :: pixelbdeffectinfo(:,:,:)
	
	type :: connect_mat_element
		integer :: crflag = -1 ! flag for cr (connected region)
		logical :: marked = .false., checked = .false.
	end type
	type(connect_mat_element), allocatable :: connect_mat(:,:,:)

	! Infomation of CR
	type :: cr
		integer :: numofpixel=0
		integer(2) :: binnum=-1
		real(dl) :: sumrho=0.0,meanx=0.0,meany=0.0,meanz=0.0
		logical :: has_bd_effect = .false., split_by_bin = .false.
	end type
	type(cr), allocatable :: crinfolist(:)

	real, allocatable :: saddle_rho(:,:,:)
	logical, allocatable :: check_saddle_rho(:,:,:)
	
contains

  !------------------------------------------
  ! stucount
  !------------------------------------------
  	subroutine strucount(RSD, AP, cs, pixelsize)
  		! Dummy
  		integer, intent(in) :: RSD, AP
  		type(chisq_settings), intent(in) :: cs
		real(dl), intent(in) :: pixelsize ! I think it may be better to input size rather than # of pixel 
		! Local
		integer :: n, coneclev
		real(dl), allocatable :: pos_list(:,:),rho_list(:),distance_list(:)
		integer, allocatable :: ixiyizlist(:,:)
		real(dl) :: rho_cri_list(cs%nbins_rhoav)
		real(dl) :: rhocricoef(nfitrhocri+1)
		logical, parameter :: savememory = .true.
		! TESTING
		integer :: ix,iy,iz, i,j,k, nowcr
		logical, allocatable :: written_already(:,:,:)
		
		if(cs%print_info) print *, '  (strucount) Initializing density field...'

	  	call strucount_init(RSD, AP, cs, pixelsize, pos_list,rho_list,distance_list,ixiyizlist,n)

		if(savememory) then
			if(cs%print_info) write(*,'(A)') '   (strucount) Deallocating gb_cell_mat to save memory...'
			deallocate(gb_cell_mat)
		endif

		do coneclev = 0,0
			print *, '##############################'
			if(cs%print_info) write(*,'(A,i2)') '   (strucount) structure based on connection level ', coneclev
		  	call strucount_conecmat(coneclev,cs,pos_list,rho_list,distance_list,ixiyizlist,n)
		enddo
		
		! TESTING: output and check the info of CRs
		open(unit=13,file='Test/CRinfo.txt')
		allocate(written_already(gb_n_cellx,gb_n_celly,gb_n_cellz))
		written_already = .false.
		do ix=1,gb_n_cellx
		do iy=1,gb_n_celly
		do iz=1,gb_n_cellz
!			if(written_already(ix,iy,iz)) cycle
			if(connect_mat(ix,iy,iz)%marked) then
				nowcr = connect_mat(ix,iy,iz)%crflag
				if(crinfolist(nowcr)%has_bd_effect) then
					write(13,'(5i6)') ix,iy,iz,nowcr,1
					written_already(ix,iy,iz) = .true.
					do i=max(ix-1,1),min(ix+1,gb_n_cellx)
  					do j=max(iy-1,1),min(iy+1,gb_n_celly)
  					do k=max(iz-1,1),min(iz+1,gb_n_cellz)
!						if(written_already(i,j,k)) cycle
  						if(abs(i-ix)+abs(j-iy)+abs(k-iz)>1) cycle ! only "direct nearby" pixel
  						if(pixelbdeffectinfo(i,j,k)) then
  							write(13,'(5i6)') i,j,k,-1,3
							written_already(i,j,k) = .true.
						endif
  					enddo
  					enddo
  					enddo
				else
					write(13,'(5i6)') ix,iy,iz,nowcr,0
					written_already(ix,iy,iz) = .true.
				endif
			endif
		enddo
		enddo
		enddo
		close(13)
	end subroutine strucount

	
  !------------------------------------------
  ! Initialization of stucount
  !  (calc gb_cell_mat, rho_list, drho_list, pos_list)
  !------------------------------------------
	subroutine strucount_init(RSD, AP, cs, pixelsize, pos_list,rho_list,distance_list,ixiyizlist,n)
		! DUMMY
		integer, intent(in) :: RSD, AP
		type(chisq_settings) :: cs
		real(dl), intent(in) :: pixelsize
		real(dl), allocatable, intent(out) :: pos_list(:,:), rho_list(:), distance_list(:)
		integer, allocatable, intent(out) :: ixiyizlist(:,:)
		integer, intent(out) :: n
		! Local
		real(dl) :: changenuminx, rl_num_in_x, boundary_rmin,boundary_rmax
		real(dl) :: x,y,z,rho,rsq, numpixelmean,numpixelvar
		real(dl), allocatable :: drho_list(:,:)
		integer :: numzerorho,i,ix,iy,iz
		! TESTING 
		real(dl) :: rmin, rmax, rho_av_list(nbinrhocri), rho_er_list(nbinrhocri), &
			binned_r_list(nbinrhocri), rho_var_list(nbinrhocri), rho_coef(nfitrhocri+1), realm
		
		if(cs%print_info) then
		        write(*,'(A,i4,i4)')  '   (strucount_init) Counting structures using RSD, AP = ', RSD, AP
		        write(*,'(A,f6.3,A)') '                    High-rho criteria: ', real(cutsigma), ' sigma.'
		endif

		call init_mult_lists(RSD, AP, cs%print_info)

		! fixed grid range. x,y,z start from 0 !!!NOTICE
		gridrangeset = use_fixgridrange
		fixgridxmin = gbxmin;	fixgridxmax = gbxmax
		fixgridymin = gbymin;	fixgridymax = gbymax
		fixgridzmin = gbzmin;	fixgridzmax = gbzmax
		rl_num_in_x = (fixgridxmax - fixgridxmin) / pixelsize
		changenuminx = rl_num_in_x / dble(cs%num_in_x) - 1.0

		if(cs%print_info) then
			if(gridrangeset.eq.use_fixgridrange) then
				write(*,'(A,f12.3,1x,f12.3)'), '   (strucount_init) Use fix grid range: xrange = ', &
					real(fixgridxmin), real(fixgridxmax)
			endif
			print *, '  (strucount_init) Use rl_num_in_x = ', real(rl_num_in_x)
		endif
		
		! calculating rho at each grid point; Notice that whether using num density has been included in cs
		call init_mult_lists(RSD, AP, cs%print_info)
		call do_cell_init(RSD, AP, real(cs%num_in_x)*(1.0_dl+changenuminx), cs%print_info)
		call grid_rho_drho_list(cs%smnum, cs%print_info, pos_list, rho_list, drho_list, ixiyizlist)
		! seperate the sample into different bins; find out the middle value of rho in each bin
		n = size(rho_list)
		allocate(distance_list(n))
		do i = 1, n
			distance_list(i) = sqrt(pos_list(1,i)**2.0+pos_list(2,i)**2.0+pos_list(3,i)**2.0)
		enddo
		
		! TESTING
!		print *, 'OUTPUTING RAW RHO LIST!'
		! List of rho after removing 0-rho
!		testfile = trim(adjustl(teststr))//'_raw_rholist.txt'
!		call output_1d(testfile,rho_list,n)
		print *, 'OUTPUTING CLEAN RHO LIST!'
		testfile = trim(adjustl(teststr))//'_clean_raw_rholist.txt'
		open(unit=1,file=testfile)
		do i = 1,n
			if(rho_list(i) .ne. 0.0_dl) then
				write(1,*) rho_list(i)
			endif
		enddo
		close(1)
		! list of r after removing 0-rho
		print *, 'OUTPUTING CLEAN DISTANCE LIST!'
		testfile = trim(adjustl(teststr))//'_clean_distancelist.txt'
		open(unit=1,file=testfile)
		do i = 1,n
			if(rho_list(i) .ne. 0.0_dl) then
				write(1,*) distance_list(i)
			endif
		enddo
		close(1)
		! binned values of mean/rms of rho after removing 
		call find_min_max(distance_list, n, rmin, rmax)
		call eqvl_binned_quan(rho_list, distance_list, rmin, rmax, n, nbinrhocri, &
		  		rho_av_list, rho_er_list, binned_r_list, rho_var_list)	
		print *, 'OUTPUTING BINNED MEAN / RMS OF RAW RHO!'
		testfile = trim(adjustl(teststr))//'rawrhoav.txt'
		open(unit=1,file=testfile)
		do i = 1,nbinrhocri
			realm = rho_var_list(i) / (rho_er_list(i))**2.0 + 1
			write(1,'(4(e14.7,1x))') binned_r_list(i), rho_av_list(i), sqrt(rho_var_list(i)), &
				sqrt(2.0/(realm-1.0))*rho_var_list(i)
		enddo
		close(1)
		call poly_fit(binned_r_list,rho_av_list,rho_coef,nbinrhocri,nfitrhocri)
		testfile4 = trim(adjustl(teststr))//'rawrhoav_polyfit.txt'
		open(unit=1,file=testfile4)
		do i = 1, nfitrhocri
			testx = rmin + (rmax-rmin)*dble(i-1)/(dble(nfitrhocri-1))
			testy = poly(testx,rho_coef,nfitrhocri)
			write(1,*) testx, testy
		enddo
		close(1)

		! END TESTING
		! How many pixels have zero rho
		numzerorho = 0
		do i = 1, n
			if(rho_list(i).eq. 0.0_dl) numzerorho = numzerorho + 1
		enddo
		if(cs%print_info) print *, '  (strucount_init) # / ratio of 0-rho pixel = ', numzerorho, real(numzerorho)/real(n)
		
		! list of rhos after removoing evolving effect
		call medrho_mlist(distance_list, rho_list, n, cs%print_info, nbinrhocri, 1)
		call get_val_dval_list(cs%smnum, cs%print_info, pos_list, rho_list, drho_list, n, save_in_cellmat = .true.)
!		call medrho_mlist(distance_list, rho_list, n, cs%print_info, 5, 1) ! For test: coefficient shall be close to 1 or 0

	  	! CCC bdeffect info
  		if(allocated(pixelbdeffectinfo)) deallocate(pixelbdeffectinfo)
		allocate(pixelbdeffectinfo(gb_n_cellx,gb_n_celly,gb_n_cellz))
		do ix=1,gb_n_cellx
		do iy=1,gb_n_celly
		do iz=1,gb_n_cellz
 			pixelbdeffectinfo(ix,iy,iz) = gb_cell_mat(ix,iy,iz)%rhodrho_has_er
  		enddo
  		enddo
  		enddo
	  end subroutine strucount_init

  !------------------------------------------
  ! stucount (conecmat part)
  !------------------------------------------
	  subroutine strucount_conecmat(coneclev, cs, pos_list,rho_list,distance_list,ixiyizlist,n,&
	  	opt_chisqcrA,opt_chisqmpA,opt_binnedcrratio,opt_binnedcrratioer,opt_binnedmpratio,opt_binnedmpratioer)
	  	! DUMMY
		integer, intent(in) :: coneclev, n
		type(chisq_settings), intent(in) :: cs
		real(dl), intent(in) :: pos_list(3,n), rho_list(n), distance_list(n)
		integer, intent(in) :: ixiyizlist(3,n)
		real(dl), intent(out), optional :: opt_chisqcrA,opt_chisqmpA
		real(dl), intent(out), optional :: opt_binnedcrratio(cs%nbins_rhoav), opt_binnedcrratioer(cs%nbins_rhoav), &
						   opt_binnedmpratio(cs%nbins_rhoav), opt_binnedmpratioer(cs%nbins_rhoav)
		! Local
		logical, parameter :: cr_bdsplit_weight = .true.
		real(dl) :: rmin,rmax,rmincube,totvol,deltavol,nowrcube
		real(dl) :: x,y,z,rho,rhocri,rsq, numpixelmean,numpixelvar
		real(dl) :: meancrratio,chisqcrA,chisqcrB, meanmpratio,chisqmpA
		real(dl) :: time1, time2
		real(dl) :: binnedcrratio(cs%nbins_rhoav), binnedcrratioer(cs%nbins_rhoav), &
			binnedmpratio(cs%nbins_rhoav), binnedmpratioer(cs%nbins_rhoav)
		real(dl) :: binnedpixelnum(cs%nbins_rhoav), binnedmpnum(cs%nbins_rhoav)
		real(dl) :: binned_cr_num(cs%nbins_rhoav), binned_bdcr_weinum(cs%nbins_rhoav) ! # of cr may be digit if we consider weights.
		integer :: ix,iy,iz, i, binnum, nowcrflag,nowcr, numcr, nummarked, numzerorho, sumdiff
		! In each bin check # of CR which: has bd effect; split by bin division; either has bd or be splitted (logical or)
		! Testing
		integer, parameter :: nfit = 100
		real(dl) :: rhocricoef(cs%nbins_rhoav), realm
		real(dl) :: rho_av_list(nbinrhocri), rho_er_list(nbinrhocri), &
			binned_r_list(nbinrhocri), rho_var_list(nbinrhocri), rho_cri_list(nbinrhocri)
		real(dl) :: rhocri_coef(nfitrhocri+1), rho_coef(nfitrhocri+1)
		integer :: testi,j,k

! criteria of rho at each bin
		call find_min_max(distance_list, n, rmin, rmax)
		rmincube=rmin**3.0; totvol=rmax**3.0-rmincube; deltavol=totvol/dble(cs%nbins_rhoav);
		
		! fit using more points
		call eqvl_binned_quan(rho_list, distance_list, rmin, rmax, n, nbinrhocri, &
		  		rho_av_list, rho_er_list, binned_r_list, rho_var_list)	
		call eqvl_fracrilist(crifrac, rho_list, distance_list, n, rmin, rmax, nbinrhocri, rho_cri_list, cs%print_info)
		call poly_fit(binned_r_list,rho_cri_list,rhocri_coef,nbinrhocri,nfitrhocri)
		call poly_fit(binned_r_list,rho_av_list,rho_coef,nbinrhocri,nfitrhocri)
		
		print *, 'OUTPUTING INFORMATION OF RHO_CRI & <RHO>!'
		testfile1 = trim(adjustl(teststr))//'rhocri.txt'
		open(unit=1,file=testfile1)
		do i = 1, nbinrhocri
			write(1,*) binned_r_list(i), rho_cri_list(i)
		enddo
		close(1)
		testfile2 = trim(adjustl(teststr))//'rhocri_polyfit.txt'
		open(unit=1,file=testfile2)
		do i = 1, nfit
			testx = rmin + (rmax-rmin)*dble(i-1)/(dble(nfit-1))
			testy = poly(testx,rhocri_coef,nfitrhocri)
			write(1,*) testx, testy
		enddo
		close(1)
		
		! TESTING
		print *, 'OUTPUTING RAW/CLEAN OF NORMALIZED RHO LIST!'
		! List of rho after removing 0-rho
		testfile = trim(adjustl(teststr))//'_norm_rholist.txt'
		call output_1d(testfile,rho_list,n)
		testfile = trim(adjustl(teststr))//'_clean_norm_rholist.txt'
		open(unit=1,file=testfile)
		do i = 1,n
			if(rho_list(i) .ne. 0.0_dl) then
				write(1,*) rho_list(i)
			endif
		enddo
		close(1)
		
		!CCCCC
		! evolving of <rho>
		print *, 'OUTPUTING BINNED MEAN / VAR OF RHO!'
		testfile3 = trim(adjustl(teststr))//'rhoav.txt'
		open(unit=1,file=testfile3)
		do i = 1, nbinrhocri
			realm = rho_var_list(i) / (rho_er_list(i))**2.0 + 1
			write(1,'(4(e14.7,1x))') binned_r_list(i), rho_av_list(i), sqrt(rho_var_list(i)), &
				sqrt(2.0/(realm-1.0))*rho_var_list(i)
		enddo
		close(1)
		testfile4 = trim(adjustl(teststr))//'rhoav_polyfit.txt'
		open(unit=1,file=testfile4)
		do i = 1, nfit
			testx = rmin + (rmax-rmin)*dble(i-1)/(dble(nfit-1))
			testy = poly(testx,rho_coef,nfitrhocri)
			write(1,*) testx, testy
		enddo
		close(1)
		

		!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC

		write(*,'(24x,A,<nbinrhocri>(f9.3))') 'mean of rho   : ', real(rho_av_list)
		write(*,'(24x,A,<nbinrhocri>(f9.3))') 'rms  of rho   : ', sqrt(real(rho_var_list))
		
		if(cs%print_info) then
			write(*,'(3x,A,2f10.4,2e15.7)') '(strucount_conecmat) Binning sets: rmin/rmax/totvol/deltavol = ',&
				rmin,rmax,totvol,deltavol
			write(*,'(22x,A,i2,A,i2,A)') 'Fractional criterias of rho based on', nfitrhocri, &
				'th-order poly-fitted ', nbinrhocri, ' equal-volume bins: '
			write(*,'(24x,A,<nbinrhocri>f9.2)')  'binned-distance: ', binned_r_list
			write(*,'(24x,A,<nbinrhocri>(f9.2))') 'mean of rho   : ', real(rho_av_list)
			write(*,'(24x,A,<nbinrhocri>(f9.2))') 'rms  of rho   : ', sqrt(real(rho_var_list))
			write(*,'(24x,A,<nbinrhocri>(f9.2))') 'mean+sigma*rms: ', real(rho_av_list)+cutsigma*sqrt(real(rho_var_list))
			write(*,'(24x,A,<nbinrhocri>f9.2)')   'rho-cris      : ', rho_cri_list
			write(*,'(24x,A,<nfitrhocri+1>f9.2)') 'Coeffs        : ', rhocri_coef
		endif

		! CCCCC saddle points. used to check confused nearby pixels. (Not finished...)
		if(strucount_checksaddle) then
			print *, 'ERROR (strucount)!!! saddle point routines not finished yet.'; stop
			if(cs%print_info) print *, '  (strucount_init) Calculating rho at saddle points...'
			call calc_saddle_mat(rmin,rmax,cs%nbins_rhoav,rho_cri_list,cs%print_info)
		else
			if(cs%print_info) print *, '  (strucount_init) Do not check saddle points...'
	  	endif
	
		!CCCCC
		! Count # of connected regions in each bin
		! Initialize connect mat: mark selected pixels; count # of pixels / marked pixels in each bin
		if(cs%print_info) print *, '  (strucount_conecmat) Constructing connect matrix...'
  		if(allocated(connect_mat)) deallocate(connect_mat); allocate(connect_mat(gb_n_cellx,gb_n_celly,gb_n_cellz))
		binnedpixelnum = 0
		binnedmpnum = 0
		nummarked = 0
  		do i = 1, n
	  		binnum = max(min(int((distance_list(i)**3.0-rmincube)/deltavol)+1,cs%nbins_rhoav),1)
	  		binnedpixelnum(binnum) = binnedpixelnum(binnum) + 1
	  		rhocri = poly(distance_list(i),rhocri_coef,nfitrhocri)
	  		if(rho_list(i) .ge. rhocri) then
!	  		if(rho_list(i) .ge. rho_cri_list(binnum)) then
				nummarked = nummarked+1
				binnedmpnum(binnum) = binnedmpnum(binnum) + 1
				ix=ixiyizlist(1,i);iy=ixiyizlist(2,i);iz=ixiyizlist(3,i);
	  			connect_mat(ix,iy,iz)%marked = .true.
	  		endif
  		enddo

  		! Search for cr regions; count # of cr in each bin (binned_cr_num)
		if(cs%print_info) print *, '  (strucount_conecmat) Searching for connected regions (CR)...'
		call conmat_cr_scan(coneclev,numcr)
		
		! Info of crs: estimating totrho, position of each cr
		if(cs%print_info) print *, '  (strucount_conecmat) Analyzing info of CRs: # = ', numcr
		call strucount_crinfo(cs,pos_list,rho_list,distance_list,ixiyizlist,rho_cri_list,rhocri_coef,n,numcr,binned_cr_num,binned_bdcr_weinum)
  					
		! chisqcrA based on ratio of cr regions
		if(cs%print_info) write(*,'(A)')' #######################################################################################'
		write(*,'(A,i3,A,f11.2,A)') '   (strucount_conecmat) Results of structount (level=',coneclev,', #=',sum(binned_cr_num),'):'
		write(*,'(32x,A,f11.2)') 'binnum   #-of-CR   realative ratio of CR (divided by pixel # in this bin)'
		do i =1,cs%nbins_rhoav
			binnedcrratio(i) = binned_cr_num(i) / real(binnedpixelnum(i))
			binnedcrratioer(i) = sqrt(1.0/dble(binned_cr_num(i)) + 1.0/dble(binnedpixelnum(i))) * binnedcrratio(i)
			write(*,'(34x,i4,f10.2,10x,f8.5,A,f9.6)') i,binned_cr_num(i),real(binnedcrratio(i)),' +/-',real(binnedcrratioer(i))
		enddo
		meancrratio = sum(binnedcrratio)/dble(cs%nbins_rhoav); chisqcrA=0;
		do i =1,cs%nbins_rhoav
			chisqcrA = chisqcrA + ((binnedcrratio(i)-meancrratio)/binnedcrratioer(i))**2.0
		enddo
		write(*,'(24x,A,f16.5)') 'chisq of CR-ratio-evolving: = ', chisqcrA
		
		! weighted cr #: boundary effect reduction
		do i=1, cs%nbins_rhoav
			binned_cr_num(i)=binned_cr_num(i)-binned_bdcr_weinum(i)!0.5*real(binned_bdcr_num(i))
		enddo
		write(*,'(A,i3,A,f11.2,A)') '   (strucount_conecmat) Results after boundary correction (level=',coneclev,&
			', #=',sum(binned_cr_num),'):'
		write(*,'(32x,A,f11.2)') 'binnum   #-of-CR   realative ratio of CR (divided by pixel # in this bin)'
		do i =1,cs%nbins_rhoav
			binnedcrratio(i) = binned_cr_num(i) / real(binnedpixelnum(i))
			binnedcrratioer(i) = sqrt(1.0/dble(binned_cr_num(i)) + 1.0/dble(binnedpixelnum(i))) * binnedcrratio(i)
			write(*,'(34x,i4,f10.2,10x,f8.5,A,f9.6)') i,binned_cr_num(i),real(binnedcrratio(i)),' +/-',real(binnedcrratioer(i))
		enddo
		meancrratio = sum(binnedcrratio) / dble(cs%nbins_rhoav)
		chisqcrB = 0
		do i =1,cs%nbins_rhoav
			chisqcrB = chisqcrB + ((binnedcrratio(i)-meancrratio)/binnedcrratioer(i))**2.0
		enddo
		write(*,'(24x,A,f16.8)') 'chisq of CR-ratio-evolving = ', chisqcrB

		! chisqsmpA based on ratio of marked pixels
		! # / ratio (divided by # of pixels) of mp (marked pixels) in each bin
		write(*,'(A,i10,A)') '   (strucount_conecmat) Results of mp (marked pixels) (#=',int(sum(binnedmpnum)+0.1),'):'
		write(*,'(32x,A,f11.2)') 'binnum   #-of-mp   realative ratio of mp (divided by pixel # in this bin)'
		do i =1,cs%nbins_rhoav
			binnedmpratio(i) = binnedmpnum(i) / real(binnedpixelnum(i))
			binnedmpratioer(i) = sqrt(1.0/dble(binnedmpnum(i)) + 1.0/dble(binnedpixelnum(i))) * binnedmpratio(i)
			write(*,'(34x,i4,i10,10x,f8.5,A,f9.6)') i,int(binnedmpnum(i)+0.1),real(binnedmpratio(i)),' +/-',real(binnedmpratioer(i))
		enddo
		meanmpratio = sum(binnedmpratio) / dble(cs%nbins_rhoav)
		chisqmpA = 0
		do i =1,cs%nbins_rhoav
			chisqmpA = chisqmpA + ((binnedmpratio(i)-meanmpratio)/binnedmpratioer(i))**2.0
		enddo
		write(*,'(24x,A,f16.8)') 'chisq of mp-ratio-evolving = ', chisqmpA
		if(cs%print_info) write(*,'(A)')' #######################################################################################'
		
		!!CCC
		if(present(opt_chisqcrA)) opt_chisqcrA=chisqcrA
		if(present(opt_chisqmpA)) opt_chisqmpA=chisqmpA
		if(present(opt_binnedcrratio)) opt_binnedcrratio = binnedcrratio
		if(present(opt_binnedcrratioer)) opt_binnedcrratioer = binnedcrratioer
		if(present(opt_binnedmpratio)) opt_binnedmpratio = binnedmpratio
		if(present(opt_binnedmpratioer)) opt_binnedmpratioer = binnedmpratioer

  		if(cs%print_info) print *, '  (strucount_conecmat) Checking consistencies in connect_mat...'
		call strucount_check(cs,coneclev)

  		if(cs%print_info) print *, '  (strucount_conecmat) Done.'
  	end subroutine strucount_conecmat

  !------------------------------------------
  ! Analyzing info of CR
  !------------------------------------------
	subroutine strucount_crinfo(cs,pos_list,rho_list,distance_list,ixiyizlist,rho_cri_list,rhocri_coef,n,numcr,&
		binned_cr_num,binned_bdcr_weinum)
		! Dummy
		type(chisq_settings), intent(in) :: cs
		real(dl), intent(in) :: pos_list(3,n),rho_list(n),distance_list(n),rhocri_coef(nfitrhocri+1),rho_cri_list(cs%nbins_rhoav)
		integer, intent(in) :: ixiyizlist(3,n), n, numcr
		real(dl), intent(out) :: binned_cr_num(cs%nbins_rhoav), binned_bdcr_weinum(cs%nbins_rhoav)
		! Local
		integer :: i,j,k, ix,iy,iz, nowcr,binnum
		real(dl) :: rhocri, binned_bdcr_num(cs%nbins_rhoav), binned_nobdcr_num(cs%nbins_rhoav), &
			binnedcrsplitnum(cs%nbins_rhoav), binnedcrbdsplitnum(cs%nbins_rhoav)
		real(dl) :: binned_nobdcr_size(cs%nbins_rhoav), binned_bdcr_size(cs%nbins_rhoav)
		real(dl) :: rmin,rmax,rmincube,totvol,deltavol, numpixelmean,numpixelvar, x,y,z, nowrcube

		if(cs%print_info) write(*,'(A,i9)') '   (strucount_crinfo) Initializing crinfolist: len = ', numcr
		call find_min_max(distance_list, n, rmin, rmax)
		rmincube=rmin**3.0; totvol=rmax**3.0-rmincube; deltavol=totvol/dble(cs%nbins_rhoav);
		if(allocated(crinfolist)) deallocate(crinfolist); allocate(crinfolist(numcr))
		do i = 1, n
			binnum = max(min(int((distance_list(i)**3.0-rmincube)/deltavol)+1,cs%nbins_rhoav),1)
	  		rhocri = poly(distance_list(i),rhocri_coef,nfitrhocri)
	  		if(rho_list(i) .ge. rhocri) then
!	  		if(rho_list(i) .ge. rho_cri_list(binnum)) then
				ix=ixiyizlist(1,i);iy=ixiyizlist(2,i);iz=ixiyizlist(3,i);				
	  			nowcr = connect_mat(ix,iy,iz)%crflag
	  			if(nowcr.le.0.or.nowcr.gt.numcr) then
	  				print *, 'ERROR (strucount_crinfo)!!! Strange out-of-range mark: ',nowcr,numcr; stop;
	  			endif
	  			crinfolist(nowcr)%numofpixel = crinfolist(nowcr)%numofpixel+1
	  			crinfolist(nowcr)%sumrho = crinfolist(nowcr)%sumrho + rho_list(i)
	  			crinfolist(nowcr)%meanx = crinfolist(nowcr)%meanx + rho_list(i)*pos_list(1,i)
	  			crinfolist(nowcr)%meany = crinfolist(nowcr)%meany + rho_list(i)*pos_list(2,i)
	  			crinfolist(nowcr)%meanz = crinfolist(nowcr)%meanz + rho_list(i)*pos_list(3,i)
	  		endif
  		enddo
  		binned_cr_num=0;numpixelmean = 0.0;numpixelvar = 0.0;
  		do nowcr = 1, numcr
	  		numpixelmean = numpixelmean + dble(crinfolist(nowcr)%numofpixel)
	  		numpixelvar = numpixelvar + dble(crinfolist(nowcr)%numofpixel**2.0)
  			crinfolist(nowcr)%meanx = crinfolist(nowcr)%meanx / crinfolist(nowcr)%sumrho;
  			crinfolist(nowcr)%meany = crinfolist(nowcr)%meany / crinfolist(nowcr)%sumrho;
  			crinfolist(nowcr)%meanz = crinfolist(nowcr)%meanz / crinfolist(nowcr)%sumrho;
  			nowrcube=sqrt(crinfolist(nowcr)%meanx**2.0d0+crinfolist(nowcr)%meany**2.0d0+crinfolist(nowcr)%meanz**2.0d0)**3.0
!  			write(*,'(A,i7,5(e14.7,1x))'), 'TESTING: nowcr, meanxyz, rcube,r = ',nowcr,crinfolist(nowcr)%meanx,&
 ! 				crinfolist(nowcr)%meany,crinfolist(nowcr)%meanz,nowrcube,nowrcube**(1.0/3.0)
			binnum = max(min(int((nowrcube-rmincube)/deltavol)+1,cs%nbins_rhoav),1)
			crinfolist(nowcr)%binnum = binnum
			binned_cr_num(binnum) = binned_cr_num(binnum) + 1
		enddo
		numpixelmean = numpixelmean/dble(numcr)
		numpixelvar = numpixelvar/dble(numcr) - numpixelmean**2.0
		if(cs%print_info) write(*,'(24x,A,f10.3,f10.3)') 'mean / sqrt(var) of # of pixels in cr = ', numpixelmean, sqrt(numpixelvar)
		
		! Check the boundary effect, splitted-into-bins effect in CRs;
		!  Get weighted # of cr (pixel has bd punished by weight = 0.5)
		if(cs%print_info) write(*,'(A)') '   (strucount_crinfo) Checking boundary and split-into-bins effect in CRs...'
		do ix=1,gb_n_cellx
		do iy=1,gb_n_celly
		do iz=1,gb_n_cellz
			if(.not.connect_mat(ix,iy,iz)%marked) then
				continue
			else
				nowcr = connect_mat(ix,iy,iz)%crflag
  				if(ix.eq.1.or.ix.eq.gb_n_cellx.or.iy.eq.1.or.iy.eq.gb_n_celly.or.iz.eq.1.or.iz.eq.gb_n_cellz) then
  					crinfolist(nowcr)%has_bd_effect = .true.
				else
					do i=max(ix-1,1),min(ix+1,gb_n_cellx)
  					do j=max(iy-1,1),min(iy+1,gb_n_celly)
  					do k=max(iz-1,1),min(iz+1,gb_n_cellz)
  						if(abs(i-ix)+abs(j-iy)+abs(k-iz)>1) cycle ! only "direct nearby" pixel
  						if(pixelbdeffectinfo(i,j,k)) crinfolist(nowcr)%has_bd_effect = .true.
  					enddo
  					enddo
  					enddo
  				endif
  				call cell_pos(ix,iy,iz,x,y,z)
  				nowrcube=(x**2.0+y**2.0+z**2.0)**(1.5_dl)
  				binnum = max(min(int((nowrcube-rmincube)/deltavol)+1,cs%nbins_rhoav),1)
  				if(binnum.ne.int(crinfolist(nowcr)%binnum)) then
  					crinfolist(nowcr)%split_by_bin = .true.
  				endif
			endif
		enddo
		enddo
		enddo
		binned_bdcr_num=0.0;binned_nobdcr_num=0.0;binnedcrsplitnum=0.0; binnedcrbdsplitnum=0.0; 
		binned_nobdcr_size=0.0;binned_bdcr_size=0.0;binned_bdcr_weinum=0.0;
  		do nowcr = 1, numcr
  			binnum = crinfolist(nowcr)%binnum
  			if(crinfolist(nowcr)%has_bd_effect) then
  				binned_bdcr_num(binnum)=binned_bdcr_num(binnum)+1
  				binned_bdcr_size(binnum)=binned_bdcr_size(binnum)+crinfolist(nowcr)%numofpixel
  			else
  				binned_nobdcr_num(binnum)=binned_nobdcr_num(binnum)+1
  				binned_nobdcr_size(binnum)=binned_nobdcr_size(binnum)+crinfolist(nowcr)%numofpixel
  			endif
  			if(crinfolist(nowcr)%split_by_bin)  binnedcrsplitnum(binnum)=binnedcrsplitnum(binnum)+1
  			if(crinfolist(nowcr)%has_bd_effect.or.crinfolist(nowcr)%split_by_bin)  &
  				binnedcrbdsplitnum(binnum)=binnedcrbdsplitnum(binnum)+1
  		enddo
  		binned_bdcr_size = binned_bdcr_size / dble(binned_bdcr_num)
  		binned_nobdcr_size = binned_nobdcr_size / dble(binned_nobdcr_num)

!		binned_bdcr_weinum = dble(binned_bdcr_num)*binned_bdcr_size/binned_nobdcr_size
		binned_bdcr_weinum = 0.0_dl
  		do nowcr = 1, numcr
  			binnum = crinfolist(nowcr)%binnum
  			if(crinfolist(nowcr)%has_bd_effect) then
  				! do boundary-correction only for small pixels
  				if(crinfolist(nowcr)%numofpixel < binned_nobdcr_size(binnum)) then
  					binned_bdcr_weinum(binnum) = binned_bdcr_weinum(binnum) + (1-crinfolist(nowcr)%numofpixel/binned_nobdcr_size(binnum))
  				endif
  			endif
		enddo		
		
		write(*,'(A)') '   (strucount_crinfo) Boundary effects of CRs at differnt bins:'
		write(*,'(26x,A)') 'binnum   #-of-cr   #-has_bd  <#-pixel>-of-has_bd  <#-pixel>-of-no_bd    correction'
		do i=1 ,cs%nbins_rhoav
			write(*,'(24x,i6,2i10,3f19.2)') i,int(binned_cr_num(i)+0.1),int(binned_bdcr_num(i)+0.1),binned_bdcr_size(i),binned_nobdcr_size(i),binned_bdcr_weinum(i)
		enddo		
	end subroutine strucount_crinfo

  !------------------------------------------
  ! Check of strucount result
  !------------------------------------------
	subroutine strucount_check(cs,coneclev)
		! Dummy
		type(chisq_settings), intent(in) :: cs
		integer, intent(in) :: coneclev
		! Local
		integer :: ix,iy,iz, i,j,k, sumdiff
		if(cs%print_info) then
			write(*,'(1x,A)') '  (strucount_check) Checking connect matrix:'
			write(*,'(1x,A)') '                      All marked pixels shall has assigned crflag;'
			write(*,'(1x,A)') '                      All unmarked pixels shall has no assigned crflag;'
			write(*,'(1x,A)') '                      Nearby, in same bin marked pixels shall be in same CR.'
		endif
  		do ix=1,gb_n_cellx
  		do iy=1,gb_n_celly
  		do iz=1,gb_n_cellz
  			if(.not.connect_mat(ix,iy,iz)%marked) continue
  			if(connect_mat(ix,iy,iz)%marked .and. (connect_mat(ix,iy,iz)%crflag.eq.-1)) then
  				print *, 'ERROR (strucount_conecmat)!!!! Marked pixels shall has crflag!'
				print *, 'ix,iy,iz = ',ix,iy,iz
  				print *, 'binnum,crflag,marked,checked = ', connect_mat(ix,iy,iz)
  				stop
  			elseif((.not.connect_mat(ix,iy,iz)%marked) .and. (connect_mat(ix,iy,iz)%crflag.ne.-1)) then
				print *, 'ERROR (strucount_conecmat)!!!! Unmarked pixels shall not have crflag!'
  				stop
  			endif
  			if(connect_mat(ix,iy,iz)%marked) then
	  			do i = max(1,ix-1), min(ix+1,gb_n_cellx)
				do j = max(1,iy-1), min(iy+1,gb_n_celly)
				do k = max(1,iz-1), min(iz+1,gb_n_cellz)
					sumdiff = (abs(i-ix)+abs(j-iy)+abs(k-iz))
					if(coneclev.eq.0) then
						goto 7512
					elseif(coneclev.eq.1.and.sumdiff.le.2) then
						goto 7512
					elseif(coneclev.eq.2.and.sumdiff.le.1) then
						goto 7512
					else
						cycle
					endif
7512					continue
					if(connect_mat(i,j,k)%marked.and.connect_mat(ix,iy,iz)%crflag.ne.connect_mat(i,j,k)%crflag) then
						print *, 'ERROR (strucount_conecmat)!!! Nearby same-bin pixels have different crflag: ', &
							connect_mat(ix,iy,iz)%crflag,connect_mat(i,j,k)%crflag
						print *, 'Pixel A: ', i,j,k,connect_mat(i,j,k)
						print *, 'Pixel B: ', ix,iy,iz,connect_mat(ix,iy,iz)
					endif
				enddo
				enddo
				enddo
			endif
  		enddo
  		enddo
  		enddo	
  	end subroutine strucount_check
  


  !------------------------------------------
  ! recursive subroutine which searches for
  !  a connected region
  ! Marked, nearby pixels will be assigned same crflag
  ! It does two things:
  !  (1). If it's not satisfied or already correctly assigned, return
  !  (2). If it's satisfied and not assigned; then:
  !       (a). Assign crflag to it;
  !       (b). Recursivly do the same procedure for nearby pixels
  !------------------------------------------
	recursive subroutine conecregion(ix,iy,iz,crflag,coneclev)
		! Dummy
		integer, intent(in) :: ix,iy,iz,coneclev
		integer, intent(in) :: crflag
		!local
		integer :: i,j,k,sumdiff

		! return if already checked or not marked
		if(connect_mat(ix,iy,iz)%checked .or. .not.connect_mat(ix,iy,iz)%marked)  return

		connect_mat(ix,iy,iz)%checked = .true.
		if(connect_mat(ix,iy,iz)%crflag .ne. -1) then
			write(*,*) '  (conecregion) Strange! Crflag for unchecked shall be -1!', connect_mat(ix,iy,iz)%crflag
			stop
		endif
		connect_mat(ix,iy,iz)%crflag = crflag

!		Check most nearby pixels
		do i = max(1,ix-1), min(ix+1,gb_n_cellx)
		do j = max(1,iy-1), min(iy+1,gb_n_celly)
		do k = max(1,iz-1), min(iz+1,gb_n_cellz)
			sumdiff = (abs(i-ix)+abs(j-iy)+abs(k-iz))
			if(coneclev.eq.0) then
				call conecregion(i,j,k,crflag,coneclev)
			elseif(coneclev.eq.1) then
				if(sumdiff .le. 2) call conecregion(i,j,k,crflag,coneclev)
			elseif(coneclev.eq.2) then
				if(sumdiff .le. 1) call conecregion(i,j,k,crflag,coneclev)
			else
				print *, 'ERROR (conecregion)!!! coneclev must be one of 0,1,2: ', coneclev
				stop
			endif
		enddo
		enddo
		enddo
	end subroutine conecregion
	
  !------------------------------------------
  ! Scan connect_mat and assigning crflags
  !------------------------------------------
	subroutine conmat_cr_scan(coneclev,numcr)
		! Dummy 
		integer, intent(in) :: coneclev
		integer, intent(out) :: numcr
		! Local
		integer :: nowcrflag, ix,iy,iz
  		nowcrflag = 0
  		do ix=1,gb_n_cellx
  		do iy=1,gb_n_celly
  		do iz=1,gb_n_cellz
  			if(connect_mat(ix,iy,iz)%marked.eq..false.) then
  				connect_mat(ix,iy,iz)%checked = .true.; cycle
  			elseif(connect_mat(ix,iy,iz)%checked.eq..true.) then
  				if(connect_mat(ix,iy,iz)%crflag .eq. -1) then
  					print *, 'ERROR (strucount)!!! Unchecked, marked pixel shall have crflag!'
  					print *, 'ix,iy,iz,connect_mat = ',ix,iy,iz,connect_mat(ix,iy,iz)
  				endif
  				cycle
  			else
  				nowcrflag = nowcrflag + 1
	  			call conecregion(ix,iy,iz,nowcrflag,coneclev)
  			endif
		enddo
		enddo
		enddo
		numcr = nowcrflag
	end subroutine conmat_cr_scan
	
  !------------------------------------------
  ! Using saddle point to determine whether 
  !  treated as nearby pixels
  !------------------------------------------	
	logical function saddled_nearby(ix1,iy1,iz1,ix2,iy2,iz2)
		integer :: ix1,iy1,iz1,ix2,iy2,iz2
		! TBA!!!!!!!!!!!!!!!!!
		saddled_nearby = .true.
	end function saddled_nearby
	
  !------------------------------------------
  ! Lists of criteria based on fraction
  !------------------------------------------  	
  	subroutine eqvl_fracrilist(frac, quan_list, distance_list, n, rmin, rmax, nbins, cri_list, printinfo)
  		! Dummy
  		real(dl), intent(in) :: frac, quan_list(n), distance_list(n), rmin, rmax
  		integer, intent(in) :: n, nbins
  		real(dl), intent(out) :: cri_list(nbins)
  		logical, intent(in) :: printinfo
  		! Local
		integer :: ibins, i,j, nownum, nowfracnum
		real(dl) :: rmincube, totvol, deltavol, nowrcube1, nowrcube2, rcube, r1, r2
		real(dl), allocatable :: nowquanlist(:)
		
		if(printinfo) write(*,'(A,f7.4,1x,i9,1x,i3)') '   (eqvl_fracrilist) Getting fractional criterias: Frac, n, nbins = ', &
					frac, n, nbins
  		rmincube = rmin**3.0
  		totvol = rmax**3.0 - rmincube
  		deltavol = totvol / dble(nbins)

  		do ibins = 1, nbins
  			nowrcube1 = rmincube + deltavol*(ibins-1)
  			nowrcube2 = nowrcube1 + deltavol
  			nownum = 0
  			do i = 1, n
  				rcube = distance_list(i)**3.0
  				if(rcube < nowrcube2 .and. rcube > nowrcube1) nownum = nownum+1
  			enddo
  			allocate(nowquanlist(nownum))
  			j = 0
  			do i = 1, n
  				rcube = distance_list(i)**3.0
  				if(rcube < nowrcube2 .and. rcube > nowrcube1) then 
  					j = j+1
	  				nowquanlist(j) = quan_list(i)
	  			endif
  			enddo
			if(j.ne.nownum) then
				print *, 'ERROR (eqvl_fracrilist)!!! # from two scans not match: ', j, nownum
				stop
			endif
			nowfracnum = int(nownum * frac)
			call ltquan2(nowquanlist,nownum,nowfracnum,cri_list(ibins),3)
			if(printinfo) then
				write(*,'(22x,i2,A,i9,A,f9.3,A,f9.3,A,f10.5)') ibins,'th bin: # = ', nownum,'; r =',&
					(nowrcube1)**(1.0/3.0),' to',(nowrcube2)**(1.0/3.0),'; criteria = ',cri_list(ibins)
			endif
			call count_num_smaller(nowquanlist,nownum,cri_list(ibins),j)
			if(printinfo) then
				write(*,'(22x,A,i9,i9)') 'Check cri: expected, obtained # smaller than cri = ', nowfracnum, j
				if(abs(nowfracnum - j) > 3) then
					print *, 'ERROR (eqvl_fracrilist)!!! # smaller than cri not right: ', nowfracnum, j
					stop
				endif
			endif
			deallocate(nowquanlist)
		enddo
	end subroutine eqvl_fracrilist


  !------------------------------------------
  ! Calculate saddle mats for usage
  !------------------------------------------ 	
	subroutine calc_saddle_mat(rmin,rmax,nbins,rho_cri_list,printinfo)
		! Dummy
		real(dl), intent(in) :: rmin,rmax,rho_cri_list(nbins)
		integer, intent(in) :: nbins
		logical, intent(in) :: printinfo
		! Local
		integer :: ix,iy,iz,i,j,k,binnum
		real(dl) :: rho,x,y,z,rsq,rmincube,totvol,deltavol
		! Used for test
		logical, parameter :: testsaddlemat = .false.
		
		if(allocated(saddle_rho)) deallocate(saddle_rho)
		if(allocated(check_saddle_rho)) deallocate(check_saddle_rho)
  		allocate(saddle_rho(gb_n_cellx-1,gb_n_celly-1,gb_n_cellz-1), &
  			check_saddle_rho(gb_n_cellx-1,gb_n_celly-1,gb_n_cellz-1))
		saddle_rho = 0.0_dl
  		do ix=1,gb_n_cellx-1; 
  		do iy=1,gb_n_celly-1; 
  		do iz=1,gb_n_cellz-1
  			call nb_fixmd_list_saddle(ix,iy,iz,rho); 
  			saddle_rho(ix,iy,iz) = rho
  		enddo	
  		enddo
  		enddo
	  		
  		! XXXXXXXXXX Some Testings
  		if(testsaddlemat) then
			print *, 'Please Input some ix,iy,iz for test...'
			print *, 'Maximal valus of ix,iy,iz: ', gb_n_cellx-1, gb_n_celly-1, gb_n_cellz-1
			read(*,*) ix,iy,iz
	  		do while(ix>0 .and. iy>0 .and. iz>0)
	  			if(gb_cell_mat(ix,iy,iz)%rhodrho_has_er) print *, '  Warning, this pixel has bd effect.'
				x = gbgridxmin + dble(ix)*gbdeltax
		  		y = gbgridymin + dble(iy)*gbdeltay
		  		z = gbgridzmin + dble(iz)*gbdeltaz
	  			write(*,'(8x,A,4f12.3)'), 'pos/rho at saddle point:', x,y,z,saddle_rho(ix,iy,iz)
				rho = 0.0
	  			do i=ix,ix+1
	  			do j=iy,iy+1
	  			do k=iz,iz+1
	  				call cell_pos(i,j,k,x,y,z)
	  				if(gb_cell_mat(i,j,k)%rhodrho_has_er) print *, '  Warning, pixel has bd effect:', i,j,k
		  			write(*,'(8x,A,4f12.3)'), 'pos/rho at nearby point:', x,y,z,gb_cell_mat(i,j,k)%rhodrhos(1)
		  			rho = rho + gb_cell_mat(i,j,k)%rhodrhos(1)
		  		enddo
		  		enddo
		  		enddo
	 			write(*,'(8x,A,f12.3)'), '<rho> at nearby cells:', rho/8.0
	 			print *, 'Please Input some ix,iy,iz for test...'
				print *, 'Maximal valus of ix,iy,iz: ', gb_n_cellx-1, gb_n_celly-1, gb_n_cellz-1
				read(*,*) ix,iy,iz
			enddo
		endif

		! check whether saddle points satisfying the criteria
		rmincube=rmin**3.0; totvol=rmax**3.0-rmincube; deltavol=totvol/dble(nbins);
		if(printinfo) print *, '  (strucount) Checking saddle points...'
  		do ix=1,gb_n_cellx-1 
  		do iy=1,gb_n_celly-1 
  		do iz=1,gb_n_cellz-1
			x = gbgridxmin + dble(ix)*gbdeltax
	  		y = gbgridymin + dble(iy)*gbdeltay
	  		z = gbgridzmin + dble(iz)*gbdeltaz
  			rsq =x*x+y*y+z*z
  			binnum = max(min(int((rsq**1.5_dl-rmincube)/deltavol)+1,nbins),1)
			if(saddle_rho(ix,iy,iz) > rho_cri_list(binnum)) then
				check_saddle_rho(ix,iy,iz) = .true.
			else
				check_saddle_rho(ix,iy,iz) = .false.
			endif
  		enddo	
  		enddo
  		enddo
  	end subroutine
  	
  	
  !------------------------------------------
  ! estimating rho and gradient rho based on
  !  cubic spline kernel; fixed radius;
  ! Used for saddle points
  !------------------------------------------
  	subroutine nb_fixmd_list_saddle(ix,iy,iz,rho)
		! DUMMY ARGUMENTS
		integer, intent(in) :: ix,iy,iz
!  		integer, intent(out) :: num
  		real(dl), intent(out) :: rho
  		! LOCAL VARIABLES
  		integer, parameter :: max_num = 100000 !maximal # of neary halos 
  		real(dl) :: x,y,z,distance_array(max_num), xyz_mass_array(4,max_num)
  		integer :: di, i,j,k,l, nownum,nowindex, num
  		real(dl) :: r0(3), h, nowr, mass

		x = gbgridxmin + dble(ix)*gbdeltax
  		y = gbgridymin + dble(iy)*gbdeltay
  		z = gbgridzmin + dble(iz)*gbdeltaz

		di = floor(gb_fixmd / gbdeltax)

		if(has_boundary_effect(x,y,z,gb_fixmd,1.0_dl)) then
			rho = 0.0
			return
		endif
		
		num = 0
		r0(1)=x; r0(2)=y; r0(3)=z;
		
		do i = max(1,ix-di), min(ix+di+1,gb_n_cellx)
		do j = max(1,iy-di), min(iy+di+1,gb_n_celly)
		do k = max(1,iz-di), min(iz+di+1,gb_n_cellz)
			do l = 1, gb_cell_mat(i,j,k)%halo_num
				nowindex = gb_cell_mat(i,j,k)%list(l)
				nowr = distance(gb_xyz_list(1:3,nowindex),r0)
				if(nowr < gb_fixmd) then
					if(gb_do_seg_cut.and.nowr>gb_seg_cut_dist) cycle ! Apply seg_cut
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

		rho = 0
		if(num .le. gb_min_smnum) then
			return
		endif	
			
  		h = gb_fixmd / 2.0

		do i = 1, num
			nowr = distance_array(i)
			mass = xyz_mass_array(4,i)
!			print *, i, mass !TESTING			
			rho = rho + mass*w_kernel(nowr, h)
		enddo
!		stop !TESTING
	end subroutine nb_fixmd_list_saddle
end module ap_structure_count

