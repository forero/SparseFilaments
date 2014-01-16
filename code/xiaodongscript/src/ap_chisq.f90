
module ap_chisq

use ap_grad_fields
implicit none

!!! IMPORTANT SETTINGS
! settings for velocity correction
	logical  :: gb_do_vcor = .true.
	real(dl), public :: gb_vcor_fixmd = 40.0_dl
	logical, parameter, public :: gb_vcor_do_seg_cut = .true. ! Enforce segment cut in vel-correction
	real(dl), public :: gb_vcor_seg_cut_dist = 5.0_dl
	integer :: gb_num_vcor = 1

	logical :: gb_chisq_initied = .false.
	logical :: bf_dotsbe
	integer, public :: gb_chisq_method = 3
	integer, public :: gb_num_changenuminx = 15
	real(dl), public :: gb_amp_changenuminx = 0.05_dl
	logical :: output_mu_info = .false.
	logical :: use_mu_corect = .true.
	real(dl), public :: mucorect1=0.0, mucorect2=0.0
	! 1: (<abs(mu)>-0.5)**2.0 / var(abs(mu))
	! 2: histogram, nbins.
	
contains
	
  !------------------------------------------
  ! estimating gradient of fields of 
  !  rho/delta/normed delta at grid points
  !------------------------------------------
	subroutine gd_mldprho_chi2s(RSD, AP, cs, changenuminx, chisqlist, dfchisqlist, multdfchisqlist, numchisq)
		! Dummy
		integer, intent(in) :: RSD, AP, numchisq
		type(chisq_settings) :: cs
		real(dl), intent(in) :: changenuminx
		real(dl), intent(out) :: chisqlist(numchisq), dfchisqlist(numchisq)
		real(dl), intent(out) :: multdfchisqlist(nbinchisq,numchisq)
		! Local
		real(dl) :: rmin,rmax
		integer :: i,j,k, idrop, n1,n2, num_bins, num_bins_test, num
		real(dl), allocatable :: pos_list0(:,:), rho_list0(:), drho_list0(:,:), drho_mu_data(:), &
			pos_list(:,:), rho_list(:), drho_list(:,:),  absdrholist(:), reflist(:)
		integer, allocatable :: markdrop(:)
		logical :: print_time = .false.
		real(dl) :: tmp, time0, time1, time2, time3, time4, time5, time6, time7
		!variables used to estimate mu at different r bins
		real(dl), allocatable :: r_list(:), muavlist(:), muerlist(:), ravlist(:), muvarlist(:), &
			muavlist_test(:), muerlist_test(:), ravlist_test(:), muvarlist_test(:)
		!tmp variables used for outputing info of gradient field
		real(dl) :: xyz1,xyz2,r1,r2,tsbedropratio,nowx,nowy,nowz,nowrsq,rmaxsq,invrsq,rhomean, rhovar, tmpvar
		character(len=char_len) :: str1, tsbefile, file1, file2
	
		if(numchisq.ne.cs%numdrop) then
			print *, 'ERROR (gd_mldprho_chi2s)! size of chisqlist mismathces than ', cs%numdrop
			stop
		endif

		call cpu_time(time0)
		
		if(cs%print_info) then
			write(*,'(A,f9.4,A,i3)') '   (gd_mldprho_chi2s) Estimating chisq: changenuminx = ', real(changenuminx), &
				'; # of drop = ', cs%numdrop
			write(*,'(22x,A,i2,A,<nbinchisq>(i2,1x))') 'Settings of bins (', nbinchisq, ' options): ', nbinchisqlist
		endif

!		TESTING theoretical velocity
		call init_mult_lists(RSD, AP, cs%print_info)
		call do_cell_init(RSD, AP, real(cs%num_in_x)*(1.0_dl+changenuminx), cs%print_info)		
		
		write(str1,'(f15.2)')  real(cs%num_in_x)*(1.0_dl+changenuminx)
		basestr_om_cube_RSD = trim(adjustl(basestr_om))//'_'//trim(adjustl(str1))//'cube_'
		if(RSD.eq.0) then
			basestr_om_cube_RSD = trim(adjustl(basestr_om_cube_RSD))//'NORSD'
		elseif(gb_do_vcor .eq. .false.) then
			basestr_om_cube_RSD = trim(adjustl(basestr_om_cube_RSD))//'WithRSD'
		else
			basestr_om_cube_RSD = trim(adjustl(basestr_om_cube_RSD))//'WithRSDVcor'
		endif
		
		if(cs%print_info) then
			write(*,'(A)') '  (gd_mldprho_chi2s) basestr_om_cube_RSD = '//trim(adjustl(basestr_om_cube_RSD))
		endif
		
		call grid_rho_drho_list(cs%smnum, cs%print_info, pos_list0, rho_list0, drho_list0)
!		print *, '(gd_mldprho_chi2s) Test B...'
		call cpu_time(time1)	
		
		call get_mean_var(rho_list0, rhomean, rhovar)
		n1 = size(rho_list0)
		do i = 1, n1
			rho_list0(i) = rho_list0(i) - 1
		enddo
		if(cs%print_info) then
			write(*,'(A,i8,2(1x,e14.7,1x))') '   (gd_mldprho_chi2s) # / MEAN / VAR of deltas = ', &
				 n1, real(rhomean), real(rhovar)
		endif
		
!		print *, '(gd_mldprho_chi2s) Test C...'
		allocate(absdrholist(n1),markdrop(n1),reflist(n1))
		markdrop = 0
		do i = 1, n1
			absdrholist(i)=sqrt(drho_list0(1,i)**2.0+drho_list0(2,i)**2.0+drho_list0(3,i)**2.0)
		enddo

		! converting rho to delta!!!! 
		!  This has been done in init_mult_lists
!		allocate(r_list(n1))
!		do i = 1, n1
!			r_list(i) = sqrt(pos_list0(1,i)**2.0+pos_list0(2,i)**2.0+pos_list0(3,i)**2.0)
!		enddo
!		call avg_mlist(r_list, rho_list0, n1, cs, 1)
!		call medrho_mlist(r_list, rho_list0, n1, cs%print_info, nbinrhocri, 1)
!		call get_val_dval_list(cs%smnum, cs%print_info, pos_list0, rho_list0, drho_list0, n1)
!		deallocate(r_list)
		
		num_bins = 2
		allocate(muavlist(num_bins), muerlist(num_bins), ravlist(num_bins), muvarlist(num_bins))
		num_bins_test = 2
		allocate(muavlist_test(num_bins_test), muerlist_test(num_bins_test),&
			ravlist_test(num_bins_test), muvarlist_test(num_bins_test)) ! for tests
		if(cs%print_info) write(*,'(A)') '   (gd_mldprho_chi2s) Estimating multi-dropped chisqs...'
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
				
				n2 = 0
				do i = 1, n1
					if(markdrop(i).eq.0) n2 = n2+1
				enddo
				if(cs%print_info.or.gradfieldprintinfo) then
					write(*,'(20x,A,f6.2,1x,f6.2,A,f6.2,1x,f6.2)') &
						'  Drop rho: low/high  =     ', &
						cs%lowdropvalratio(idrop), cs%highdropvalratio(idrop), &
						';  Drop drho: low/high =     ', &
						cs%lowdropdvalratio(idrop), cs%highdropdvalratio(idrop)
					write(*,'(20x,A,i7,5x,f8.5)'), '  Left_# / drop_ratio  in drop:        ', n2, (n1-n2)/dble(n1)
				endif

				call drop_pixels2(pos_list0, rho_list0, drho_list0, markdrop, pos_list, rho_list, drho_list)
				call get_mu_from_gradient_list(pos_list, drho_list, drho_mu_data)

				allocate(r_list(size(rho_list)))
				do i = 1, size(r_list)
					r_list(i) = sqrt(pos_list(1,i)**2.0+pos_list(2,i)**2.0+pos_list(3,i)**2.0)
				enddo
				call find_min_max(r_list, size(r_list), rmin, rmax)
				
				! Output data for test: mu, rho, r, xyz, 
				if(dotsbe) then
					! Only do tsbe for 1th and 5th dropping...
					if((idrop.eq.1.or.idrop.eq.5))  then
						write(str1,'(f5.2)') gb_omegam
						tsbefile = trim(adjustl(tsbestr))//'om'//trim(adjustl(str1))
						write(str1,'(f5.2)') cs%highdropvalratio(idrop)
						tsbefile = trim(adjustl(tsbefile))//'_droprat'//trim(adjustl(str1))
						if(RSD.eq.0) then
							tsbefile = trim(adjustl(tsbefile))//'_NORSD'
						elseif(gb_do_vcor.eq..false.) then
							tsbefile = trim(adjustl(tsbefile))//'_WithRSD'
						else
							tsbefile = trim(adjustl(tsbefile))//'_WithRSD_Vcor'
						endif
						! Save mu into file
						file1 = trim(adjustl(tsbefile))//'_mu.txt'
						if(cs%print_info) write(*,'(21x,A)') 'mu  data saved: '//trim(adjustl(file1))
						call output_1d(file1,drho_mu_data,size(drho_mu_data))
						! Save rho into file
						file1 = trim(adjustl(tsbefile))//'_rho.txt'
						if(cs%print_info) write(*,'(21x,A)') 'rho data saved: '//trim(adjustl(file1))
						call output_1d(file1,rho_list,size(rho_list))
						! Save r into file
						file1 = trim(adjustl(tsbefile))//'_r.txt'
						if(cs%print_info) write(*,'(21x,A)') 'r   data saved: '//trim(adjustl(file1))
						call output_1d(file1,r_list,size(r_list))
						! Save xyz into file
						file1 = trim(adjustl(tsbefile))//'_xyz.txt'
						if(cs%print_info) write(*,'(21x,A)') 'xyz data saved: '//trim(adjustl(file1))
						call output_2d(file1,pos_list)
						! Save drho/dxyz to file
						file1 = trim(adjustl(tsbefile))//'_drho.txt'
						if(cs%print_info) write(*,'(21x,A)') 'drho data saved: '//trim(adjustl(file1))
						call output_2d(file1,drho_list)
					endif
				endif

				! chisq of mu
				chisqlist(idrop) = chisq_of_mu_data2(drho_mu_data)
			
				! chisq from differantial of mu
				do i = 1, size(drho_mu_data)
					drho_mu_data(i) = abs(drho_mu_data(i))
				enddo
			  	call eqvl_binned_quan(drho_mu_data, r_list, rmin, rmax, size(drho_mu_data), num_bins, &	
  					muavlist, muerlist, ravlist, muvarlist)
				if(cs%print_info) & write(*,'(22x,A,2e14.7)') '  <mu> at two vlumes:            ', &
							muavlist(1), muavlist(2)
				if(idrop.eq.1.and.use_mu_corect.and.RSD.eq.1) then
					if(cs%print_info) write(*,'(22x,A,2e14.7)') 'Correcting RSD in mu!!! Correction:', &
						mucorect1, mucorect2
					dfchisqlist(idrop) = abs(muavlist(2) - muavlist(1) - (mucorect2-mucorect1))**2.0 &
						/ ((muerlist(1)**2.0 + (muerlist(2)**2.0)))
					tmpvar = abs(muavlist(2) - muavlist(1))**2.0 &
						/ ((muerlist(1)**2.0 + (muerlist(2)**2.0)))
					if(cs%print_info) then
						write(*,'(22x,A,2e14.7)') '  Chisq before/after correction: ', &
							tmpvar, dfchisqlist(idrop)
					endif
				else
					dfchisqlist(idrop) = abs(muavlist(2) - muavlist(1))**2.0 &
						/ ((muerlist(1)**2.0 + (muerlist(2)**2.0)))
				endif
				call dfchisqs(nbinchisq, nbinchisqlist, drho_mu_data, r_list, n2, multdfchisqlist(1:nbinchisq,idrop))
				if(cs%print_info) then
					write(*,'(25x,A,f13.5)')                  'chisq from ((<mu>-0.5)/sigma)^2     :', &
						real(chisqlist(idrop))
					write(*,'(25x,A,<nbinchisq+1>(f7.3,1x))') 'chisq from comparing bins           :', &
							real(dfchisqlist(idrop)), multdfchisqlist(1:nbinchisq,idrop)
				endif
				! output info of mu
				call eqvl_binned_quan(drho_mu_data, r_list, rmin, rmax, size(drho_mu_data), num_bins_test, &	
  					muavlist_test, muerlist_test, ravlist_test, muvarlist_test)
				if(output_mu_info) then
					write(str1,*) idrop
	  				open(unit=1035,file=trim(adjustl(basestr_om_cube_RSD))//'_drop'&
	  					//trim(adjustl(str1))//'_muinfo.txt')
	  				do i = 1, num_bins_test
	  					write(1035,'(9(e14.7,1x))') gb_omegam, gb_w, real(cs%num_in_x)*(1.0_dl+changenuminx), &
	  						de_zfromintpl(ravlist_test(i)), ravlist_test(i), muavlist_test(i), &
	  						muvarlist_test(i), muerlist_test(i), (muvarlist_test(i)/muerlist_test(i))**2.0+1
	  				enddo
	  				close(1035)
	  			endif
				deallocate(r_list,pos_list,rho_list,drho_list,drho_mu_data)
			endif
		enddo
		deallocate(muavlist,muerlist,ravlist,muvarlist,muavlist_test,muerlist_test,ravlist_test,muvarlist_test)
		deallocate(absdrholist,markdrop,reflist)
		call cpu_time(time2)		


		if(print_time .or. gradfieldprintinfo .or. cs%print_info) then
			write(*,'(A,f12.5)') '   (gd_mldprho_chi2s) Time used in grid_rho_drho_list: ', real(time1-time0)
			write(*,'(A,f12.5)') '   (gd_mldprho_chi2s) Time used in drop: ', real(time2-time1)
		endif
		gradfieldprintinfo = .false.
		if(cs%print_info) write(*,'(A)') '   (gd_mldprho_chi2s) Done.'
	end subroutine gd_mldprho_chi2s


  !------------------------------------------
  ! calculating chisq based
  !------------------------------------------
	  subroutine gf_mldprho_chi2s(omegam, w, h, cs, chisqlist, dfchisqlist, multdfchisqlist, numchisq, calc_comvr)
		! dummy
		real(dl), intent(in) :: omegam, w, h
		integer, intent(in) :: numchisq
		type(chisq_settings), intent(inout) :: cs
		real(dl), intent(out) :: chisqlist(numchisq), dfchisqlist(numchisq), multdfchisqlist(nbinchisq,numchisq)
		logical, intent(in) :: calc_comvr 
		! local
		integer :: i, j, k, RSD, AP, numchanges, num
		real(dl), allocatable :: changenuminxlist(:)
		real(dl) :: tmpchisqlist(numchisq), tmpdfchisqlist(numchisq), tmpmultdfchisqlist(nbinchisq,numchisq)
		! variables used in r-coorection
		real(dl) :: vxth,vyth,vzth,vxcomp,vycomp,vzcomp,vloscomp, ratio,rcorrecth, time1,time2
		real(dl) :: vlosmean, vlosthmean, vlosvar, vlosthvar, vlosratiomean, vlosratiovar, vloscompmean, vloscompvar, &
			vloscompratiomean, vloscompratiovar, rdiffmean1, rdiffvar1, rdiffmean2, rdiffvar2, samesign1, samesign2, samesign3
		integer :: i_vcor, num_has_vcorer, num_no_vcorer
		logical :: touchbdflag, vcorect_er_flag, outputinfo = .false. !!CHECK
		character(len=char_len) :: str1, file1

		if(cs%print_info) then
			write(*,'(A)') '   (gf_mldprho_chi2s) Averaged-chisq estimated from multiple adjustments of grid: '
			print *, '                      # of chisqs   = ',  2*gb_num_changenuminx+1
			if(gb_num_changenuminx.eq.0) then
				print *, '                      maximal ratio = ', 0
			else
				print *, '                      maximal ratio = ', real(gb_amp_changenuminx)
			endif
		endif

		if (numchisq.ne.cs%numdrop) then
			print *, 'ERROR (gf_mldprho_chi2s)! size of chisqlist mismatches with ', cs%numdrop
			stop
		endif

		if (.not. gb_chisq_initied) then
			call cosmo_funs_init()
			call read_in_halo_data()
			call init_halo_info()
			gb_chisq_initied = .true.
		endif

		! take the first cosmology as the input cosmology
		AP = 1
		call init_AP_cosmo(omegam,w,h,AP,cs%print_info)

		if(cs%has_RSD) then
			RSD = 1
		else
			RSD = 0
		endif

		! velocity correction 
		if(RSD .eq. 1 .and. gb_do_vcor) then
			do i_vcor = 1, gb_num_vcor
				vlosmean = 0.0
				vlosthmean = 0.0
				vlosvar = 0.0 
				vlosthvar =0.0 
				vlosratiomean = 0.0 
				vlosratiovar = 0.0
				vloscompmean = 0.0 
				vloscompvar = 0.0
				vloscompratiomean = 0.0
				vloscompratiovar = 0.0
				rdiffmean1 = 0.0
				rdiffvar1 = 0.0 
				rdiffmean2 = 0.0 
				rdiffvar2 = 0.0
				samesign1 = 0.0
				samesign2 = 0.0
				samesign3 = 0.0
				teststr1= '15-100'
				call init_mult_lists(RSD, AP, cs%print_info)
				call do_cell_init(RSD, AP, dble(cs%num_in_x), cs%print_info)
				if(cs%print_info) then
					write(*,'(22x,A,i1,A,i1,A)') 'Correcting RSD... (',i_vcor,'th of ',gb_num_vcor,' iteration)'
					write(*,'(22x,A,f8.2)') ' fixmd = ', gb_vcor_fixmd
					write(*,'(22x,A,L5)')   ' do_seg_cut = ', gb_vcor_do_seg_cut
					write(*,'(22x,A,f8.2)') ' seg_cut_dist = ', gb_vcor_seg_cut_dist
					if(outputinfo) then
						write(*,'(22x,A)') 'Detail Information will be output to file...'
					endif
				endif
				num_has_vcorer = 0
				num_no_vcorer = 0
				if(outputinfo) then
					write(str1,*) i_vcor
					file1=trim(adjustl(teststr1))//'vxth'//trim(adjustl(str1))//'.txt';open(unit=1,file=file1)
					file1=trim(adjustl(teststr1))//'vyth'//trim(adjustl(str1))//'.txt';open(unit=2,file=file1)
					file1=trim(adjustl(teststr1))//'vzth'//trim(adjustl(str1))//'.txt';open(unit=3,file=file1)
					file1=trim(adjustl(teststr1))//'vlosth'//trim(adjustl(str1))//'.txt';open(unit=4,file=file1)
					file1=trim(adjustl(teststr1))//'rcorth'//trim(adjustl(str1))//'.txt';open(unit=5,file=file1)
					file1=trim(adjustl(teststr1))//'vxcomp'//trim(adjustl(str1))//'.txt';open(unit=6,file=file1)
					file1=trim(adjustl(teststr1))//'vycomp'//trim(adjustl(str1))//'.txt';open(unit=7,file=file1)
					file1=trim(adjustl(teststr1))//'vzcomp'//trim(adjustl(str1))//'.txt';open(unit=8,file=file1)
					file1=trim(adjustl(teststr1))//'erflag'//trim(adjustl(str1))//'.txt';open(unit=9,file=file1)
					file1=trim(adjustl(teststr1))//'slice'//trim(adjustl(str1))//'.txt';open(unit=10,file=file1)
				endif

				call cpu_time(time1)
				do i = 1, num_halo
!				do i = num_halo/2, num_halo/2 + 10
					ratio = halo_info(i)%r_AP_RSD_VCOR(AP) / halo_info(i)%r
					call vth_fixmd(halo_info(i)%x*ratio,halo_info(i)%y*ratio,halo_info(i)%z*ratio,&
						num,vxth,vyth,vzth,gb_vcor_fixmd,gb_vcor_do_seg_cut,gb_vcor_seg_cut_dist,touchbdflag,&
						vcorect_er_flag,vxcomp,vycomp,vzcomp)
					call cpu_time(time2)
					if(cs%print_info.and.time2 > time1 + 5) then
						write(*,'(25x,A,f8.2)')  'ratio ', real(i)/real(num_halo)
						time1 = time2
					endif
!					write(*,'(L2,f7.4,1x,3(f9.3,1x),5x,3(f9.3,1x),5x,3(f9.3,1x))') touchbdflag, halo_info(i)%z_real, &
!						halo_info(i)%vx,halo_info(i)%vy,halo_info(i)%vz, &
!						vx,vy,vz, vxcomp,vycomp,vzcomp
					halo_info(i)%vxth = vxth
					halo_info(i)%vyth = vyth
					halo_info(i)%vzth = vzth
					halo_info(i)%vth_er = touchbdflag !.or. vcorect_er_flag ! don't so strict: vcorect_er is 2rd er
					halo_info(i)%vlosth = (halo_info(i)%x*halo_info(i)%vxth + halo_info(i)%y*halo_info(i)%vyth &
						+ halo_info(i)%z*halo_info(i)%vzth) / halo_info(i)%r
					vloscomp = (halo_info(i)%x*vxcomp + halo_info(i)%y*vycomp &
						+ halo_info(i)%z*vzcomp) / halo_info(i)%r
					rcorrecth = (halo_info(i)%vlosth/Hz(dble(halo_info(i)%z_obs))*gb_h)
!					write(*,'(10x,4(f10.3,1x))') halo_info(i)%vlosth, halo_info(i)%vlos, rcorrecth, &
!						halo_info(i)%r_AP_RSD(AP)-halo_info(i)%r_AP(AP)

					if(outputinfo) then
						write(1,*) vxth
						write(2,*) vyth
						write(3,*) vzth
						write(4,*) halo_info(i)%vlosth
						write(5,*) rcorrecth
						write(6,*) vxcomp
						write(7,*) vycomp
						write(8,*) vzcomp
						if(touchbdflag) then
							write(9,*) -1.0
						else
							write(9,*) 1.0
						endif
						if(-30.0<halo_info(i)%x.and.halo_info(i)%x<30.0 .and. halo_info(i)%y>0.0) then
							write(10,'(8(e14.7,1x))') halo_info(i)%y,halo_info(i)%z,halo_info(i)%vy,&
								halo_info(i)%vz,halo_info(i)%vyth,halo_info(i)%vzth,vycomp,vzcomp
						endif
					endif
							
					if(.not. halo_info(i)%vth_er) then
						rdiffmean1 = rdiffmean1 + abs(halo_info(i)%r_AP(AP)-halo_info(i)%r_AP_RSD_VCOR(AP))
						rdiffvar1  = rdiffvar1 + (halo_info(i)%r_AP(AP)-halo_info(i)%r_AP_RSD_VCOR(AP))**2.0
						halo_info(i)%r_AP_RSD_VCOR(AP) = halo_info(i)%r_AP_RSD(AP) - rcorrecth
						num_no_vcorer = num_no_vcorer + 1
						rdiffmean2 = rdiffmean2 + abs(halo_info(i)%r_AP(AP)-halo_info(i)%r_AP_RSD_VCOR(AP))
						rdiffvar2  = rdiffvar2 + (halo_info(i)%r_AP(AP)-halo_info(i)%r_AP_RSD_VCOR(AP))**2.0
						vlosmean = vlosmean + abs(halo_info(i)%vlos)
						vlosthmean = vlosthmean + abs(halo_info(i)%vlosth)
						vloscompmean = vloscompmean + abs(vloscomp)
						vlosratiomean = vlosratiomean + abs((halo_info(i)%vlosth/halo_info(i)%vlos))
						vloscompratiomean = vloscompratiomean + abs(halo_info(i)%vlosth / vloscomp)
						vlosvar = vlosvar + halo_info(i)%vlos**2.0
						vlosthvar = vlosthvar + halo_info(i)%vlosth **2.0
						vloscompvar = vloscompvar + vloscomp**2.0
						vlosratiovar = vlosratiovar + (halo_info(i)%vlosth/halo_info(i)%vlos)**2.0
						vloscompratiovar = vloscompratiovar + (halo_info(i)%vlosth / vloscomp)**2.0
						if(halo_info(i)%vlos*halo_info(i)%vlosth>0) then
							samesign1 = samesign1  + 1.0
						endif
						if(vloscomp*halo_info(i)%vlosth>0) then
							samesign2 = samesign2  + 1.0
						endif
						if(halo_info(i)%vlos*vloscomp>0) then
							samesign3 = samesign3  + 1.0
						endif
					else
						num_has_vcorer = num_has_vcorer + 1
					endif
				enddo
				vlosmean = vlosmean / dble(num_no_vcorer)
				vlosthmean = vlosthmean / dble(num_no_vcorer)
				vloscompmean = vloscompmean / dble(num_no_vcorer)
				vlosratiomean = vlosratiomean / dble(num_no_vcorer)
				vloscompratiomean = vloscompratiomean / dble(num_no_vcorer)
				rdiffmean1 = rdiffmean1 / dble(num_no_vcorer)
				rdiffmean2 = rdiffmean2 / dble(num_no_vcorer)
				vlosvar = vlosvar / dble(num_no_vcorer) - vlosmean**2.0
				vlosthvar = vlosthvar / dble(num_no_vcorer) - vlosthmean**2.0
				vloscompvar = vloscompvar / dble(num_no_vcorer) - vloscompmean**2.0
				vlosratiovar = vlosratiovar / dble(num_no_vcorer) - vlosratiomean**2.0
				vloscompratiovar = vloscompratiovar / dble(num_no_vcorer) - vloscompratiomean**2.0
				rdiffvar1 = rdiffvar1 / dble(num_no_vcorer) - rdiffmean1**2.0
				rdiffvar2 = rdiffvar2 / dble(num_no_vcorer) - rdiffmean2**2.0
				if(cs%print_info) then
					write(*,'(22x,i7,A,f5.2,A)') num_has_vcorer, '(',dble(num_has_vcorer)/dble(num_halo)*100.0, &
						'%) halos have vel-correction error.'
					write(*,'(22x,A,f10.2,A,f15.2,A)') ' |vlos|       = ', vlosmean , &
						'; (rms = ', sqrt(vlosvar), ')'
					write(*,'(22x,A,f10.2,A,f15.2,A)') ' |vlosth|     = ', vlosthmean , &
						'; (rms = ', sqrt(vlosthvar), ')'
					write(*,'(22x,A,f10.2,A,f15.2,A)') ' |vloscomp|   = ', vloscompmean , &
						'; (rms = ', sqrt(vloscompvar), ')'
					write(*,'(22x,A,f10.2,A,f15.2,A)') ' |vlosrat|    = ', vlosratiomean , &
						'; (rms = ', sqrt(vlosratiovar), ')'
					write(*,'(22x,A,f10.3,A,f15.2,A)') ' |vloscomprat|= ', vloscompratiomean , &
						'; (rms = ', sqrt(vloscompratiovar), ')'
					write(*,'(22x,A,f10.3,A,f15.2,A)') ' |delta r| before vcor= ', rdiffmean1 , &
						'; (rms = ', sqrt(rdiffvar1), ')'
					write(*,'(22x,A,f10.3,A,f15.2,A)') ' |delta r| after  vcor= ', rdiffmean2 , &
						'; (rms = ', sqrt(rdiffvar2), ')'
					write(*,'(22x,A,i7,A,f6.2,A)') ' Sign-consistent of vlos & vlosth:     ', &
						int(samesign1), '(',samesign1/ dble(num_no_vcorer)*100.0,'%)'
					write(*,'(22x,A,i7,A,f6.2,A)') ' Sign-consistent of vloscomp & vlosth: ', &
						int(samesign2), '(',samesign2/ dble(num_no_vcorer)*100.0,'%)'
					write(*,'(22x,A,i7,A,f6.2,A)') ' Sign-consistent of vlos & vloscomp:    ', &
						int(samesign3), '(',samesign3/ dble(num_no_vcorer)*100.0,'%)'

				endif
				if(outputinfo) then
					close(1);close(2);close(3);close(4);close(5);close(6);close(7);close(8);close(9);close(10);
				endif
			enddo
			
			! info of halos
			if(outputinfo) then
				open(unit=1,file=trim(adjustl(teststr1))//'x.txt')
				open(unit=2,file=trim(adjustl(teststr1))//'y.txt')
				open(unit=3,file=trim(adjustl(teststr1))//'z.txt')
				open(unit=4,file=trim(adjustl(teststr1))//'raptor.txt')
				open(unit=5,file=trim(adjustl(teststr1))//'vx.txt')
				open(unit=6,file=trim(adjustl(teststr1))//'vy.txt')
				open(unit=7,file=trim(adjustl(teststr1))//'vz.txt')
				open(unit=8,file=trim(adjustl(teststr1))//'vlos.txt')
				open(unit=9,file=trim(adjustl(teststr1))//'vcor.txt')
				do i = 1, num_halo
!					if(.not.gb_vcorect_er_list(i)) then
						write(1,*) halo_info(i)%x
						write(2,*) halo_info(i)%y
						write(3,*) halo_info(i)%z
						write(4,*) halo_info(i)%r_AP(AP) / halo_info(i)%r
						write(5,*) halo_info(i)%vx
						write(6,*) halo_info(i)%vy
						write(7,*) halo_info(i)%vz
						write(8,*) halo_info(i)%vlos
						write(9,*) halo_info(i)%r_AP_RSD(AP)-halo_info(i)%r_AP(AP)
!					endif
				enddo
				close(1);close(2);close(3);close(4);close(5);close(6);close(7);close(8);close(9);
			endif
			
		endif

		!get the chisqs
		numchanges = 2*gb_num_changenuminx+1
		allocate(changenuminxlist(numchanges))
		if(cs%print_info) write(*,'(A)') '   (gf_mldprho_chi2s) Estimating multiple chisqs from different adjustments...'
		changenuminxlist(1) = 0.0_dl
		do i = 1, gb_num_changenuminx
			changenuminxlist(i+1) = gb_amp_changenuminx*i/(gb_num_changenuminx+0.0)
			changenuminxlist(2*gb_num_changenuminx-i+2) = -1.0_dl *  gb_amp_changenuminx*i/(gb_num_changenuminx+0.0)
		enddo
		chisqlist = 0.0_dl
		dfchisqlist = 0.0_dl
		multdfchisqlist = 0.0_dl
		do i = 1, numchanges
			if(i.eq.1) then
				bf_dotsbe = dotsbe
			else
				dotsbe = .false.
			endif
			if(cs%print_info) write(*,'(A,f7.3,A)') '   (gf_mldprho_chi2s) chisqs with adjustment ratio ', &
				changenuminxlist(i),':'
			tmpchisqlist = 0.0
			tmpdfchisqlist = 0.0
			tmpmultdfchisqlist = 0.0
			call gd_mldprho_chi2s(RSD, AP, cs, changenuminxlist(i),tmpchisqlist,tmpdfchisqlist,tmpmultdfchisqlist,cs%numdrop)
			print *, 'Calling gd_mldprho_chi2s done.'
			do j = 1, cs%numdrop
				chisqlist(j) = chisqlist(j) + tmpchisqlist(j)/dble(numchanges)
				dfchisqlist(j) = dfchisqlist(j) + tmpdfchisqlist(j)/dble(numchanges)
				do k = 1, nbinchisq
					multdfchisqlist(k,j) = multdfchisqlist(k,j) + tmpmultdfchisqlist(k,j)/dble(numchanges)
				enddo
			enddo
		enddo
		deallocate(changenuminxlist)
		dotsbe = bf_dotsbe
		if(cs%print_info) then
			write(*,'(A)') '   (gf_mldprho_chi2s) Done.'
		endif
	end subroutine gf_mldprho_chi2s
end module ap_chisq
