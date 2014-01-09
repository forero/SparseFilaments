!####################################
! Useful halo data files:
!
!	halo_data_file = '../../data/input/MDR1_shelldata_sep15.dat'
!	halo_data_file = '../../data/input/MDR1_shelldata_sep10.dat'
!	halo_data_file = '../../data/input/MDR1_shelldata_sep8.dat'
!	halo_data_file = '../../data/input/HR3LC0_1eighth.dat'
!	halo_data_file = '../../data/input/HR3LCQua'//trim(adjustl(LClabelChar))//'_57w_600to1787.dat'; local_pecu_v = .true.
!	halo_data_file = '../../data/input/HR3LC0_NMC_600to1787.dat'
!	halo_data_file = '../../data/input/HR3LCQua'//trim(adjustl(LClabelChar))//'_57w_1200to1787.dat'

!	halo_data_file = '../../data/input/NewHR3LC00_xyzge0.dat'
!	halo_data_file = '../../data/input/NewHR3LC00.dat'			
!	halo_data_file = '../../data/input/NewHR3LC00_xyzge0_rge500.dat'

!	halo_data_file = '../../data/input/NewCBPHR3LCx.dat'
!	halo_data_file = '../../data/input/NewHR3LC00_xyzge0_rge500.dat'
!	halo_data_file = '../../data/input/NewCBPHR3LC00_xyge0_rge500.dat'

!####################################
! scan chisqs
!####################################

program ap_main

use mpi
use ap_chisq
use ap_structure_count

	implicit none

	type(chisq_settings) :: cs
	logical :: auto_num_in_x
	character(len=char_len) :: halodataname, inputfile, outputname, &
		rhofile, rhoRSDfile, dltfile, rhoRSDvcorfile, ndltfile, nrhoRSDvcorfile, str, str2
	integer :: i, j, n, num_omw, tmpnum
	real(dl), allocatable :: om_w_list(:,:)
	! chis
	real(dl), allocatable :: rhochi(:), rhoRSDchi(:), rhodfchi(:), rhoRSDdfchi(:), &
		dltchi(:), dltRSDchi(:), dltdfchi(:), dltRSDdfchi(:), &
		ndltchi(:), ndltRSDchi(:), ndltdfchi(:), ndltRSDdfchi(:), multdfchi(:,:)
!	real(dl) :: rhochi(gb_numwei), rhoRSDchi(gb_numwei)
	real(dl) :: ommin, ommax, wmin, wmax, dom, dw, om, w, time1, time2, timebegin, timeend
	! tmp variables used for testing
	real(dl) :: drmean, drvar, ratio
	real(dl), allocatable :: tmp(:), drho_mu_data(:), ddlt_mu_data(:), dnormed_dlt_mu_data(:)
	real(dl) :: chisq1, chisq2, chisq3
	integer :: i_LC, datalab
	character(len=char_len) :: str1, filename
	! mpi variables
	character(len=char_len) :: specnamestr, mpibaseoutputname
	logical :: alive, use_mpi = .true.
	integer :: ierr, nproc, myid, i1,i2
	real(dl), allocatable :: sendrho(:),sendrhoRSD(:),sendrhoRSDvcor(:), recivrho(:),recivrhoRSD(:),recivrhoRSDvcor(:)

	! Initialize MPI
	call mpi_init(ierr)
	call mpi_comm_size(mpi_comm_world,nproc,ierr)
	call mpi_comm_rank(mpi_comm_world,myid,ierr)
!	call system('sleep 12000')
     do i_LC = 0, 26
	! Settings of data file
	call mpi_barrier(mpi_comm_world,ierr)
	write(LClabelChar,*) i_LC;	
	if(i_LC.le.9) then
		halodataname = 'NewCBPHR3LC0'//trim(adjustl(LClabelChar))//'_zge0_rge500'
	else
		halodataname = 'NewCBPHR3LC'//trim(adjustl(LClabelChar))//'_zge0_rge500'
	endif
	halo_data_file = '../../data/input/'//trim(adjustl(halodataname))//'.dat'
	local_pecu_v = .true.

	! Settings of cut
	gb_use_fixmd = .true.
	gb_fixmd = 25.0_dl 		!!SETTING
	gb_min_smnum = 0 ! minimal (gb_min_smnum + 1) halos to get gradient; .le. gb_min_smnum will be skipped
	gb_keepzerorho = .false.

	! using cuts in smooth kernel (ignore nearby halos)
	gb_do_seg_cut = .true.
	gb_seg_cut_dist = 8.0_dl 	!!SETTING

	! Settings of vel correct
	gb_do_vcor = .true.
	gb_num_vcor = 2
	gb_vcor_fixmd = 30.0_dl 	!!SETTING
	gb_vcor_seg_cut_dist = 10.0_dl  !!SETTING
	
	! Settings of mu correct
	use_mu_corect = .false.
	mucorect1 = 0.015551; mucorect2 = 0.016204 ! 27 mocks, fixmd 0-25, 2 bins result, diff: 0.00065
!	mucorect1 = 0.015562_dl; mucorect2 = 0.016190_dl ! 27 mocks, fixmd 0-25: diff 0.00063
!	mucorect1 = 0.0170273; mucorect2 = 0.0172741 ! result from mock 0
!	mucorect1 = 0.0157746877778_dl; mucorect2 = 0.0164351766667_dl ! 18 mocks, fixmd 10-30
!	mucorect1 = 0.017098_dl; mucorect2 = 0.017746_dl ! 27 mocks, fixmd 0-25, 2 vcor 10-30: diff 0.00065 (vcor makes them larger?!)
	! Settings of method
	if(use_mpi) then
		cs%smnum = 30			!!CHECK
		cs%num_in_x = 240		!!CHECK
		use_num_density = .true.	!!CHECK
		specnamestr = '_use_num_dens_2bins'	!!CHECK
		gb_num_changenuminx = 0 !!CHECK
		gb_amp_changenuminx = 0.05_dl
		num_omw = 1; ommin = 0.26_dl; ommax = 0.26_dl !!CHECK
!		num_omw = 60; ommin = -0.1_dl; ommax = 0.785_dl
!		num_omw = 32; ommin = 0.02_dl; ommax = 0.64_dl
		output_mu_info = .true. !!CHECK!!!MUST BE FALSE IN CLUSTER
		allocate(om_w_list(2,num_omw))
		do i = 1, num_omw
			if(i.eq.1) then !!Specifically treat first omegam. So num_omw will not give Nan Omegam.
				om_w_list(1,i) = ommin
			else
				om_w_list(1,i) = ommin + dble(ommax-ommin)/dble(num_omw-1)*dble(i-1)
			endif
			om_w_list(2,i) = -1.0_dl
		enddo
		cs%numdrop = 6 !!CHECK
		allocate(cs%dropval(cs%numdrop),  cs%lowdropvalratio(cs%numdrop),  cs%highdropvalratio(cs%numdrop), &
			cs%dropdval(cs%numdrop), cs%lowdropdvalratio(cs%numdrop), cs%highdropdvalratio(cs%numdrop))
		! Our conventional setting of 11 dropping ratios : 0,10,20,30,...,80,90,95.
		do i = 1, cs%numdrop
			cs%dropval(i) = .true.;	cs%lowdropvalratio(i) = 0.0
			cs%dropdval(i) = .true.; cs%lowdropdvalratio(i) = 0.0
			cs%highdropvalratio(i) = min((i-1)*0.10_dl,0.95_dl) !!CHECK
			cs%highdropdvalratio(i) = min((i-1)*0.10_dl,0.95_dl)!!CHECK
		enddo
		!!CHECK !!MUST BE CHECKED ON CLUSTER
		mpibaseoutputname     = 'MPI_chisqlike_muinfo_2014Jan07/MPI_chisqlike_iLC'//trim(adjustl(LClabelChar))
		str = 'mkdir -p '//trim(adjustl(mpibaseoutputname)); call system(str)
		mpibaseoutputname     = trim(adjustl(mpibaseoutputname))//'/'//trim(adjustl(halodataname))//'_om_smnum'
                write(str,*) cs%smnum
                mpibaseoutputname      = trim(adjustl(mpibaseoutputname))//trim(adjustl(str))//'_'
                write(str,*) cs%num_in_x
                mpibaseoutputname      = trim(adjustl(mpibaseoutputname))//trim(adjustl(str))//'cube' ! redshift cut 0.36
                !str = 'mkdir -p '//trim(adjustl(mpibaseoutputname)); call system(str);
		if(gb_use_fixmd) then
			write(str,'(f4.1)') gb_fixmd
                	mpibaseoutputname      = trim(adjustl(mpibaseoutputname))//'_fixmd'//trim(adjustl(str)) ! redshift cut 0.36
		endif
		if(gb_do_seg_cut) then
			write(str,'(f4.1)') gb_seg_cut_dist
                	mpibaseoutputname      = trim(adjustl(mpibaseoutputname))//'_segmincut'//trim(adjustl(str)) ! redshift cut 0.36
		endif
		if(gb_do_vcor) then
			gb_do_vcor = .true.
			write(str,*) gb_num_vcor
                	mpibaseoutputname      = trim(adjustl(mpibaseoutputname))//'_'//trim(adjustl(str))//'vcor' ! redshift cut 0.36
			write(str,'(f4.1)') gb_vcor_seg_cut_dist
			mpibaseoutputname      = trim(adjustl(mpibaseoutputname))//trim(adjustl(str))//'to'
			write(str,'(f4.1)') gb_vcor_fixmd 
			mpibaseoutputname      = trim(adjustl(mpibaseoutputname))//trim(adjustl(str))
		endif
		mpibaseoutputname = trim(adjustl(mpibaseoutputname))//trim(adjustl(specnamestr))
                basestr 	= mpibaseoutputname
                call system('mkdir -p '//trim(adjustl(basestr)))
                write(str,*) myid
                outputname      = trim(adjustl(mpibaseoutputname))//'_MPIrank'//trim(adjustl(str))
		call mpi_barrier(mpi_comm_world,ierr)
	else
		i = iargc()	
		if(i .ne. 2) then
			print*, "ERROR! Number of arg must be 2."
			print*, "USAGE: "
			print*, "  ./post_proc inputfile outputfile"
			stop
		endif
		call getarg(1,inputfile)
		call getarg(2,outputname)
		basestr = outputname
		print *, 'basestr = ', trim(adjustl(basestr))
		open(unit=1,file=inputfile)
		read(1,*) cs%smnum, cs%num_in_x, auto_num_in_x, cs%check_boundary, cs%print_info, cs%cb_adjust_ratio, &
			cs%remov_dist_ratio, use_num_density, cs%nbins_rhoav, cs%use_intpl_rho, cs%numdrop, &
			gb_num_changenuminx, gb_amp_changenuminx, num_omw, gb_use_fixmd, gb_fixmd, gb_do_seg_cut, gb_seg_cut_dist
		allocate(cs%dropval(cs%numdrop),  cs%lowdropvalratio(cs%numdrop),  cs%highdropvalratio(cs%numdrop), &
			cs%dropdval(cs%numdrop), cs%lowdropdvalratio(cs%numdrop), cs%highdropdvalratio(cs%numdrop))
		do i = 1, cs%numdrop
			read(1,*) cs%dropval(i), cs%lowdropvalratio(i),  cs%highdropvalratio(i), &
				cs%dropdval(i), cs%lowdropdvalratio(i), cs%highdropdvalratio(i)
			cs%lowdropvalratio(i) = cs%lowdropvalratio(i) / 100.0_dl
			cs%highdropvalratio(i) = cs%highdropvalratio(i) / 100.0_dl
			cs%lowdropdvalratio(i) = cs%lowdropdvalratio(i) / 100.0_dl
			cs%highdropdvalratio(i) = cs%highdropdvalratio(i) / 100.0_dl
		enddo
		allocate(om_w_list(2,num_omw))
		do i = 1, num_omw
			read(1,*) om_w_list(1:2,i)
		enddo
		close(1)
	endif

	write(*,*) 'Settings of the programm:'
	write(*,*) '    smnum             = ', cs%smnum
	write(*,*) '    cube (x-direc)    = ', cs%num_in_x
	write(*,'(A,i4)')    '     # of cellchange   : ', gb_num_changenuminx
	write(*,'(A,e14.7)') '     amp of cellchange : ', gb_amp_changenuminx	
	write(*,*) '    num of drop       = ', cs%numdrop
	write(*,*) '    num of om/w       = ', num_omw
	print *
	write(*,*) '    check_boundary    = ', cs%check_boundary
	write(*,*) '    cb_adjust_ratio   = ', real(cs%cb_adjust_ratio)
	write(*,*) '    remov_dist_ratio  = ', real(cs%remov_dist_ratio)
	write(*,*) '    use_num_density   = ', use_num_density
	write(*,*) '    nbins_rhoav       = ', cs%nbins_rhoav
	write(*,*) '    use_intpl_rho     = ', cs%use_intpl_rho
	
	print *
	write(*,'(A,i2,A)') ' Settings of dropping (',cs%numdrop,'):'
	do i = 1, cs%numdrop
		write(*,'(A,L2,f6.3,f6.3,A,L2,f6.3,f6.3)') '    Density: ', cs%dropval(i), real(cs%lowdropvalratio(i)),  &
			real(cs%highdropvalratio(i)), ';   Density Gradient: ', &
			cs%dropdval(i), real(cs%lowdropdvalratio(i)), real(cs%highdropdvalratio(i))
	enddo
	
	write(*,'(10x,A,<num_omw>(f6.3,1x))') 'Lists of omegam: ', om_w_list(1,1:num_omw)
	write(*,'(10x,A,<num_omw>(f6.3,1x))') 'Lists of w     : ', om_w_list(2,1:num_omw)

	if(gb_use_fixmd) then
		print *
		print *, 'Applying fixed max_dist = ', gb_fixmd
	endif
	
	if(gb_do_seg_cut) then
		print *
		print *, 'Applying seg_cut_dist = ', gb_seg_cut_dist
	endif
	
	if(gb_do_vcor) then
		print *
		write(*,'(i2,A,f4.1,A,f4.1)') gb_num_vcor, ' r-corrections: vel ested in ', gb_vcor_seg_cut_dist, ' to ', gb_vcor_fixmd
	endif
	
	dotsbe = .false.
	if(dotsbe) then
		tsbestr = trim(adjustl(outputname))//'_infos/'
		str1 = 'mkdir -p '//trim(adjustl(tsbestr)); call system(str1)
		write(*,'(A)'), '  Output mu/rho/xyz data: dir='//trim(adjustl(tsbestr))
	endif
	
	print *
	call cosmo_funs_init()
	
	!!TESTING
	if(.false.) then
		gb_omegam = 0.25_dl
		gb_w = -0.5_dl
		gb_h = 0.73_dl
		call de_calc_comovr()
		do i = 1, 300, 15
			write(*,'(4e14.7,2x)') de_zi(i), de_gfz_data(i), omegamz(de_zi(i))**(6.0d0/11.0d0), omegamz(de_zi(i))**(0.618d0)
		enddo
		open(unit=1,file='Test/om0p25_wneg0p5')
		do i = 1, 300, 15
			write(1,'(3e14.7,2x)') de_zi(i), de_gfz_data(i)!, omegamz(de_zi(i))**(6.0d0/11.0d0)
		enddo
		close(1)
		stop
	endif
	
	call read_in_halo_data()
	call init_halo_info()
	gb_chisq_initied = .true.

	allocate(rhochi(cs%numdrop),rhoRSDchi(cs%numdrop),rhodfchi(cs%numdrop),rhoRSDdfchi(cs%numdrop), &
		dltchi(cs%numdrop),dltRSDchi(cs%numdrop),dltdfchi(cs%numdrop),dltRSDdfchi(cs%numdrop), &
		ndltchi(cs%numdrop),ndltRSDchi(cs%numdrop),ndltdfchi(cs%numdrop),ndltRSDdfchi(cs%numdrop), &
		multdfchi(nbinchisq,cs%numdrop),tmp(nbinchisq*cs%numdrop),&
		sendrho((2+nbinchisq)*cs%numdrop),sendrhoRSD((2+nbinchisq)*cs%numdrop),&
		sendrhoRSDvcor((2+nbinchisq)*cs%numdrop),&
		recivrho((2+nbinchisq)*cs%numdrop*nproc),recivrhoRSD((2+nbinchisq)*cs%numdrop*nproc),&
		recivrhoRSDvcor((2+nbinchisq)*cs%numdrop*nproc) )

	str = outputname
	rhofile = trim(adjustl(str))//'_rho.txt'
	rhoRSDfile = trim(adjustl(str))//'_rhoRSD.txt'
	rhoRSDvcorfile = trim(adjustl(str))//'_rhoRSD_vcor.txt'
	
	print *
	print *, 'Calculating chisqs...'	
	print *, '  Results saved in the following files:'
	write(*,'(A)'), '     '//trim(adjustl(rhofile))
	write(*,'(A)'), '     '//trim(adjustl(rhoRSDfile))
	print *

	open(unit=1,file=rhofile); open(unit=2,file=rhoRSDfile); open(unit=3,file=rhoRSDvcorfile)
	call cpu_time(timebegin)
	time1 = timebegin

	cs%print_info = .true.
	print *, 'Scanning chisqs...'
	do i =  myid+1, num_omw, nproc
		om = om_w_list(1,i)
		w  = om_w_list(2,i)

		write(*,'(A)') 'basestr = '//trim(adjustl(basestr))
		write(str1,'(f6.3)') om
		basestr_om = trim(adjustl(basestr))//'/om'//trim(adjustl(str1))
		write(*,'(A)') ' basestr_om = '//trim(adjustl(basestr_om))

		cs%has_RSD = .false.
		call gf_mldprho_chi2s(om, w, h_dft, cs, rhochi, rhodfchi, multdfchi, cs%numdrop, calc_comvr = .true.)
		write(*,'(A,i3,A,<2>(f7.2,1x))') '   Step ',i,'. Chisqs without RSD:', om, w
		write(*,'(34x,<cs%numdrop>(f7.2,1x))')	rhochi(1:cs%numdrop)
		write(*,'(A,12x,<cs%numdrop>(f7.2,1x))') '    2 bins (old method):', rhodfchi(1:cs%numdrop)
		do j = 1, nbinchisq	
			write(*,'(4x,i2,A,24x,<cs%numdrop>(f7.2,1x))') nbinchisqlist(j),' bins:', multdfchi(j,1:cs%numdrop)
			tmp((j-1)*cs%numdrop+1:j*cs%numdrop) = multdfchi(j,1:cs%numdrop)
		enddo
		sendrho(1:cs%numdrop) = rhochi(1:cs%numdrop)
		sendrho(cs%numdrop+1:2*cs%numdrop) = rhodfchi(1:cs%numdrop)
		sendrho(2*cs%numdrop+1:(2+nbinchisq)*cs%numdrop) = tmp(1:nbinchisq*cs%numdrop)
		write(1,'(<2+(2+nbinchisq)*cs%numdrop>(e14.7,1x))') om, w,rhochi(1:cs%numdrop),rhodfchi(1:cs%numdrop), tmp

		cs%has_RSD = .true.
		gb_do_vcor = .false.
		call gf_mldprho_chi2s(om, w, h_dft, cs, rhoRSDchi, rhoRSDdfchi, multdfchi, cs%numdrop, calc_comvr = .true.)
		write(*,'(A,i3,A,<2>(f7.2,1x))') '   Step ',i,'. Chisqs  with   RSD (No Vcor):', om, w
		write(*,'(34x,<cs%numdrop>(f7.2,1x))')	rhoRSDchi(1:cs%numdrop)
		write(*,'(A,12x,<cs%numdrop>(f7.2,1x))') '    2 bins (old method):', rhoRSDdfchi(1:cs%numdrop)
		do j = 1, nbinchisq	
			write(*,'(4x,i2,A,24x,<cs%numdrop>(f7.2,1x))') nbinchisqlist(j),' bins:', multdfchi(j,1:cs%numdrop)
			tmp((j-1)*cs%numdrop+1:j*cs%numdrop) = multdfchi(j,1:cs%numdrop)
		enddo
		sendrhoRSD(1:cs%numdrop) = rhoRSDchi(1:cs%numdrop)
		sendrhoRSD(cs%numdrop+1:2*cs%numdrop) = rhoRSDdfchi(1:cs%numdrop)
		sendrhoRSD(2*cs%numdrop+1:(2+nbinchisq)*cs%numdrop) = tmp(1:nbinchisq*cs%numdrop)
		write(2,'(<2+(2+nbinchisq)*cs%numdrop>(e14.7,1x))') om, w, rhoRSDchi(1:cs%numdrop), rhoRSDdfchi(1:cs%numdrop), tmp
		
		cs%has_RSD = .true.
		gb_do_vcor = .true.
		call gf_mldprho_chi2s(om, w, h_dft, cs, rhoRSDchi, rhoRSDdfchi, multdfchi, cs%numdrop, calc_comvr = .true.)
		write(*,'(A,i3,A,<2>(f7.2,1x))') '   Step ',i,'. Chisqs  with   RSD (Vcor):', om, w
		write(*,'(34x,<cs%numdrop>(f7.2,1x))')	rhoRSDchi(1:cs%numdrop)
		write(*,'(A,12x,<cs%numdrop>(f7.2,1x))') '    2 bins (old method):', rhoRSDdfchi(1:cs%numdrop)
		do j = 1, nbinchisq	
			write(*,'(4x,i2,A,24x,<cs%numdrop>(f7.2,1x))') nbinchisqlist(j),' bins:', multdfchi(j,1:cs%numdrop)
			tmp((j-1)*cs%numdrop+1:j*cs%numdrop) = multdfchi(j,1:cs%numdrop)
		enddo
		sendrhoRSDvcor(1:cs%numdrop) = rhoRSDchi(1:cs%numdrop)
		sendrhoRSDvcor(cs%numdrop+1:2*cs%numdrop) = rhoRSDdfchi(1:cs%numdrop)
		sendrhoRSDvcor(2*cs%numdrop+1:(2+nbinchisq)*cs%numdrop) = tmp(1:nbinchisq*cs%numdrop)		
		write(3,'(<2+(2+nbinchisq)*cs%numdrop>(e14.7,1x))') om, w, rhoRSDchi(1:cs%numdrop), rhoRSDdfchi(1:cs%numdrop), tmp

		call cpu_time(time2)
		if(gbtp) then
			write(*,'(A,f10.4)') '  Time used in this step: ', time2-time1
			time1 = time2
			print *
		endif
		cs%print_info = .false.
	enddo
	
	call cpu_time(timeend)
	close(1); close(2); close(3);

	if(use_mpi) then
		print *, 'Process ',myid,'done.'
		print *, '  Total time:     ',  timeend - timebegin
		print *, '  Total time / n: ', (timeend - timebegin) / (num_omw+0.0)
		call mpi_barrier(mpi_comm_world,ierr)
		call mpi_gather(sendrho,(2+nbinchisq)*cs%numdrop,MPI_DOUBLE_PRECISION,&
				recivrho,(2+nbinchisq)*cs%numdrop,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)
		call mpi_barrier(mpi_comm_world,ierr)
		call mpi_gather(sendrhoRSDvcor,(2+nbinchisq)*cs%numdrop,MPI_DOUBLE_PRECISION,&
				recivrhoRSDvcor,(2+nbinchisq)*cs%numdrop,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)
		call mpi_barrier(mpi_comm_world,ierr)
		call mpi_gather(sendrhoRSD,(2+nbinchisq)*cs%numdrop,MPI_DOUBLE_PRECISION,&
				recivrhoRSD,(2+nbinchisq)*cs%numdrop,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)
		call mpi_barrier(mpi_comm_world,ierr)
		if(myid.eq.0) then
			str = mpibaseoutputname
			rhofile = trim(adjustl(str))//'_rho.txt'
			rhoRSDfile = trim(adjustl(str))//'_rhoRSD.txt'
			rhoRSDvcorfile = trim(adjustl(str))//'_rhoRSD_vcor.txt'
			open(unit=1,file=rhofile); open(unit=2,file=rhoRSDfile); open(unit=3,file=rhoRSDvcorfile);
			do i = 1, num_omw
				om = om_w_list(1,i)
				w  = om_w_list(2,i)
				i1 = (2+nbinchisq)*cs%numdrop*(i-1)+1
				i2 = (2+nbinchisq)*cs%numdrop*i
				write(1,'(<2+(2+nbinchisq)*cs%numdrop>(e14.7,1x))') om, w,recivrho(i1:i2)
				write(2,'(<2+(2+nbinchisq)*cs%numdrop>(e14.7,1x))') om, w,recivrhoRSD(i1:i2)
				write(3,'(<2+(2+nbinchisq)*cs%numdrop>(e14.7,1x))') om, w,recivrhoRSDvcor(i1:i2)
			enddo
			close(1);close(2);close(3);
		endif
		call mpi_barrier(mpi_comm_world,ierr)
	else
		print *, 'Total time:     ',  timeend - timebegin
		print *, 'Total time / n: ', (timeend - timebegin) / (num_omw+0.0)
	endif

	deallocate(om_w_list, rhochi,rhoRSDchi,rhodfchi,rhoRSDdfchi,dltchi,&
		dltRSDchi,dltdfchi,dltRSDdfchi,ndltchi,ndltRSDchi,ndltdfchi,ndltRSDdfchi,multdfchi,tmp)
	deallocate(cs%dropval,  cs%lowdropvalratio, cs%highdropvalratio, &
			cs%dropdval, cs%lowdropdvalratio, cs%highdropdvalratio)
	deallocate(sendrho,sendrhoRSD,sendrhoRSDvcor,recivrho,recivrhoRSD,recivrhoRSDvcor)
     enddo
	call mpi_barrier(mpi_comm_world,ierr)
	call mpi_finalize(ierr)
end program ap_main
