
program ap_main_test

!use mpi
use ap_tools
use ap_chisq
use ap_grids

	implicit none
	
	real(dl) :: refDx, omegam, rl_num_in_x
	integer, parameter :: numbins = 2, numcube = 30, ndrop = 10
	character(len=char_len) :: str, outputname, file1, file2
	type(mustatype) :: mustat, NORSDmustat(numcube), RSDmustat(numcube)
	integer :: i, numom
	real(dl) :: ommin, ommax, dom
	real(dl) :: droparray(10) = (/0.0_dl, 0.1_dl, 0.2_dl, 0.3_dl, 0.4_dl, 0.5_dl, 0.6_dl, 0.7_dl, 0.8_dl, 0.9_dl/)
	real(dl) :: cubemin, cubemax, drmin, drmax, NORSDdiffmu(numcube), RSDdiffmu(numcube), &
		chisq1(ndrop), chisq2(ndrop), NORSDchisq(numcube), RSDchisq(numcube)
	logical :: calc_cos 
	
!	call getarg(1,str)
!	read(str,*) omegam
!	call random_seed()
	halo_data_file = '../../../data/input/HR3LC0_57w_600to1787.dat'
!	halo_data_file = '../../../data/input/HR3halo100w_pos1000.dat' ! Successful if use 100-150 distance cut, RSD 
!	halo_data_file = '../../../data/input/MD_fullhalos_z0_100w.dat'
	call cosmo_funs_init()
	call read_in_halo_data()
	call init_halo_info()	

	ommin = 0.0 
	ommax = 0.5
	numom = 10
	dom = (ommax - ommin) / (numom - 1.0_dl)

	cubemin=40.0_dl
	cubemax=60.0_dl
	calc_cos = .true.
	outputname = '../CICchisq/CICchisq_HR3LC0B_cube'
	write(str,*) int(cubemin)
	outputname = trim(adjustl(outputname))//trim(adjustl(str))
	write(str,*) int(cubemax)
	outputname = trim(adjustl(outputname))//'to'//trim(adjustl(str))//'_num'
	write(str,*) numcube
	outputname = trim(adjustl(outputname))//trim(adjustl(str))//'_'
	file1=trim(adjustl(outputname))//'NORSD.txt'
	file2=trim(adjustl(outputname))//'RSD.txt'

	open(unit=1,file=file1)
	open(unit=2,file=file2)
	do i = 1, numom
		omegam = ommin + (i-1)*dom
		call CICgfchisq(omegam, w_dft, h_dft, cubemin, cubemax, numcube, calc_cos, droparray, ndrop, chisq1, chisq2)
		write(1,'(<2*ndrop+1>(e14.7,1x))') omegam, chisq1(1:ndrop)
		write(2,'(<2*ndrop+1>(e14.7,1x))') omegam, chisq2(1:ndrop)
		write(*,'(<2*ndrop+1>(e14.7,1x))') omegam, chisq1(1:ndrop)
		write(*,'(<2*ndrop+1>(e14.7,1x))') omegam, chisq2(1:ndrop)
	enddo
	close(1)
	close(2)
	stop

	
	! Get the refDx from the right cosmology
	call init_pixels(RSD=0, AP=0, rl_num_in_x=50.5_dl, &
		print_info=.false., use_num_density=.true., do_xyzmassinit=.true.)
	call do_grid_avg_mlist()
	call init_pixels(RSD=0, AP=0, rl_num_in_x=50.5_dl, &
		print_info=.false., use_num_density=.false., do_xyzmassinit=.true.)
	refDx=gbxmax-gbxmin

	! Testing mu histograms (sorry we cannot find anything special ...)
	drmin = 42.0_dl
	drmax = 44.0_dl
	if(.true.) then
		rl_num_in_x = 60.0_dl
		call calc_mustat(refDx, 0.26_dl, drmin, drmax, rl_num_in_x)
		call calc_mustat(refDx, 0.05_dl, drmin, drmax, rl_num_in_x)
		call calc_mustat(refDx, 1.00_dl, drmin, drmax, rl_num_in_x)
	endif
	
	stop
	
	! Scanning different parameters; Check performance
	drmin = 50.0_dl
	drmax = 150.0_dl
	cubemin=40.0_dl
	cubemax=60.0_dl
	print *
	open(unit=1,file ='testomegam.txt')
	omegam = 0.00_dl
	do i = 1, 50
!		print *, omegam
		call tpCormlmu(omegam, w_dft, h_dft, cubemin, cubemax, numcube, &
			refDx, drmin, drmax, NORSDdiffmu, RSDdiffmu, numbins, NORSDmustat, RSDmustat)
		print *, omegam, sum(log(NORSDdiffmu)), sum(log(RSDdiffmu)) 
		write(1,'(<2*numcube+1>(e14.7,1x))') omegam, NORSDdiffmu, RSDdiffmu
		omegam = omegam + 0.02
	enddo
	close(1)
	
	drmin = 0.0_dl
	drmax = 150.0_dl
	cubemin=40.0_dl
	cubemax=60.0_dl
	print *
	open(unit=1,file ='testomegam_0to150.txt')
	omegam = 0.00_dl
	do i = 1, 50
		call tpCormlmu(omegam, w_dft, h_dft, cubemin, cubemax, numcube, &
			refDx, drmin, drmax, NORSDdiffmu, RSDdiffmu, numbins, NORSDmustat, RSDmustat)
		print *, omegam, sum(log(NORSDdiffmu)), sum(log(RSDdiffmu)) 
		write(1,'(<2*numcube+1>(e14.7,1x))') omegam, NORSDdiffmu, RSDdiffmu
		omegam = omegam + 0.02
	enddo
	close(1)
	
	drmin = 0.0_dl
	drmax = 80.0_dl
	cubemin=40.0_dl
	cubemax=60.0_dl
	print *
	open(unit=1,file ='testomegam_0to80.txt')
	omegam = 0.00_dl
	do i = 1, 50
		call tpCormlmu(omegam, w_dft, h_dft, cubemin, cubemax, numcube, &
			refDx, drmin, drmax, NORSDdiffmu, RSDdiffmu, numbins, NORSDmustat, RSDmustat)
		print *, omegam, sum(log(NORSDdiffmu)), sum(log(RSDdiffmu)) 
		write(1,'(<2*numcube+1>(e14.7,1x))') omegam, NORSDdiffmu, RSDdiffmu
		omegam = omegam + 0.02
	enddo
	close(1)

	drmin = 100.0_dl
	drmax = 150.0_dl
	print *
	open(unit=1,file ='testomegam_100to150.txt')
	omegam = 0.00_dl
	do i = 1, 50
		call tpCormlmu(omegam, w_dft, h_dft, cubemin, cubemax, numcube, &
			refDx, drmin, drmax, NORSDdiffmu, RSDdiffmu, numbins, NORSDmustat, RSDmustat)
		print *, omegam, sum(log(NORSDdiffmu)), sum(log(RSDdiffmu)) 
		write(1,'(<2*numcube+1>(e14.7,1x))') omegam, NORSDdiffmu, RSDdiffmu
		omegam = omegam + 0.02
	enddo
	close(1)
	
	drmin = 150.0_dl
	drmax = 200.0_dl
	cubemin=40.0_dl
	cubemax=60.0_dl
	print *
	open(unit=1,file ='testomegam_150to200.txt')
	omegam = 0.00_dl
	do i = 1, 50
		call tpCormlmu(omegam, w_dft, h_dft, cubemin, cubemax, numcube, &
			refDx, drmin, drmax, NORSDdiffmu, RSDdiffmu, numbins, NORSDmustat, RSDmustat)
		print *, omegam, sum(log(NORSDdiffmu)), sum(log(RSDdiffmu)) 
		write(1,'(<2*numcube+1>(e14.7,1x))') omegam, NORSDdiffmu, RSDdiffmu
		omegam = omegam + 0.02
	enddo
	close(1)

	drmin = 50.0_dl
	drmax = 150.0_dl
	cubemin=40.0_dl
	cubemax=50.0_dl
	print *
	open(unit=1,file ='testomegam_cube40to50.txt')
	omegam = 0.00_dl
	do i = 1, 50
		call tpCormlmu(omegam, w_dft, h_dft, cubemin, cubemax, numcube, &
			refDx, drmin, drmax, NORSDdiffmu, RSDdiffmu, numbins, NORSDmustat, RSDmustat)
		print *, omegam, sum(log(NORSDdiffmu)), sum(log(RSDdiffmu)) 
		write(1,'(<2*numcube+1>(e14.7,1x))') omegam, NORSDdiffmu, RSDdiffmu
		omegam = omegam + 0.02
	enddo
	close(1)

	drmin = 50.0_dl
	drmax = 150.0_dl
	cubemin=50.0_dl
	cubemax=60.0_dl
	print *
	open(unit=1,file ='testomegam_cube50to60.txt')
	omegam = 0.00_dl
	do i = 1, 50
		call tpCormlmu(omegam, w_dft, h_dft, cubemin, cubemax, numcube, &
			refDx, drmin, drmax, NORSDdiffmu, RSDdiffmu, numbins, NORSDmustat, RSDmustat)
		print *, omegam, sum(log(NORSDdiffmu)), sum(log(RSDdiffmu)) 
		write(1,'(<2*numcube+1>(e14.7,1x))') omegam, NORSDdiffmu, RSDdiffmu
		omegam = omegam + 0.02
	enddo
	close(1)

	drmin = 50.0_dl
	drmax = 150.0_dl
	cubemin=60.0_dl
	cubemax=70.0_dl
	print *
	open(unit=1,file ='testomegam_cube60to70.txt')
	omegam = 0.00_dl
	do i = 1, 50
		call tpCormlmu(omegam, w_dft, h_dft, cubemin, cubemax, numcube, &
			refDx, drmin, drmax, NORSDdiffmu, RSDdiffmu, numbins, NORSDmustat, RSDmustat)
		print *, omegam, sum(log(NORSDdiffmu)), sum(log(RSDdiffmu)) 
		write(1,'(<2*numcube+1>(e14.7,1x))') omegam, NORSDdiffmu, RSDdiffmu
		omegam = omegam + 0.02
	enddo
	close(1)

	halo_data_file = '../../../data/input/HR3LC1_57w_600to1787.dat'
	drmin = 50.0_dl
	drmax = 150.0_dl
	cubemin=40.0_dl
	cubemax=60.0_dl
	print *
	open(unit=1,file ='testomegam_LC1.txt')
	omegam = 0.00_dl
	do i = 1, 50
		call tpCormlmu(omegam, w_dft, h_dft, cubemin, cubemax, numcube, &
			refDx, drmin, drmax, NORSDdiffmu, RSDdiffmu, numbins, NORSDmustat, RSDmustat)
		print *, omegam, sum(log(NORSDdiffmu)), sum(log(RSDdiffmu)) 
		write(1,'(<2*numcube+1>(e14.7,1x))') omegam, NORSDdiffmu, RSDdiffmu
		omegam = omegam + 0.02
	enddo
	close(1)
	
	halo_data_file = '../../../data/input/HR3halo100w_pos1000.dat'
	drmin = 50.0_dl
	drmax = 150.0_dl
	cubemin=40.0_dl
	cubemax=60.0_dl
	print *
	open(unit=1,file ='testomegam_HR3A.txt')
	omegam = 0.00_dl
	do i = 1, 50
		call tpCormlmu(omegam, w_dft, h_dft, cubemin, cubemax, numcube, &
			refDx, drmin, drmax, NORSDdiffmu, RSDdiffmu, numbins, NORSDmustat, RSDmustat)
		print *, omegam, sum(log(NORSDdiffmu)), sum(log(RSDdiffmu)) 
		write(1,'(<2*numcube+1>(e14.7,1x))') omegam, NORSDdiffmu, RSDdiffmu
		omegam = omegam + 0.02
	enddo
	close(1)
		
	stop
	
contains

	subroutine calc_mustat(refDx, omegam, drmin, drmax, rl_num_in_x)
		!dummy
		real(dl), intent(in) :: refDx, omegam, rl_num_in_x, drmin, drmax
		!local
		character(len=char_len) :: str, file1, file2
		integer, parameter :: n = 10000, numchisq = 10000
		integer :: i, countmu1, countmu2, RSD
		real(dl) :: mu1,mu2, diffmu,diffmu1,diffmu2,diffmu3, muchisq1,muchisq2,muchisq3, sumdiffmu,sumlogdiffmu, &	
			deltax1,deltax2,deltaratio, &
			sumdiffmuA,sumlogdiffmuA, sumdiffmuB,sumlogdiffmuB
		integer, parameter :: numcube = 40
		real(dl) :: w, h, cubemin, cubemax, NORSDdiffmu(numcube), RSDdiffmu(numcube)
		integer :: numbins 
		type(mustatype) :: mustat, RSDmustat
	
		cubemin=40.0_dl
		cubemax=60.0_dl
		numbins = 10

		call tpCormu(omegam, w_dft, h_dft, rl_num_in_x, refDx, drmin, drmax, diffmu1, diffmu2, .true., numbins, mustat, RSDmustat)
	
		print *, 'diffmu = ', diffmu1, diffmu2
		write(str,'(f5.2)') omegam
		file1 = '../testmuhist/om'//trim(adjustl(str))//'_muhist.txt'
		file2 = '../testmuhist/om'//trim(adjustl(str))//'_RSDmuhist.txt'	
		open(unit=1,file=file1)
		open(unit=2,file=file2)
		write(1,'(<numbins>(e14.7,1x))') dble(mustat%nummuarray(1:numbins))
		write(1,'(<numbins>(e14.7,1x))') mustat%weimuarray(1:numbins)
		write(2,'(<numbins>(e14.7,1x))') dble(RSDmustat%nummuarray(1:numbins))
		write(2,'(<numbins>(e14.7,1x))') RSDmustat%weimuarray(1:numbins)
		close(1)
		close(2)
	end subroutine calc_mustat	
	
end program ap_main_test
