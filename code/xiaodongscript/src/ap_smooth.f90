

!####################################
!This module does smooth
!####################################
module ap_smooth
use ap_cell

	implicit none

!!! IMPORTANT SETTINGS
! settings of the programe	
	real(dl) :: dft_ra_ratio = 4.0d0 ! 4.0d0 TESTING
	integer, public :: gb_min_smnum = 0 ! minimal 4 (gb_min_smnum + 1) halos to get gradient; .le. gb_min_smnum will be skipped
	logical, public :: gb_keepzerorho = .true. ! do not skip pixels with num_halo <= gb_min_smnum but keep them (may set rho as 0)

	
contains

  !------------------------------------------
  ! N nearest neighbours selection
  ! select out halos in the cells nearby the
  !  given position. make sure the number is 
  !  far larger than the required
  !------------------------------------------
  	subroutine nb_select(x,y,z,num,given_ra_ratio,selected_list)
  		real(dl), intent(in) :: x,y,z
  		real(dl), optional, intent(in) :: given_ra_ratio
  		integer, intent(in) :: num
  		integer, allocatable, intent(out) :: selected_list(:)
  		integer :: i,j,k,l1,l2,ix,iy,iz,requirednum
  		integer :: di, i1,i2, j1,j2, k1,k2, now_num
  		real(dl) :: dev, ra_ratio, ratio
  		ix = max(min(int((x-gbgridxmin)/gbdeltax)+1,gb_n_cellx),1)
  		iy = max(min(int((y-gbgridymin)/gbdeltay)+1,gb_n_celly),1)
  		iz = max(min(int((z-gbgridzmin)/gbdeltaz)+1,gb_n_cellz),1)
  		dev = max(abs(abs(gbgridxmin + (ix-1)*gbdeltax - x)/gbdeltax - 0.5), &
  			abs(abs(gbgridymin + (iy-1)*gbdeltay - y)/gbdeltay - 0.5), &
  			abs(abs(gbgridzmin + (iz-1)*gbdeltaz - z)/gbdeltaz - 0.5) )
  		if(present(given_ra_ratio)) then
  			ra_ratio = given_ra_ratio
  		else
  			ra_ratio = dft_ra_ratio
  		endif
  		ratio = ra_ratio * (1.0 + dev * 1.5)
  		requirednum = max(int(num*ratio), 40) ! at least 40 halos!!!
!  		print *, 'ix, iy, iz, dev, ratio, requirednum = ', ix, iy, iz, dev, ratio, requirednum
		di = 0
		do while(1.eq.1)
			call ilist(ix,di,1,gb_n_cellx,i1,i2)
			call ilist(iy,di,1,gb_n_celly,j1,j2)
			call ilist(iz,di,1,gb_n_cellz,k1,k2)
			now_num = 0
			do i = i1, i2
			do j = j1, j2
			do k = k1, k2
				now_num = now_num + gb_cell_mat(i,j,k)%halo_num
			enddo
			enddo
			enddo
			if(now_num .ge. requirednum) exit
			di = di + 1
		enddo
		allocate(selected_list(now_num))
		l1 = 1
		do i = i1, i2
		do j = j1, j2
		do k = k1, k2
			do l2 = 1, gb_cell_mat(i,j,k)%halo_num
				selected_list(l1) = gb_cell_mat(i,j,k)%list(l2)
				l1 = l1 + 1
			enddo
		enddo
		enddo
		enddo
	end subroutine nb_select

  !------------------------------------------
  ! estimating rho and gradient rho based on
  !  cubic spline kernel; fixed radius
  !------------------------------------------
  	subroutine nb_fixmd_list(x,y,z,num,rho,drhodx,drhody,drhodz,fixmd,dosegcut,segcutdist,touchbdflag,vcorect_er_flag)!Testing
!  	subroutine nb_fixmd_list(x,y,z,num,rho,drhodx,drhody,drhodz,dosegcut,touchbdflag)!Testing
		! DUMMY ARGUMENTS
  		real(dl), intent(in) :: x,y,z, fixmd,segcutdist
  		integer, intent(out) :: num
  		real(dl), intent(out) :: rho,drhodx,drhody,drhodz
  		logical, intent(in) :: dosegcut
  		logical, intent(out):: touchbdflag, vcorect_er_flag
  		! LOCAL VARIABLES
  		integer, parameter :: max_num = 100000 !maximal # of neary halos 
  		real(dl) :: distance_array(max_num), xyz_mass_array(4,max_num)
  		integer :: i,j,k,l, imin,imax,jmin,jmax,kmin,kmax, nownum,nowindex
  		real(dl) :: r0(3), h, nowr, mass, dweight

  		imin = int((x-fixmd-gbgridxmin)/gbdeltax +1.0)
  		imax = int((x+fixmd-gbgridxmin)/gbdeltax +1.0)
  		jmin = int((y-fixmd-gbgridymin)/gbdeltay +1.0)
  		jmax = int((y+fixmd-gbgridymin)/gbdeltay +1.0)
  		kmin = int((z-fixmd-gbgridzmin)/gbdeltaz +1.0)
  		kmax = int((z+fixmd-gbgridzmin)/gbdeltaz +1.0)
  		
		num = 0
		r0(1)=x; r0(2)=y; r0(3)=z;
		
		if(imin<1.or.imin>gb_n_cellx.or.jmin<1.or.jmax>gb_n_celly.or.kmin<1.or.kmax>gb_n_cellz) then
			touchbdflag = .true.
		else
			touchbdflag = .false.
		endif
		
		vcorect_er_flag = .false.
		do i = max(1,imin), min(gb_n_cellx,imax)
		do j = max(1,jmin), min(gb_n_celly,jmax)
		do k = max(1,kmin), min(gb_n_cellz,kmax)
!		print *, 'i,j,k=',i,j,k
!		print *, '# of halo = ', gb_cell_mat(i,j,k)%halo_num
			do l = 1, gb_cell_mat(i,j,k)%halo_num
				nowindex = gb_cell_mat(i,j,k)%list(l)
				nowr = distance(gb_xyz_list(1:3,nowindex),r0)
				if(nowr < fixmd) then
					if(dosegcut.and.nowr<segcutdist) cycle ! Apply seg_cut
					num = num+1
!					print *, 'num,max_num = ', num,max_num
					distance_array(num) = nowr
					xyz_mass_array(1:3, num) = gb_xyz_list(1:3,nowindex)
					xyz_mass_array(4,num) = gb_mass_list(nowindex)
					if(halo_info(nowindex)%vth_er) vcorect_er_flag = .true.
				endif
			enddo
		enddo
		enddo
		enddo

		if(num > max_num) then
			print *, 'ERROR (nb_fixmd_list): # of halos overflow: ', num, max_num
			stop
		endif

		rho = 0; drhodx=0;drhody=0;drhodz=0;		
		if(num .le. gb_min_smnum) then
			return
		endif			
			
  		h = fixmd / 2.0
		do i = 1, num
			nowr = distance_array(i)
			mass = xyz_mass_array(4,i)
			rho = rho + mass*w_kernel(nowr, h)
			dweight = der_w_kernel(nowr,h)
!			print *, 'i,nowr,mass,rho,weight,dweight = ', i,nowr,mass,rho,w_kernel(nowr,h),dweight
			drhodx = drhodx + mass*(x-xyz_mass_array(1,i)) / nowr * dweight
			drhody = drhody + mass*(y-xyz_mass_array(2,i)) / nowr * dweight
			drhodz = drhodz + mass*(z-xyz_mass_array(3,i)) / nowr * dweight
		enddo
!		stop !TESTING
	end subroutine nb_fixmd_list


  !------------------------------------------
  ! estimating rho and gradient rho based on
  !  cubic spline kernel; fixed radius
  !------------------------------------------
  	subroutine vth_fixmd(x,y,z,num,vxth,vyth,vzth,fixmd,dosegcut,segcutdist,touchbdflag,vcorect_er_flag,vxcomp,vycomp,vzcomp)!Testing
		! DUMMY ARGUMENTS
  		real(dl), intent(in) :: x,y,z, fixmd,segcutdist
  		integer, intent(out) :: num
  		real(dl), intent(out) :: vxth,vyth,vzth
  		logical, intent(in) :: dosegcut
  		real(dl), optional, intent(out) :: vxcomp,vycomp,vzcomp
  		logical, intent(out):: touchbdflag,vcorect_er_flag ! will be true if touches boundary or vcorert er (indicating error)!
  		! LOCAL VARIABLES
  		integer, parameter :: max_num = 100000 !maximal # of neary halos 
  		real(dl) :: distance_array(max_num), xyz_mass_array(4,max_num)
  		integer :: ix,iy,iz, i,j,k,l, imin,imax,jmin,jmax,kmin,kmax, nownum,nowindex, indexarray(max_num)
  		real(dl) :: r0(3), h, nowr,nowrcube, mass, weight,sumweight, r,redshift,Hubble,growthfactor,frac,vol,dvol

		call cell_index(ix,iy,iz,x,y,z)
		
  		imin = int((x-fixmd-gbgridxmin)/gbdeltax +1.0_dl)
  		imax = int((x+fixmd-gbgridxmin)/gbdeltax +1.0_dl)
  		jmin = int((y-fixmd-gbgridymin)/gbdeltay +1.0_dl)
  		jmax = int((y+fixmd-gbgridymin)/gbdeltay +1.0_dl)
  		kmin = int((z-fixmd-gbgridzmin)/gbdeltaz +1.0_dl)
  		kmax = int((z+fixmd-gbgridzmin)/gbdeltaz +1.0_dl)
  		
		num = 0
		r0(1)=x; r0(2)=y; r0(3)=z;
		
		if(imin<1.or.imin>gb_n_cellx.or.jmin<1.or.jmax>gb_n_celly.or.kmin<1.or.kmax>gb_n_cellz) then
			touchbdflag = .true.
		else
			touchbdflag = .false.
		endif
		
		vcorect_er_flag = .false.
		do i = max(1,imin), min(gb_n_cellx,imax)
		do j = max(1,jmin), min(gb_n_celly,jmax)
		do k = max(1,kmin), min(gb_n_cellz,kmax)
			do l = 1, gb_cell_mat(i,j,k)%halo_num
				nowindex = gb_cell_mat(i,j,k)%list(l)
				nowr = distance(gb_xyz_list(1:3,nowindex),r0)
				if(nowr < fixmd) then
					if(dosegcut.and.nowr<segcutdist) cycle 
					num = num+1
					distance_array(num) = nowr
					indexarray(num) = nowindex
					xyz_mass_array(1:3, num) = gb_xyz_list(1:3,nowindex)
					xyz_mass_array(4,num) = gb_mass_list(nowindex)
					if(halo_info(nowindex)%vth_er) vcorect_er_flag = .true.
				endif
			enddo
		enddo
		enddo
		enddo

		if(num > max_num) then
			print *, 'ERROR (nb_fixmd_list): # of halos overflow: ', num, max_num
			stop
		endif

		vxth=0;vyth=0;vzth=0;		
		if(num .le. gb_min_smnum) then
			return
		endif			
		if(present(vxcomp)) then
			vxcomp=0;vycomp=0;vzcomp=0;sumweight=0.0
		endif
  		h = fixmd / 2.0_dl
		vol = (4.0_dl*const_pi)/3.0_dl*fixmd**3.0
		dvol = vol/dble(num)
		do i = 1, num
			nowr = distance_array(i)
			if(nowr < 1.0e-4) cycle
			nowrcube = nowr**3.0
			mass = xyz_mass_array(4,i)
			weight = w_kernel(nowr, h)
			vxth = vxth + (mass-dvol)*(xyz_mass_array(1,i)-x) / nowrcube
			vyth = vyth + (mass-dvol)*(xyz_mass_array(2,i)-y) / nowrcube
			vzth = vzth + (mass-dvol)*(xyz_mass_array(3,i)-z) / nowrcube
!			print *, 'i,nowr,vx,vy,vz,mass,weight,nowrcube=',i,nowr,vx,vy,vz,mass,weight,nowrcube
			if(present(vxcomp)) then
				vxcomp = vxcomp + gb_vel_list(1,indexarray(i))*weight
				vycomp = vycomp + gb_vel_list(2,indexarray(i))*weight
				vzcomp = vzcomp + gb_vel_list(3,indexarray(i))*weight
				sumweight = sumweight+weight
			endif
!			print *, 'i,nowr,mass,rho,weight,dweight = ', i,nowr,mass,rho,w_kernel(nowr,h),dweight
		enddo
		if(present(vxcomp)) then
			vxcomp = vxcomp / sumweight
			vycomp = vycomp / sumweight
			vzcomp = vzcomp / sumweight
		endif

		r = sqrt(x*x+y*y+z*z)
		redshift = de_zfromintpl(dble(r))
		Hubble = Hz(dble(redshift))
		growthfactor = de_gfz_intpl(dble(redshift))
		frac = Hubble * growthfactor / (4.0_dl*const_pi) / (1.0_dl+redshift)!*vol
!		write(*,'(A,8(f10.3,1x))') 'redshift, Hubble, growthfactor, frac, vx,vy,vz = ', &
!			redshift, Hubble, growthfactor, frac, vx,vy,vz 
!		print *, 'vx,vy,vz,frac = ', vx,vy,vz,frac
		vxth = vxth*frac
		vyth = vyth*frac
		vzth = vzth*frac
!		stop !TESTING
	end subroutine vth_fixmd


  !------------------------------------------
  ! estimating rho and gradient rho based on
  !  cubic spline kernel
  !------------------------------------------
  	subroutine nb_list0(x,y,z,num,rho,drhodx,drhody,drhodz,max_dist)
		! DUMMY ARGUMENTS
  		real(dl), intent(in) :: x,y,z
  		integer, intent(in) :: num
  		real(dl), intent(out) :: rho,drhodx,drhody,drhodz,max_dist
  		! LOCAL VARIABLES
  		integer, allocatable :: index_array(:)
  		real(dl), allocatable :: distance_array(:), xyz_mass_array(:,:), tmp(:), smda(:)
  		integer, allocatable :: tmpindex(:), smlablist(:)
  		integer :: i,j,n,nowindex
  		real(dl) :: h, r0(3), r, mass, dweight

  		call nb_select(x,y,z,num,selected_list=index_array)
  		n=size(index_array)

  		allocate(tmpindex(n), distance_array(n))
		r0(1)=x; r0(2)=y; r0(3)=z;
		do i = 1, n
			distance_array(i) = distance(gb_xyz_list(1:3,index_array(i)),r0)
			tmpindex(i) = index_array(i)
		enddo

		allocate(smlablist(num),smda(num))
		call ltlablist2(distance_array,n,num,smda,smlablist)

		allocate(xyz_mass_array(4,num))
		
		do i = 1, num
			nowindex = tmpindex(smlablist(i))
			xyz_mass_array(1:3,i) = gb_xyz_list(1:3,nowindex)
			xyz_mass_array(4,i) = gb_mass_list(nowindex)
!		  		print *, i, real(xyz_mass_array(1:5,i))
 		enddo

  		max_dist = maxval(smda)
  		h = max_dist / 2.0
		rho = 0; drhodx=0;drhody=0;drhodz=0;
		do i = 1, num
			r = distance_array(smlablist(i))
			mass = xyz_mass_array(4,i)
!			print *, i, mass !TESTING			
			rho = rho + mass*w_kernel(r, h)
			dweight = der_w_kernel(r,h)
			drhodx = drhodx + mass*(x-xyz_mass_array(1,i)) / r * dweight
			drhody = drhody + mass*(y-xyz_mass_array(2,i)) / r * dweight
			drhodz = drhodz + mass*(z-xyz_mass_array(3,i)) / r * dweight
		enddo
!		stop !TESTING
	end subroutine nb_list0


	
end module ap_smooth	
	

