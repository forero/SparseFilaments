!############################3333333333333333333333
!
! definitions, functions related with cosmology
!
!############################3333333333333333333333

module ap_cosmo_funs
use ap_tools
implicit none

	integer,  parameter	:: de_num_intpl = 1500

	real(dl),  parameter	:: de_basenumber = 1.0_dl + 1.0_dl/512.0_dl
	real(dl),  parameter	:: de_logbasenumber = log(de_basenumber)
	real(dl) :: de_maxintplz

	
	!Common used array, saving the interpolating data.
	! Frequently used in many dark energy models, so put them here as common variables.
	! Array for redshift.
	real(dl) :: de_zdata(de_num_intpl)
	! Array for ez, defined as H(z)/H(0)
	real(dl) :: de_comovr_data(de_num_intpl)

	logical :: cosmo_funs_inited = .false.

	real(dl) :: gb_omegam, gb_w, gb_h

contains

  !------------------------------------------
  ! Initialization for cosmo_funs
  !------------------------------------------
	subroutine cosmo_funs_init()
		integer :: i
		
		write(*,*) "Initializing cosmo_funs..."
		write(*,*) " Redshift interpolating numbers = ", de_num_intpl

		!Initialize the redshift data
		do i = 1, de_num_intpl
			de_zdata(i) = de_zi(i)
		enddo
		
		!Range of the interpolation (maximal redshift)
		de_maxintplz = de_zdata(de_num_intpl)
		write(*,*) " Maximal redshift in interpolating = ", de_maxintplz
				
		cosmo_funs_inited = .true.
	end subroutine cosmo_funs_init

  !------------------------------------------
  ! Constructing the comovr table
  !------------------------------------------
  	subroutine de_calc_comovr()
  		integer :: i
  		de_comovr_data(1) = 0
  		do i = 2, de_num_intpl
  			de_comovr_data(i) = de_comovr_data(i-1) + &
		  		simpson(inv_Hz,de_zdata(i-1),de_zdata(i),1)*const_c*gb_h
		enddo
	end subroutine de_calc_comovr
  !------------------------------------------
  ! get comov_r from interploation
  !------------------------------------------
  	real(dl) function de_get_comovr(z)
  		integer :: i, i1, i2
  		real(dl) :: z
  		i = de_iz(z)
  		call ilist(i,1,1,de_num_intpl,i1,i2)
  		de_get_comovr = intpl_vl(z, de_zdata(i1), de_comovr_data(i1), &
  			de_zdata(i1+1), de_comovr_data(i1+1),&
  			de_zdata(i2), de_comovr_data(i2))
	end function de_get_comovr

  !------------------------------------------
  ! Get z from index
  !------------------------------------------
	real(dl) function de_zi(i)
		integer :: i
		de_zi = de_basenumber**DBLE(i-1) - 1.0_dl
	end function de_zi

  !------------------------------------------
  ! Get index from z
  !------------------------------------------
	integer function de_iz(z)
		real(dl) :: z, temp
		de_iz = Ceiling(log(1.0_dl+z)/de_logbasenumber + 1.0_dl)
		if(de_iz < 2) de_iz = 2
	end function de_iz

  !------------------------------------------
  ! Hubble parameter in unit of km/s/Mpc
  !------------------------------------------	
  	real(dl) function Hz(z)
  		real(dl) :: z
  		Hz = 100.0*gb_h*&
  		sqrt(gb_omegam*(1.0+z)**3.0 + (1.0-gb_omegam)*(1.0+z)**(3.0*(1.0+gb_w)))
  	end function Hz
  	real(dl) function inv_Hz(z)
  		real(dl) :: z
  		inv_Hz = 100.0*gb_h*sqrt(gb_omegam*(1.0+z)**3.0 + (1.0-gb_omegam)*(1.0+z)**(3.0*(1.0+gb_w)))
  		inv_Hz = 1.0 / inv_Hz
  	end function inv_Hz
  	
  !------------------------------------------
  ! comoving r in unit of Mpc/h
  !------------------------------------------	
	real(dl) function comov_r(z)
  		real(dl) :: z
  		integer :: n
  		n = max(4*ceiling(z / 0.125),4)
  		comov_r = simpson(inv_Hz,0.0_dl,z,n)*const_c*gb_h
  	end function comov_r

  !------------------------------------------
  !   	 get z according to comoving_r 
  !------------------------------------------	
	real(dl) function get_z(gv_comov_r, gv_zl, gv_zr)
		real(dl) :: gv_comov_r, zl, zr
		real(dl), optional :: gv_zl, gv_zr
		if(present(gv_zl)) then
			zl = gv_zl
			else
			zl = 0.0
		endif
		if(present(gv_zr)) then
			zr = gv_zr
			else
			zr = 3.0
		endif
		get_z = findroot(comov_r, gv_comov_r, zl, zr)
	end function get_z
end module ap_cosmo_funs

