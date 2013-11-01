EXE = ../getcor ../select_obj ../post_proc ../post_proc_BDMW ../calc_grad_chisqs ./test
OBJS = getcor.o select_obj.o post_proc.o post_proc_BDMW.o ap_tools.o  ap_mu_tools.o ap_cosmo_funs.o ap_settings_init.o ap_smooth.o ap_grids.o ap_grad_fields.o ap_chisq.o ap_main_calc_grad_chisqs.o ap_main_test.o
F90C = ifort -heap-arrays
F90FLAGS = 

default: $(EXE)

$(EXE): $(OBJS)
	$(F90C) -o ../getcor getcor.o
	$(F90C) -o ../select_obj select_obj.o
	$(F90C) -o ../post_proc post_proc.o
	$(F90C) -o ../post_proc_BDMW post_proc_BDMW.o
#	$(F90C) -o ../scan_chisqs ap_tools.o ap_mu_tools.o ap_cosmo_funs.o ap_settings_init.o ap_smooth.o ap_grad_fields.o ap_chisq.o ap_main_calc_grad_chisqs2.o
	$(F90C) -o ../calc_grad_chisqs ap_tools.o ap_mu_tools.o ap_cosmo_funs.o ap_settings_init.o ap_smooth.o ap_grad_fields.o ap_chisq.o ap_main_calc_grad_chisqs.o

	$(F90C) -o ./test ap_tools.o ap_mu_tools.o ap_cosmo_funs.o ap_settings_init.o ap_smooth.o ap_grids.o ap_grad_fields.o ap_chisq.o ap_main_test.o  
#	$(F90C) -o ../scan_chisqB ap_tools.o ap_mu_tools.o ap_cosmo_funs.o ap_settings_init.o ap_smooth.o ap_grad_fields.o ap_chisq.o ap_main_scan_chisq_B.o

ap_cosmo_funs.o: ap_tools.o
ap_mu_tools.o: ap_tools.o
ap_settings_init.o: ap_cosmo_funs.o 
ap_smooth.o: ap_settings_init.o
ap_grids.o: ap_smooth.o
ap_grad_fields.o: ap_smooth.o ap_mu_tools.o
ap_chisq.o: ap_grad_fields.o
ap_main_calc_grad_chisqs.o:  ap_chisq.o
ap_main_test.o: ap_chisq.o ap_grids.o
#ap_main.o: ap_settings_init.o
#ap_main.o: ap_tools.o
#ap_main.o: ap_cosmo_funs.o


%.o: %.f90
	$(F90C) -c -heap-arrays $*.f90 $(F90FLAGS)

clean :
	rm -rf *.o *.mod $(EXE)
