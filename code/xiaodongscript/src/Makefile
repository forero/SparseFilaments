EXE = ../getcor ../select_obj ../post_proc
OBJS = getcor.o select_obj.o post_proc.o
F90C = ifort
F90FLAGS = 

default: $(EXE)

$(EXE): $(OBJS)
	$(F90C) -o ../getcor getcor.o
	$(F90C) -o ../select_obj select_obj.o
	$(F90C) -o ../post_proc post_proc.o

%.o: %.f90
	$(F90C) -c $*.f90 $(F90FLAGS)

clean :
	rm -rf *.o *.mod $(EXE)
