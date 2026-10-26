lAmitex = ../../libAmitex

#amitex_fftp-config : script permettant de connaitre l'emplacement des modules (includes) et des librairies
AMITEX_FFTP_CONFIG = $(lAmitex)/bin/amitex_fftp-config

# Amitex flags
F90FLAGS = `$(AMITEX_FFTP_CONFIG) --fcflags` 

LIBS = `$(AMITEX_FFTP_CONFIG) --plugin-libs`

F90 = `$(AMITEX_FFTP_CONFIG) --mpifc`

OBJ = test_field_divgrad.o

OPTIM = -O3

#----- DEBUG GCC
#DEBUG=-Wall -g -fcheck=all
#DEBUG=-Wall -g -fcheck=all -fbounds-check -ffpe-trap=zero,overflow

#----- DEBUG INTEL
DEBUG=-g -check all -fpe0 -warn all,noexternal -traceback -debug extended
DEBUG=-warn all
DEBUG=

test_field_divgrad: $(OBJ)
	$(F90) -o $@ $(DEBUG) $(OBJ) $(LIBS)     
#        $(F90) -o $@ $(STD) $(ALGO) $(LIBS)
#        $(F90) -o $@ $(STD) $(ALGO) $(LIBS) $(OMP)
clean:
	rm -f *.o test_field_divgrad *.log *.vtk

%.o : %.f90
	$(F90) $(F90FLAGS) $(DEBUG) $(OPTIM) -c $<
