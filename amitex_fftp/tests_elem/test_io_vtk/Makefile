lAmitex  =../../libAmitex

#amitex_fftp-config : script permettant de connaitre l'emplacement des modules (includes) et des librairies
AMITEX_FFTP_CONFIG = $(lAmitex)/bin/amitex_fftp-config

# Amitex flags
F90FLAGS = `$(AMITEX_FFTP_CONFIG) --fcflags` 

LIBS = `$(AMITEX_FFTP_CONFIG) --plugin-libs`

F90 = `$(AMITEX_FFTP_CONFIG) --mpifc`

OBJ = test_io_vtk.o

OPTIM = -O3

#DEBUG=-Wall -g -fcheck=all
DEBUG=
#Intel debug options
DEBUG=-g -check all -fpe0 -warn all,noexternal -traceback -debug extended
DEBUG=-warn all
DEBUG=

test_io_vtk: $(OBJ)
	$(F90) -o $@ $(OBJ) $(LIBS)

clean:
	rm -f *.o *.vtk test_io_vtk *.log

%.o : %.f90
	$(F90) $(F90FLAGS) $(DEBUG) $(OPTIM) -c $<

