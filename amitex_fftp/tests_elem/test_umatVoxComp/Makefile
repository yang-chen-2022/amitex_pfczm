lAmitex  =../../libAmitex

#amitex_fftp-config : script permettant de connaitre l'emplacement des modules (includes) et des librairies
AMITEX_FFTP_CONFIG = $(lAmitex)/bin/amitex_fftp-config

# Amitex flags
F90FLAGS = `$(AMITEX_FFTP_CONFIG) --fcflags` 

LIBS = `$(AMITEX_FFTP_CONFIG) --plugin-libs`

F90 = `$(AMITEX_FFTP_CONFIG) --mpifc`

#debug gcc
#DEBUG=-Wall

#debug intel
#DEBUG=-g -check all -fpe0 -warn -traceback -debug extended
#DEBUG=-warn all
DEBUG=


OBJ1 = test_elasiso.o
OBJ2 = test_umatvoigt.o
OBJ3 = test_umatreuss.o
OBJ4 = test_CC.o
OBJ5 = test_umatvoigtCC.o
OBJ6 = test_umatreussCC.o

all: test_elasiso test_umatvoigt test_umatreuss test_CC test_umatvoigtCC test_umatreussCC

test_elasiso: $(OBJ1)
	$(F90) -o $@ $(OBJ1) $(LIBS) $(FoX) $(DEBUG)

test_umatvoigt: $(OBJ2)
	$(F90) -o $@ $(OBJ2) $(LIBS) $(FoX) $(DEBUG)

test_umatreuss: $(OBJ3)
	$(F90) -o $@ $(OBJ3) $(LIBS) $(FoX) $(DEBUG)

test_CC : $(OBJ4)
	$(F90) -o $@ $(OBJ4) $(LIBS) $(FoX) $(DEBUG)

test_umatvoigtCC: $(OBJ5)
	$(F90) -o $@ $(OBJ5) $(LIBS) $(FoX) $(DEBUG)

test_umatreussCC: $(OBJ6)
	$(F90) -o $@ $(OBJ6) $(LIBS) $(FoX) $(DEBUG)

clean:
	rm -f *.o test_elasiso test_umatvoigt test_umatreuss test_CC test_umatvoigtCC test_umatreussCC

%.o : %.f90
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<
