#=============================================================
#
# Use amitex_fftp-config tool :
# -----------------------------
# To access help, once AMITEX installed :
#    - source (amitex)/env_amitex.sh
#    - amitex_fftp-config --help
#
#=============================================================

# Chemin vers amitex
lAmitex=$(AMITEX_PATH)/libAmitex

#amitex_fftp-config : script permettant de connaitre l'emplacement des modules (includes) et des librairies
AMITEX_FFTP_CONFIG = $(lAmitex)/bin/amitex_fftp-config

#acces aux interfaces des modules fortran de 2decompfft et amitex (field_mod)
INCLUDE = `$(AMITEX_FFTP_CONFIG) --include` 

# option de compilation d'AMITEX
F90FLAGS = `$(AMITEX_FFTP_CONFIG) --fcflags`

#debug options
DEBUG=

#pour edition de lien du cas-test
LIBS = `$(AMITEX_FFTP_CONFIG) --plugin-libs`

F90 = `$(AMITEX_FFTP_CONFIG) --mpifc` 
CPP=g++
#CPP=icpc

all:
	$(CPP) -std=c++17 -c caiFilters.cpp
	$(F90) -c caifilters_mod.f90 $(INCLUDE) $(F90FLAGS) $(DEBUG)
	$(F90) -c test_field_filterF.f90 $(INCLUDE) $(F90FLAGS) $(DEBUG)
	$(F90) -o test_field_filterF *.o $(LIBS) -lstdc++  -ldl
	
example:
	$(CPP) -std=c++17 -c caiFilters.cpp
	$(F90) -c caifilters_mod.f90 $(INCLUDE) $(F90FLAGS) $(DEBUG)
	$(F90) -c example.f90
	$(F90) -o example *.o $(LIBS) $(FoX) -lstdc++  -ldl

clean:
	rm *.o *.mod example test_field_filterF *.log *.vtk
