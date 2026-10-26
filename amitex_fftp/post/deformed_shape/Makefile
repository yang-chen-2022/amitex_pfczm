#==================================================================
# ADJUST : the name of the executable : MAIN 
#
# If you want to add new own modules :
#    add them in MOD_user
#    add a rule in the Makefile for compiling the new added module
#==================================================================
#
# amitex_fftp-config tool :
# -------------------------
# To access help, once AMITEX installed :
#    - source (amitex)/env_amitex.sh
#    - amitex_fftp-config --help
#
#=======================================================================

# Path to amitex (path_to_amitex_fftp/libAmitex)
lAmitex=$(AMITEX_PATH)/libAmitex

#=============================================================
#                                                    VARIABLES

# Executable name
MAIN = deformedShape

# "User" modules 
MOD_user = deformed_shape_mod.o 

#amitex_fftp-config : script permettant de connaitre l'emplacement des modules (includes) et des librairies
AMITEX_FFTP_CONFIG = $(lAmitex)/bin/amitex_fftp-config

# Amitex flags
F90FLAGS = `$(AMITEX_FFTP_CONFIG) --fcflags` 

LIBS = `$(AMITEX_FFTP_CONFIG) --plugin-libs`

F90 = `$(AMITEX_FFTP_CONFIG) --mpifc`

#=============================================================
#                                           LINK -> EXECUTABLE
all: $(MAIN) 

$(MAIN): $(MAIN).o
	$(F90) -o $@ $(MAIN).o $(MOD_user) $(LIBS)

#=============================================================
#                                         OPTIMISATION OPTIONS
OPTIM = -O3

#=============================================================
#                                                DEBUG OPTIONS
#debug intel
#DEBUG=-g -check all -fpe0 -warn -traceback -debug extended
#DEBUG=-g -check all
#DEBUG=-warn all

#debug gcc
#DEBUG=-Wall
#DEBUG=-g -fcheck=all
#
DEBUG=

#=============================================================
#                                                 COMPILATIONS

$(MAIN).o : $(MAIN).f90 deformed_shape_mod.o
	$(F90) $(F90FLAGS) $(DEBUG) $(OPTIM) -c $<

%.o : %.f90
	$(F90) $(F90FLAGS) $(DEBUG) $(OPTIM) -c $<

#=============================================================
#                                                        CLEAN
# supprime les .o .a et .mod
clean:
	rm -f *.o *.mod *.log $(MAIN)



