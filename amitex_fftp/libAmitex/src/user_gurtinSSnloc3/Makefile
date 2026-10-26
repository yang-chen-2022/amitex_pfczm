#==================================================================
# ADJUST : the name of the executable : MAIN 
#
# If you want to add new own modules :
#    add them in MOD_user
#    add a rule in the Makefile for compiling the new added module
#==================================================================
#
# Use amitex_fftp-config tool :
# -----------------------------
# To access help, once AMITEX installed :
#    - source (amitex)/env_amitex.sh
#    - amitex_fftp-config --help
#
#=============================================================

lAmitex=$(AMITEX_PATH)/libAmitex

#=============================================================
#                                                    VARIABLES

# Executable name
MAIN = amitex_fftp

# "User" modules 
MOD_user = non_local_user_mod.o

#amitex_fftp-config : script permettant de connaitre l'emplacement des modules (includes) et des librairies
AMITEX_FFTP_CONFIG = $(lAmitex)/bin/amitex_fftp-config

# Amitex flags
F90FLAGS = `$(AMITEX_FFTP_CONFIG) --fcflags` 

LIBS = `$(AMITEX_FFTP_CONFIG) --plugin-libs`

F90 = `$(AMITEX_FFTP_CONFIG) --mpifc`

#=============================================================
#                                           LINK -> EXECUTABLE
all: $(MAIN) 

$(MAIN): $(MOD_user)
	$(F90) -o $@ $(MOD_user) $(LIBS)

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

non_local_user_mod.o: non_local_user_mod.f90 
	$(F90) $(F90FLAGS) $(DEBUG) $(OPTIM) -c $<

#=============================================================
#                                                        CLEAN
# supprime les .o .a et .mod
clean:
	rm -f *.o *.mod $(MAIN)



