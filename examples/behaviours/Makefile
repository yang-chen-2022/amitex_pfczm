# Compilateur 
FC=gfortran

# Chemin vers la librairie 2deomp_fft
lDecomp=../../amitex_fftp/lib_extern/2decomp_fft-1.5.847

# options :
# NO_US pour supprime les ajout de sous-tirets (incompatible avec l'appel via umat_read)
# FASTMATH (pour gnu) accelere certaines fonctions... mais ATTENTION, comportement incertain.
# DEF_INT : pour assurer la compatibilite avec l'interface umat d'mfront et de cast3m,
#           les entiers par defaut sont codes sur 8 bits

FASTMATHgnu=-ffast-math
ifeq ($(FC),ifort)
	NO_US=-assume nounderscore
	FASTMATH=
        DEF_INT=-i8
endif
ifeq ($(FC),gfortran)
	NO_US=-fno-underscoring
	FASTMATH=$(FASTMATHgnu)
        DEF_INT=-fdefault-integer-8
endif

## On reprend le fichier Makefile.inc de la librairie 2decomp_fft
include $(lDecomp)/src/Makefile.inc

#donne accès aux modules fortran de 2decompfft
INCLUDE = -I$(lDecomp)/include

# fichiers sources et objets correspondants
SRC = $(wildcard *.f90) $(wildcard *.F)
OBJ=$(SRC:.f90=.o) $(SRC:.F=.o)


#compilation de la librairie libUmatAmitex.so
all : $(OBJ)
	$(F90) -o libUmatAmitex.so -shared *.o

#regle generique
# Compilation des sources (fortran)
%.o : %.f90
#	$(F90) $(INCLUDE) -O3 $(NO_US) -fPIC -DPIC $< -o $@ -c
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEF_INT) $(NO_US) $(FASTMATH) -fPIC -DPIC -c $<

# Compilation des sources (fortran)
%.o : %.F
#	$(F90) $(INCLUDE) -O3 $(NO_US) -fPIC -DPIC $< -o $@ -c
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEF_INT) $(NO_US) $(FASTMATH) -fPIC -DPIC -c $<

clean:
	rm -f *.o
	rm -f *.so

#===================================================================
# avant utilisation des librairies dynamiques :
#
## On reprend le fichier Makefile.inc de la librairie 2decomp_fft
#include $(lDecomp)/src/Makefile.inc
#
##donne accès aux modules fortran de 2decompfft
#INCLUDE = -I$(lDecomp)/include
#
#SRC = $(wildcard *.f90)
#OBJ=$(SRC:.f90=.o)
#
#all : $(OBJ)
#
##regle generique
## Compilation des sources (fortran)
#%.o : %.f90
#	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) -c $<
#
#
#clean:
#	rm -f *.o
