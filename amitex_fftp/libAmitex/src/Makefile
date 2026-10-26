#=======================================================================
#
#  COMPILE :
#
#   * amitex_fftp : "executable" file in (amitex)/libAmitex/bin
#
#   * libamitex.a(.so) : static (dynamic) library in (amitex)/libAmitex/lib
#                        include everythig but the main program 'amitex_fftp'
#
#   * libamitex_plugin.so : dynamic library in (amitex)/libAmitex/lib
#                           main program 'amitex_fftp'
#
#=======================================================================
#
# amitex_fftp-config tool :
# -------------------------
# To access help, once AMITEX installed :
#    - source (amitex)/env_amitex.sh
#    - amitex_fftp-config --help
#
#=======================================================================

lAmitex=..
AMITEX_CFG=$(lAmitex)/bin/amitex_fftp-config

# Compilateur
F90=$(shell $(AMITEX_CFG) --mpifc)
F90FLAGS := $(shell $(AMITEX_CFG) --fcflags --fcflags-opt) 

CC=$(shell $(AMITEX_CFG) --cc)
CFLAGS := $(shell $(AMITEX_CFG) --fcflags-opt) 

#nom fichier principal
MAIN = amitex_fftp

#fichiers objets de modules (+read_geom)
MOD  = simu_mod.o amitex_mod.o resolution_mod.o NL_base_mod.o sortie_std_mod.o error_mod.o material_mod.o param_algo_mod.o loading_mod.o io_amitex_mod.o io2_amitex_mod.o green_mod.o algo_functions_mod.o linear_mod.o read_geom.o field_mod.o non_local_mod.o amitex_lib.o amitex_mpi_mod.o

#fichiers objets de modules "user"
MOD_user = amitex_user_mod.o resolution_user_mod.o standard_user_mod.o non_local_user_mod.o

#fichiers objets des fonction C utiles pour l'appel a des umat externes
UMAT_C = umat_load.o umat_call.o umatD_call.o

#acces aux interfaces des modules
INCLUDE = -I../include

# for bib partagée
OPTIONS += -fPIC
CFLAGS += -fPIC

#fichiers sources et objets des comportements materiaux
LOI_src = $(wildcard ./materiaux/*.f90)
LOI = $(LOI_src:.f90=.o)

# --libs-deps-only : to link with with libraries FoX and 2decomp
LIBS=$(shell $(AMITEX_CFG) --libs-deps-only)

LIBSamitex = -luser_functions

#=============================================================
#                                     CREATION DE L'EXECUTABLE 
#                                      (ET DES LOIS MATERIAUX)
#                                         PUIS DE LA LIBRAIRIE

all: $(MAIN) LOI lib bin includes

# Creation de l'executable a partir des fichiers .o et .mod dans src  
$(MAIN): $(MOD) $(MOD_user) $(UMAT_C) $(MAIN).o LIBUSER_FUNCTIONS
	$(F90) -o $@ *.o $(LIBS) $(F90FLAGS) $(LIBSamitex) -ldl

# Creation de la librairie dynamique user_functions (module utilisateur d'interruption du programme)
# et suppression du .o avant creation executable et librairies libamitex 
LIBUSER_FUNCTIONS: user_functions_mod.o
	$(F90) -o libuser_functions.so -shared user_functions_mod.o
	rm user_functions_mod.o
	mv libuser_functions.so ../lib/

# Move executable in ../bin
bin: $(MAIN) 
	mv $(MAIN) ../bin/

# Create and move libamitex.(a|so) in ../lib
lib: $(MAIN) 
	ar qc libamitex.a $(MOD) $(MOD_user) $(UMAT_C)
	mv libamitex.a ../lib/
	$(F90) -shared -o libamitex.so  $(MOD) $(MOD_user) $(UMAT_C) $(LIBS)  $(LIBSamitex) -ldl
	mv libamitex.so ../lib/
	$(F90) -shared -o libamitex_plugin.so $(MAIN).o -lamitex $(LIBS)  $(LIBSamitex) -ldl
	mv libamitex_plugin.so ../lib/

# Copy .mod files in ../include  #(+ remove ./include/*_user_mod.mod)
includes: $(MAIN)
	cp *.mod ../include/
	#rm ../include/*_user_mod.mod

#=============================================================
#                                             OPTIONS DE DEBUG
#debug intel
#DEBUG=-g -check all -fpe0 -warn -traceback -debug extended
#To remove "warning #8889: Explicit interface or EXTERNAL declaration is required"
#   for all MPI calls and calls to read_geom, read_geom_composite/continuous 
#DEBUG=-g -check all -fpe0 -warn all,noexternal -traceback -debug extended
#ci-dessous propose sur le web pour debuggue
#DEBUG=-O0 -fpe-all=0 -traceback -ftrapuv -fp-stack-check -check all -debug
#DEBUG=-g -check all
#DEBUG=-warn all
#DEBUG=-ftrapuv  # Initialize stack local variables to an unusual value to aid error detection.

#debug gcc
#DEBUG=-Wall
#ci-dessous propose sur le web pour debuggue
#DEBUG=-Wall -fbounds-check -ffpe-trap=invalid,zero,overflow
# supprimer invalid tant qu'on utilise des NaN dans le code!
#DEBUG=-Wall -fbounds-check -ffpe-trap=zero,overflow
#DEBUG=-g -fcheck=all
# un mix :
#DEBUG=-Wall -g -fcheck=all -fbounds-check -ffpe-trap=zero,overflow
DEBUG=


#=============================================================
#                                                 COMPILATIONS

# Compilation du programme principal
$(MAIN).o: $(MOD) $(MOD_user) $(MAIN).f90
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $(MAIN).f90 

# compilation des modeles (dans le dossier materiaux)
LOI :
	cd ./materiaux && $(MAKE) 
	cd ./materiauxK && $(MAKE) 


# Compilation des modules avec dependances
user_functions_mod.o : user_functions_mod.f90
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -fPIC -c $<

error_mod.o : error_mod.f90 amitex_mpi_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

field_mod.o : field_mod.f90 amitex_mpi_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

linear_mod.o : linear_mod.f90 error_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

io_amitex_mod.o: io_amitex_mod.f90  error_mod.o field_mod.o amitex_mpi_mod.o
	$(F90) $(INCLUDE) $(modFox) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

amitex_mod.o: amitex_mod.f90 error_mod.o io_amitex_mod.o amitex_mpi_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

param_algo_mod.o: param_algo_mod.f90 error_mod.o io_amitex_mod.o amitex_mpi_mod.o
	$(F90) $(INCLUDE) $(modFox) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

material_mod.o: material_mod.f90 io_amitex_mod.o error_mod.o amitex_mod.o param_algo_mod.o linear_mod.o algo_functions_mod.o amitex_mpi_mod.o
	$(F90) $(INCLUDE) $(modFox) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

loading_mod.o: loading_mod.f90 material_mod.o error_mod.o param_algo_mod.o io_amitex_mod.o amitex_mpi_mod.o
	$(F90) $(INCLUDE) $(modFox) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

green_mod.o: green_mod.f90 error_mod.o param_algo_mod.o field_mod.o amitex_mod.o amitex_mpi_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

algo_functions_mod.o: algo_functions_mod.f90  linear_mod.o param_algo_mod.o error_mod.o amitex_mpi_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

io2_amitex_mod.o: io2_amitex_mod.f90  material_mod.o loading_mod.o param_algo_mod.o field_mod.o io_amitex_mod.o error_mod.o amitex_mpi_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

sortie_std_mod.o: sortie_std_mod.f90 loading_mod.o material_mod.o param_algo_mod.o amitex_mod.o error_mod.o field_mod.o amitex_mpi_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

non_local_mod.o: non_local_mod.f90 param_algo_mod.o error_mod.o field_mod.o green_mod.o amitex_mod.o material_mod.o non_local_user_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

NL_base_mod.o: NL_base_mod.f90 sortie_std_mod.o loading_mod.o material_mod.o io2_amitex_mod.o algo_functions_mod.o param_algo_mod.o green_mod.o amitex_mod.o error_mod.o field_mod.o non_local_mod.o user_functions_mod.o amitex_mpi_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<
	
resolution_mod.o: resolution_mod.f90 NL_base_mod.o sortie_std_mod.o loading_mod.o material_mod.o io2_amitex_mod.o algo_functions_mod.o param_algo_mod.o green_mod.o amitex_mod.o error_mod.o field_mod.o non_local_mod.o standard_user_mod.o user_functions_mod.o simu_mod.o amitex_mpi_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

amitex_lib.o: amitex_lib.f90 amitex_mod.o green_mod.o field_mod.o error_mod.o simu_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<  
	
simu_mod.o: simu_mod.f90 sortie_std_mod.o loading_mod.o material_mod.o param_algo_mod.o green_mod.o amitex_mod.o field_mod.o io_amitex_mod.o io2_amitex_mod.o non_local_mod.o error_mod.o NL_base_mod.o amitex_mpi_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<  

amitex_mpi_mod.o: amitex_mpi_mod.f90
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<  

# Minimalist user modules

amitex_user_mod.o: amitex_user_mod.f90 amitex_mod.o error_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

resolution_user_mod.o: resolution_user_mod.f90 amitex_user_mod.o amitex_mod.o error_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

standard_user_mod.o: standard_user_mod.f90 amitex_user_mod.o amitex_mod.o error_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

non_local_user_mod.o: non_local_user_mod.f90 amitex_mod.o error_mod.o
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<


# compilation des fichiers d'interface umat
umat_load.o : umat_load.c
	$(CC) $(CFLAGS) -c $<

umat_call.o : umat_call.c
	$(CC) $(CFLAGS) umat_call.c -o umat_call.o -c

umatD_call.o : umatD_call.c
	$(CC) $(CFLAGS) umatD_call.c -o umatD_call.o -c

# Regle generique
# Compilation des sources (fortran) autres que le programme principal
# -> error_mod, linear_mod  
%.o : %.f90
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c $<

# read_geom : 	compilé (si f90 existe) ou copié du répertoire read_geom 
#		choix du .o cohérent avec le compilateur utilisé
G_TEST := $(shell "$(F90)" --version | grep -ic "gnu")
I_TEST := $(shell "$(F90)" --version | grep -ic "intel")

ifneq ($(wildcard read_geom.f90),)
read_geom.o : read_geom.f90
#	$(F90) $(INCLUDE) $(OPTIONS) -O0 $(DEBUG) -c read_geom.f90
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) $(DEBUG) -c read_geom.f90
else
read_geom.o : error_mod.o amitex_mpi_mod.o
ifneq ($(I_TEST),0)
	cp read_geom/read_geom.o_intel15 ./read_geom.o
endif
ifneq ($(G_TEST),0)
	cp read_geom/read_geom.o_gnu482 ./read_geom.o
endif
endif

# supprime les .o .a et .mod
clean:
	rm -rf *.o *.mod ../lib/*.a ../include/*.mod ../include/*.smod ../bin/amitex_fftp ../lib/*.so ../lib/decomp2d* > /dev/null 2>&1
	rm -f *__genmod.f90  > /dev/null 2>&1
	cd ./materiaux && $(MAKE) clean 
	cd ./materiauxK && $(MAKE) clean 

