lAmitex  =../../libAmitex

#amitex_fftp-config : script permettant de connaitre l'emplacement des modules (includes) et des librairies
AMITEX_FFTP_CONFIG = $(lAmitex)/bin/amitex_fftp-config

# Amitex flags
F90FLAGS = `$(AMITEX_FFTP_CONFIG) --fcflags` 

LIBS = `$(AMITEX_FFTP_CONFIG) --plugin-libs`

F90 = `$(AMITEX_FFTP_CONFIG) --mpifc`


OBJ = conventions_2decomp_fftw.o test_DTT.o

#fft_test_r2c: $(OBJ)
#	$(F90) -o $@ $(OBJ1) $(LIBS) $(OMP)

conventions_2decomp_fftw: conventions_2decomp_fftw.o
	$(F90) -o $@ conventions_2decomp_fftw.o $(LIBS) $(F90FLAGS)

test_DTT: test_DTT.o
	$(F90) -o $@ test_DTT.o $(LIBS) $(F90FLAGS)

clean:
	rm -f *.o *.log conventions_2decomp_fftw test_DTT

%.o : %.f90
	$(F90) $(F90FLAGS) -c $<
