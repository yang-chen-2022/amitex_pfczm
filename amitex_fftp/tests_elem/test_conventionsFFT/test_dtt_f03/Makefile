#Debug gcc
F90=gfortran
DEBUG=-Wall -g -fcheck=all -fbounds-check -ffpe-trap=zero,overflow

F90=ifort
DEBUG=-g -check all -fpe0 -warn all,noexternal -traceback -debug extended
DEBUG=
F90FLAGS=$(DEBUG) 


#fftw access
INCLUDE = -I/usr/include
LIBS = -L/usr/lib64 -lfftw3f -lfftw3

OBJ = test_DTT.o

test_DTT: test_DTT.o
	$(F90) -o $@ test_DTT.o $(LIBS) $(INCLUDE) $(DEBUG) 
	./test_DTT

clean:
	rm -f *.o test_DTT

%.o : %.f90
	$(F90) $(INCLUDE) $(LIBS) $(F90FLAGS) $(DEBUG) -c $<
