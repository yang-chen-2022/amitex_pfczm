include ../../src/Makefile.inc

INCLUDE = -I../../include
LIBS = -L../../lib -l2decomp_fft $(LIBFFT)

OBJ = fft_test_r2c_mg.o

fft_test_r2c_mg: $(OBJ)
	$(F90) -o $@ $(OBJ) $(LIBS) $(OMP)

clean:
	rm -f *.o fft_test_r2c_mg

%.o : %.f90
	$(F90) $(INCLUDE) $(OPTIONS) $(F90FLAGS) -c $<
