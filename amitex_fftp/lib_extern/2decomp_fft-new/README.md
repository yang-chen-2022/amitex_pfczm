2decomp "new" can be dowloaded at 
	
	https://github.com/2decomp-fft/2decomp-fft

documentation at :

	https://github.com/2decomp-fft/2decomp-fft.github.io/tree/main
========================================================================

Below a few information to compile and test the library


INSTALLATION + TESTS
=====================

Examples for building:
https://github.com/2decomp-fft/2decomp-fft/blob/dev/.github/workflows/Build_all.yml
https://github.com/2decomp-fft/2decomp-fft/blob/dev/.github/workflows/RunnerFullCI.yml#L35

Build on AMITEX
---------------

source env_amitex_gcc810-openmpi310.sh
export FC=mpif90
rm -rf ./build
#cmake -S . -B build -DCMAKE_BUILD_TYPE=Dev -DFFT_Choice=fftw_f03 -DBUILD_TESTING=ON
#cmake -S . -B build -DFFT_Choice=fftw_f03 -DBUILD_TESTING=ON -DMPIEXEC_MAX_NUMPROCS=14
#cmake -S . -B build -DBUILD_TESTING=ON 
cmake -S . -B build -DFFT_Choice=fftw_f03 -DBUILD_TESTING=ON -DNX=256 -DNY=256 -DNZ=256 -DMPIEXEC_MAX_NUMPROCS=14
#ccmake build
cmake --build build --verbose
cmake --install build
ctest --test-dir build

ATTENTION :  - DMPIEXEC_MAX_NUMPROCS=14 modifies MAX_NUMPROCS... bug reported
               must be modified with ccmake build ([toggle][configure][generate]

LANCEMENT DES TESTS
====================

ctest --test-dir build
-> log files in build/RunTests 
-> standard output  in build/Testing/Temporary/CTestCostData.txt
-> by default tests are run on 1 proc
-> modify MPIEXEC_MAX_NUMPROCS using ccmake build ([toggle][configure][generate]

TESTING without ctest (independant of MPIEXEC_MAX_NUMPROCS)
-------------------------------------------------------------------

depuis (2decomp)/.
mpirun -n 14 build/bin/fft_c2c_x 7 2 256 256 256 20  # p_row, pcol, nx,ny,nz,ntest
mpirun -n 14 build/bin/io_var_test 7 2 256 256 256   # p_row, pcol, nx,ny,nz : ce test est long avec 256, utilise aussi des 512^3
