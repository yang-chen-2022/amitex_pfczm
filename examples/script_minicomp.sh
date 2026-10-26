#!/bin/bash
#trap 'exit' INT

source ../../amitex_fftp/env_amitex.sh

# EXECUTABLE amitex_fftp
AMITEX="../../amitex_fftp/libAmitex/src/user_pfczm/amitex_fftp"

MPIRUNalias='mpirun -np 10'

MICR="microstructures/minicomp/iUC3/vp0.01/h0.001/iUC3_h0.001_notched.vtk"
MATE="materials/mat_minicomp.xml"
LOAD="loads/load_minicomp.xml"
ALGO="algorithms/algo_Miehe2.xml"
RESU="results/minicomp"

mkdir $RESU
$MPIRUNalias $AMITEX -nm $MICR -m $MATE -c $LOAD -a $ALGO -s $RESU"/test"





