#!/bin/bash
#trap 'exit' INT

source ../../amitex_fftp/env_amitex.sh

# EXECUTABLE amitex_fftp
AMITEX="../../amitex_fftp/libAmitex/src/user_pfczm/amitex_fftp"

MPIRUNalias='mpirun -np 10'

MICR="microstructures/periodic2D/iUC1/vp0.1/h0.001/iUC1_h0.001.vtk"
MATE="materials/mat_periodic2D.xml"
LOAD="loads/load_periodic2D.xml"
ALGO="algorithms/algo_Miehe2.xml"
RESU="results/periodic2D"

mkdir $RESU
$MPIRUNalias $AMITEX -nm $MICR -m $MATE -c $LOAD -a $ALGO -s $RESU"/test"





