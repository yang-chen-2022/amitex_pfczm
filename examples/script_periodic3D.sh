#!/bin/bash
#trap 'exit' INT

source ../../amitex_fftp/env_amitex.sh

# EXECUTABLE amitex_fftp
AMITEX="../../amitex_fftp/libAmitex/src/user_pfczm/amitex_fftp"

MPIRUNalias='mpirun -np 10'

MICR="microstructures/periodic3D/iUC3/vp0.01/h0.001/iUC3_h0.001_notched.vtk"
MATE="materials/mat_periodic3D.xml"
LOAD="loads/load_periodic3D.xml"
ALGO="algorithms/algo_Miehe2.xml"
RESU="results/periodic3D"

mkdir $RESU
$MPIRUNalias $AMITEX -nm $MICR -m $MATE -c $LOAD -a $ALGO -s $RESU"/test"





