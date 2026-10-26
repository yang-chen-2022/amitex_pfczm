#!/bin/bash


mkdir mate     > /dev/null 2>&1
mkdir algo     > /dev/null 2>&1
mkdir load     > /dev/null 2>&1
mkdir micro    > /dev/null 2>&1
mkdir behavior > /dev/null 2>&1
mkdir results  > /dev/null 2>&1 
mkdir script   > /dev/null 2>&1 


if [[ $1 == -help ]];then
  echo "USAGE : amitex_case -clean remove mate, algo, load, micro, behavior, results, script"
  echo "                     -default (matincl_elasiso)"
  echo "                     -matincl_elasiso"
  echo "                     -polyx_elasaniso"
  exit 0

elif [[ $1 == -clean ]]; then
  rm -rf mate
  rm -rf algo  
  rm -rf load
  rm -rf micro
  rm -rf behavior
  rm -rf results  
  rm -rf script 
  
elif [[ $1 == -default ]]; then
  echo "AMITEX_PATH="$AMITEX_PATH
  cp $AMITEX_PATH/cas_tests/microstructures/al2p5_64/mat_al2p5_64.vtk micro
  cp $AMITEX_PATH/cas_tests/algorithmes/algo_default.xml algo
  cp $AMITEX_PATH/cas_tests/mate/mat_incl_pore1.xml mate
  cp $AMITEX_PATH/cas_tests/chargements/char_traction_1.xml load
  cp $AMITEX_PATH/cas_tests/scripts/script_matincl_elasiso.sh script

elif [[ $1 == -matincl_elasiso ]]; then
  echo "AMITEX_PATH="$AMITEX_PATH
  cp $AMITEX_PATH/cas_tests/microstructures/al2p5_64/mat_al2p5_64.vtk micro
  cp $AMITEX_PATH/cas_tests/algorithmes/algo_default.xml algo
  cp $AMITEX_PATH/cas_tests/mate/mat_incl_pore1.xml mate
  cp $AMITEX_PATH/cas_tests/chargements/char_traction_1.xml load
  cp $AMITEX_PATH/cas_tests/scripts/script_matincl_elasiso.sh script

elif [[ $1 == -polyx_elasaniso ]]; then
  echo "AMITEX_PATH="$AMITEX_PATH
  cp $AMITEX_PATH/cas_tests/microstructures/27G_65/zone_27G_65.vtk micro
  cp $AMITEX_PATH/cas_tests/mate/mat_polyx_elasaniso.xml mate
  cp -r $AMITEX_PATH/cas_tests/mate/vec_orient_ISOfdo_21^3 mate
  cp $AMITEX_PATH/cas_tests/algorithmes/algo_default.xml algo
  cp $AMITEX_PATH/cas_tests/chargements/char_traction_1.xml load
  cp $AMITEX_PATH/cas_tests/scripts/script_polyx_elasaniso.sh script

else
  echo "ERROR : no or no valid argument use -help"
  exit 0
fi



