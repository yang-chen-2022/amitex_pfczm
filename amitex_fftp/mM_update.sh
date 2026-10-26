#!/bin/bash
if [ ! -d "libAmitex" -a ! -d "lib_extern" ]; then
   echo "ERROR : mM_update must be launched in 'amitex_fftp' main directory"
   exit 0
fi

rep0=$(pwd)

#===================================================================== GO TO mM
cd ../mM/dd/trunk/bin

#======================================================================== CLEAN
make cleanall
rm makefile
rm -rf *db*

#======================================================================= UPDATE
svn update

#====================================================================== COMPILE

#-- config (create makefile)
make -f config linux

#-- compile various mM

# standard (non parallel)
make mm

# parallel + AMITEX coupling
make mMAMITEX 	

# parallel + AMITEX coupling + DEBUG
make mMAMITEX VER=db 

# parallel (MPI) 
make mmp 	

# standard and graphical
make gmm	

# tool cam (vizualize film.bin) 
make cam

# tool fil2par (generate vtk files, read in paraview)
make film2para

#===================================================================== GO BACK
cd $rep0

