#!/bin/bash

module unload darshan
module load gcc				## for c++11 compatibility
module load scons
module load cray-hdf5-parallel
module load boost
#module load gcc/4.8.2
#module load boost szip cray-petsc

TIFFDIR=/global/homes/a/asarje/local/tiff-4.0.2-edison-intel


## for TAU profiling
#export TAU_MAKEFILE=${HOME}/local/tau-2.24-edison/craycnl/lib/Makefile.tau-intel-mpi
#export TAU_OPTIONS="-optVerbose -optCompInst"
#export PATH=/global/homes/a/asarje/local/tau-2.24-edison/craycnl/bin:$PATH
#export LD_LIBRARY_PATH=/global/homes/a/asarje/local/tau-2.24-edison/craycnl/lib:$LD_LIBRARY_PATH
#export CRAY_LD_LIBRARY_PATH=/global/homes/a/asarje/local/tau-2.24-edison/craycnl/lib:$CRAY_LD_LIBRARY_PATH


#scons --with-mpi --extrapath=/global/homes/a/asarje/local/hdf5-1.8.9-edison-intel,/global/homes/a/asarje/local/tiff-4.0.2-edison-intel,/global/homes/a/asarje/local/boost_1_49_0,/usr/common/usg/szip/2.1/intel -j 24

scons --with-mpi --extrapath=$BOOST_ROOT,$TIFFDIR -j12

#scons --with-mpi --with-tau --extrapath=$BOOST_ROOT,$TIFFDIR -j12
