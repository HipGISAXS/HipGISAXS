#!/bin/bash

TIFFDIR=/global/homes/a/asarje/local/tiff-4.0.2-edison-intel

## for TAU profiling
#export TAU_MAKEFILE=${HOME}/local/tau-2.24-edison/craycnl/lib/Makefile.tau-intel-mpi
#export TAU_OPTIONS="-optVerbose -optCompInst"
#export PATH=/global/homes/a/asarje/local/tau-2.24-edison/craycnl/bin:$PATH
#export LD_LIBRARY_PATH=/global/homes/a/asarje/local/tau-2.24-edison/craycnl/lib:$LD_LIBRARY_PATH
#export CRAY_LD_LIBRARY_PATH=/global/homes/a/asarje/local/tau-2.24-edison/craycnl/lib:$CRAY_LD_LIBRARY_PATH

CC=cc CXX=CC scons --with-mpi --extrapath=$BOOST_ROOT,$TIFFDIR -j12
#scons --with-mpi --with-papi --extrapath=$BOOST_ROOT,$TIFFDIR -j16
#scons --with-mpi --with-tau --extrapath=$BOOST_ROOT,$TIFFDIR -j12
