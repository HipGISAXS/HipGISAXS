#!/bin/bash

TIFFDIR=$HOME/local/tiff-4.0.2-edison-intel
#TIFFDIR=$HOME/local/tiff-4.0.4-cori
PAPIDIR=/opt/cray/papi/5.4.1.3
VTUNEDIR=${VTUNE_AMPLIFIER_XE_2016_DIR}
BOOSTDIR=${BOOST_ROOT}

ALLDIRS=$TIFFDIR,$PAPIDIR,$VTUNEDIR,$BOOSTDIR

#echo "CC=cc CXX=CC scons --with-mpi --with-papi --extrapath=$ALLDIRS -j12"
#CC=cc CXX=CC scons --with-mpi --with-papi --extrapath=$ALLDIRS -j12
echo "CC=cc CXX=CC scons --dbgo --with-mpi --with-papi --extrapath=$ALLDIRS -j12"
CC=cc CXX=CC scons --dbgo --with-mpi --with-papi --extrapath=$ALLDIRS -j12
