#!/bin/bash

#TIFFDIR=$HOME/local/tiff-4.0.2-edison-intel
#TIFFDIR=$HOME/local/tiff-4.0.4-cori
#TIFFDIR=/global/homes/a/asarje/local/tiff-4.0.2-edison-intel
TIFFDIR=/project/projectdirs/als/local/tiff-4.0.6
PAPIDIR=/opt/cray/papi/5.4.1.3
VTUNEDIR=${VTUNE_AMPLIFIER_XE_2016_DIR}
BOOSTDIR=${BOOST_ROOT}

ALLDIRS=$TIFFDIR,$PAPIDIR,$VTUNEDIR,$BOOSTDIR

#CC=cc CXX=CC scons --with-mpi --with-papi --extrapath=$ALLDIRS -j12
CC=cc CXX=CC scons --dbgo --with-mpi --with-papi --extrapath=$ALLDIRS -j12
