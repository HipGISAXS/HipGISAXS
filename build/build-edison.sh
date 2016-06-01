#!/bin/bash

#TIFFDIR=/global/homes/a/asarje/local/tiff-4.0.2-edison-intel
TIFFDIR=/project/projectdirs/als/local/tiff-4.0.6

CC=cc CXX=CC scons --with-mpi --extrapath=$BOOST_ROOT,$TIFFDIR -j12
