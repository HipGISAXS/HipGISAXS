#!/bin/bash

#TIFF_DIR=$HOME/local/tiff-4.0.4-babbage
#BOOST_DIR=$HOME/local/boost_1_59_0-babbage
TIFF_DIR=$HOME/babbage/temp_usg.hipgisaxs.git/TIFF

CXX=mpiicpc CC=mpiicc scons --with-mpi --with-mic --extrapath=$BOOST_DIR,$TIFF_DIR -j 16
