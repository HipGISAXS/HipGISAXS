#!/bin/bash

#scons --extrapath=/usr/local/openmpi,/usr/local/cuda,/usr/local/boost,/usr/local/hdf5-1.8.9-intel,/usr/local/tiff-4.0.2,/usr/local/szip-2.1 --with-mpi --with-cuda -j 32
#CC=icc CXX=icpc scons --dbgo --extrapath=/usr/local/boost,/usr/local/tiff-4.0.2 -j 32
CC=icc CXX=icpc scons --with-papi --extrapath=/usr/local/boost,/usr/local/tiff-4.0.2,/usr/local/papi-5.1.0 -j 32
