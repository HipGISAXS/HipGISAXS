#!/bin/bash

#scons --extrapath=/usr/local/openmpi,/usr/local/cuda,/usr/local/boost,/usr/local/hdf5-1.8.9,/usr/local/tiff-4.0.2,/usr/local/szip-2.1 --with-mpi --with-cuda --detail-time --detail-mem --verbose -j 32
#scons --extrapath=/usr/local/openmpi,/usr/local/cuda,/usr/local/boost,/usr/local/hdf5-1.8.9,/usr/local/tiff-4.0.2,/usr/local/szip-2.1 --with-mpi --with-cuda --detail-time -j 32
scons --extrapath=/usr/local/openmpi,/usr/local/cuda,/usr/local/boost,/usr/local/tiff-4.0.2,/usr/local/szip-2.1 --with-mpi --with-cuda -j 32
#scons --extrapath=/usr/local/openmpi,/usr/local/cuda,/usr/local/boost,/usr/local/hdf5-1.8.9,/usr/local/tiff-4.0.2,/usr/local/szip-2.1 --with-mpi --with-cuda --with-parallel-hdf5 -j 32
