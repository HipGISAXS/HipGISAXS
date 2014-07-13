#!/bin/bash

#module load scons
#module load gcc/4.8.2
#module load boost szip cray-petsc

scons --with-mpi --extrapath=/global/homes/a/asarje/local/hdf5-1.8.9-edison-intel,/global/homes/a/asarje/local/tiff-4.0.2-edison-intel,/global/homes/a/asarje/local/boost_1_49_0,/usr/common/usg/szip/2.1/intel -j 24
