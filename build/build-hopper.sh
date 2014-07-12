#!/bin/bash

module swap PrgEnv-pgi PrgEnv-gnu
module load szip cray-petsc scons

scons --with-mpi --extrapath=/global/homes/a/asarje/local/hdf5-1.8.9-edison-gnu,/global/homes/a/asarje/local/tiff-4.0.2-hopper,/global/homes/a/asarje/local/boost_1_49_0 -j 24
