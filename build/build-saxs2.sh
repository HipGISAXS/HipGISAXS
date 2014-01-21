#!/bin/bash

scons --with-mpi --with-cuda --extrapath=/usr/local/openmpi-1.6.1,/usr/local/boost,/usr/local/cuda,/usr/local/hdf5-1.8.9,/usr/local/szip-2.1 -j 12
