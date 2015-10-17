#!/bin/bash

CUDA_TOOLKIT_PATH=$CRAY_CUDATOOLKIT_DIR CXX=CC CC=cc scons --with-mpi --with-cuda --extrapath=$BOOST_DIR -j 16
