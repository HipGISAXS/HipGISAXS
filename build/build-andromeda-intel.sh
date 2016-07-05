#!/bin/bash

VTUNE_DIR=/opt/intel/vtune_amplifier_xe_2013
BOOST_DIR=/usr/local/boost
TIFF_DIR=/usr/local/tiff-4.0.2
PAPI_DIR=/usr/local/papi-5.1.0

EXTRA_DIRS=$VTUNE_DIR,$BOOST_DIR,$TIFF_DIR,$PAPI_DIR

#CC=icc CXX=icpc scons --with-papi --extrapath=$EXTRA_DIRS -j 32
CC=icc CXX=icpc scons --dbgo --with-papi --extrapath=$EXTRA_DIRS -j 32
