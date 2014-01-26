##
# Project: HipGISAXS
#
# File: Makefile
# Created: April 5, 2013
#
# Author: Abhinav Sarje <asarje@lbl.gov>
##

## location of all the required packages
MPI_DIR = /usr/local/openmpi
BOOST_DIR = /usr/local/boost
HDF5_DIR = /usr/local/hdf5-1.8.9
SZ_DIR = /usr/local/szip-2.1
ZLIB_DIR = /usr/local/zlib-1.2.7
TIFF_DIR = /usr/local/tiff-4.0.2
CUDA_DIR = /usr/local/cuda
HIPGISAXS_DIR = /home/asarje/hipgisaxs.git

## the compiler
CXX = $(MPI_DIR)/bin/mpicxx

## name of the needed binary
BINARY = hipgisaxs_fit

## some flags
CXX_FLAGS = -m64 -fopenmp -std=c++0x
CXX_FLAGS += -g

## hipgisaxs include and library linking flags
HIPGISAXS_INCL = -I $(HIPGISAXS_DIR)/include
HIPGISAXS_LDFLAGS = -L $(HIPGISAXS_DIR)/lib
HIPGISAXS_LIBS = -lhipgisaxs

## cuda include and libraries linking flags
CUDA_INCL = -I $(CUDA_DIR)/include
CUDA_LDFLAGS = -L $(CUDA_DIR)/lib64
CUDA_LIBS = -lcudart

## boost include and libraries linking flags
BOOST_INCL = -I $(BOOST_DIR)/include
BOOST_LDFLAGS = -L $(BOOST_DIR)/lib
BOOST_LIBS = -lboost_system -lboost_filesystem -lboost_timer -lboost_chrono

## hdf5 include and library linking flags
HDF5_INCL = -I $(HDF5_DIR)/include
HDF5_LDFLAGS = -L $(HDF5_DIR)/lib -L $(SZ_DIR)/lib -L $(ZLIB_DIR)/lib
HDF5_LIBS = -lhdf5 -lsz -lz

## tiff include and library linking flags
TIFF_INCL = -I $(TIFF_DIR)/include
TIFF_LDFLAGS = -L $(TIFF_DIR)/lib
TIFF_LIBS = -ltiff

## put things together
ALL_INCL = $(HIPGISAXS_INCL) $(CUDA_INCL) $(BOOST_INCL) $(HDF5_INCL) $(TIFF_INCL)
ALL_LDFLAGS = $(HIPGISAXS_LDFLAGS) $(CUDA_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(TIFF_LDFLAGS)
ALL_LIBS = $(HIPGISAXS_LIBS) $(CUDA_LIBS) $(BOOST_LIBS) $(HDF5_LIBS) $(TIFF_LIBS)

## generate the binary
$(BINARY): hipgisaxs_fit.cpp
	$(CXX) -o $@ $^ $(CXX_FLAGS) $(ALL_INCL) $(ALL_LDFLAGS) $(ALL_LIBS)

.PHONY: clean

clean:
	rm -f *.o $(BINARY)
