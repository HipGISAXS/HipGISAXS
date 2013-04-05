##
# Project: HipGISAXS
#
# File: Makefile
# Created: April 5, 2013
#
# Author: Abhinav Sarje <asarje@lbl.gov>
##

## location of all the required packages
MPI_DIR = /home/asarje/opt/openmpi-1.6.4-intel
BOOST_DIR = /home/asarje/opt/boost_1_53_0-intel
HDF5_DIR = /home/asarje/opt/hdf5-1.8.10-patch1-parallel-intel
SZ_DIR = /home/asarje/opt/szip-2.1-intel
TIFF_DIR = /home/asarje/opt/tiff-4.0.2-intel
HIPGISAXS_DIR = /home/asarje/hipgisaxs.git

## the compiler
CXX = $(MPI_DIR)/bin/mpicxx

## name of the needed binary
BINARY = test_hipgisaxs

## some flags
CXX_FLAGS = -m64 -openmp -std=c++0x

## hipgisaxs include and library linking flags
HIPGISAXS_INCL = -I $(HIPGISAXS_DIR)/src
HIPGISAXS_LDFLAGS = -L $(HIPGISAXS_DIR)/lib
HIPGISAXS_LIBS = -lhipgisaxs

## boost include and libraries linking flags
BOOST_INCL = -I $(BOOST_DIR)/include
BOOST_LDFLAGS = -L $(BOOST_DIR)/lib
BOOST_LIBS = -lboost_system -lboost_filesystem -lboost_timer -lboost_chrono

## hdf5 include and library linking flags
HDF5_INCL = -I $(HDF5_DIR)/include
HDF5_LDFLAGS = -L $(HDF5_DIR)/lib -L $(SZ_DIR)/lib
HDF5_LIBS = -lhdf5 -lsz

## tiff include and library linking flags
TIFF_INCL = -I $(TIFF_DIR)/include
TIFF_LDFLAGS = -L $(TIFF_DIR)/lib
TIFF_LIBS = -ltiff

## put things together
ALL_INCL = $(HIPGISAXS_INCL) $(BOOST_INCL) $(HDF5_INCL) $(TIFF_INCL)
ALL_LDFLAGS = $(HIPGISAXS_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(TIFF_LDFLAGS)
ALL_LIBS = $(HIPGISAXS_LIBS) $(BOOST_LIBS) $(HDF5_LIBS) $(TIFF_LIBS)

## generate the binary
$(BINARY): test_hipgisaxs.cpp
	$(CXX) -o $@ $^ $(CXX_FLAGS) $(ALL_INCL) $(ALL_LDFLAGS) $(ALL_LIBS)

.PHONY: clean

clean:
	rm -f *.o $(BINARY)
