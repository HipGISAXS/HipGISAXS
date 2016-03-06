## the following are currently defined platforms:
## * andromeda-gnu
## * andromeda-intel
## * nersc-gnu
## * nersc-intel
## * titan-gnu
## * titan-intel
## Note: all the Intel versions need to be tested first, so dont use them for now.

default:
	( $(MAKE) error )

gnu:
	(	\
		)

intel:
	(	\
		)

bragg-gnu:
	( $(MAKE) hipgisaxs install                               \
		"CXX = mpicxx"                                          \
		"CXX_FLAGS = -std=c++11 -fopenmp -DUSE_MPI"             \
		"CXX_FLAGS += -O3"                                      \
		"OPT_FLAGS = -O3"                                       \
		"DBG_FLAGS = -g"                                        \
		"MPI_LIBS = -lmpi -lmpi_cxx"                            \
		"PETSC_INCL = -I $(PETSC_DIR)/include"                  \
		"PETSC_LDFLAGS = -L $(PETSC_DIR)/lib"                   \
		"PETSC_LIBS = -lpetsc"                                  \
		"HDF5_LIBS = -lhdf5 -lz"                                \
		)

andromeda-gnu:
	( $(MAKE) hipgisaxs install                               \
		"MPI_DIR = /usr/local/openmpi"                          \
		"PETSC_DIR = /usr/local/petsc-3.6.3"                    \
		"PETSC_VERSION = \"3.6\""                               \
		"HIPGISAXS_DIR = ../.."                                 \
		"HDF5_DIR = /usr/local/hdf5-1.8.9"                      \
		"TIFF_DIR = /usr/local/tiff-4.0.2"                      \
		"BOOST_DIR = /usr/local/boost"                          \
		"SZIP_DIR = /usr/local/szip-2.1"                        \
		"CUDA_DIR = /usr/local/cuda"                            \
		"CXX = mpicxx"                                          \
		"CXX_FLAGS = -std=c++11 -fopenmp -DUSE_MPI"             \
		"CXX_FLAGS += -O3"                                      \
		"OPT_FLAGS = -O3"                                       \
		"DBG_FLAGS = -g"                                        \
		"MPI_LIBS = -lmpi -lmpi_cxx"                            \
		"PETSC_LIBS = -lpetsc"                                  \
		"HDF5_LIBS = -lhdf5 -lz"                                \
		)

# with intel don't use mpi for now as it has not been tested
andromeda-intel:
	( $(MAKE) hipgisaxs install                               \
		"PETSC_DIR = /usr/local/petsc-3.6.2"                    \
		"PETSC_VERSION = \"3.6\""                               \
		"HIPGISAXS_DIR = $(HOME)/hipgisaxs.edison.git"          \
		"HDF5_DIR = /usr/local/hdf5-1.8.9-intel"                \
		"TIFF_DIR = /usr/local/tiff-4.0.2-intel"                \
		"BOOST_DIR = /usr/local/boost"                          \
		"SZIP_DIR = /usr/local/szip-2.1"                        \
		"CUDA_DIR = /usr/local/cuda"                            \
		"CXX = icpc"                                            \
		"CXX_FLAGS = -std=c++11 -openmp"                        \
		"CXX_FLAGS += -O3"                                      \
		"OPT_FLAGS = -O3"                                       \
		"DBG_FLAGS = -g"                                        \
		"PETSC_LIBS = -lpetsc"                                  \
		"HDF5_LIBS = -lhdf5 -lz"                                \
		)

nersc-intel:
	( $(MAKE) hipgisaxs install                               \
		"PETSC_DIR = $(PETSC_DIR)"                              \
		"PETSC_VERSION = \"3.6\""                               \
		"HIPGISAXS_DIR = ../.."                                 \
		"BOOST_DIR = $(BOOST_DIR)"                              \
		"HDF5_DIR = $(HOME)/local/hdf5-1.8.9-edison-intel"      \
		"TIFF_DIR = $(HOME)/local/tiff-4.0.2-edison-intel"      \
		"TAO_DIR = $(HOME)/local/tao-2.2.0-intel/lib64/python2.6/site-packages/tao"	\
		"CXX = CC"                                              \
		"CXX_FLAGS = -std=c++11 -openmp -DUSE_MPI"              \
		"CXX_FLAGS += -O3"                                      \
		"OPT_FLAGS = -O3"                                       \
		"DBG_FLAGS = -g"                                        \
		"PETSC_LIBS = -lcraypetsc_intel_real"                   \
		"HIPGISAXS_LIBS = -lhipgisaxs -lstdc++"                 \
		"HDF5_LIBS = -lhdf5 -lz"                                \
		)

nersc-gnu:
	( $(MAKE) hipgisaxs install                               \
		"PETSC_DIR = $(PETSC_DIR)"                              \
		"PETSC_VERSION = \"3.6\""                               \
		"HIPGISAXS_DIR = ../.."                                 \
		"BOOST_DIR = $(BOOST_DIR)"                              \
		"HDF5_DIR = $(HOME)/local/hdf5-1.8.9-edison-gnu"        \
		"TIFF_DIR = $(HOME)/local/tiff-4.0.2-edison-gnu"        \
		"TAO_DIR = $(HOME)/local/tao-2.2.0-gnu/lib64/python2.6/site-packages/tao"	\
		"CXX = CC"                                              \
		"CXX_FLAGS = -std=c++11 -fopenmp -DUSE_MPI"             \
		"CXX_FLAGS += -O3"                                      \
		"OPT_FLAGS = -O3"                                       \
		"DBG_FLAGS = -g"                                        \
		"PETSC_LIBS = -lcraypetsc_gnu_real"                     \
		"HIPGISAXS_LIBS = -lhipgisaxs"                          \
		"HDF5_LIBS = -lhdf5 -lz -lsz"                           \
		)

titan-gnu:
	( $(MAKE) hipgisaxs install                               \
		"PETSC_DIR = $(PETSC_DIR)"                              \
		"PETSC_VERSION = \"3.6\""                               \
		"HIPGISAXS_DIR = ../.."                                 \
		"BOOST_DIR = $(BOOST_DIR)"                              \
		"CUDA_DIR = $(CRAY_CUDATOOLKIT_DIR)"                    \
		"CXX = CC"                                              \
		"CXX_FLAGS = -std=c++11 -fopenmp -DUSE_MPI"             \
		"CXX_FLAGS += -O3"                                      \
		"OPT_FLAGS = -O3"                                       \
		"DBG_FLAGS = -g"                                        \
		"PETSC_LIBS = -lcraypetsc_gnu_real"                     \
		"HIPGISAXS_LIBS = -lhipgisaxs -lcudart"                 \
		)

ifdef MPI_DIR
MPI_LDFLAGS = -L $(MPI_DIR)/lib
endif
ifdef PETSC_DIR
PETSC_INCL = -I $(PETSC_DIR)/include
PETSC_LDFLAGS = -L $(PETSC_DIR)/lib
endif
ifdef HDF5_DIR
HDF5_INCL = -I $(HDF5_DIR)/include
HDF5_LDFLAGS = -L $(HDF5_DIR)/lib
ifdef SZIP_DIR
HDF5_LDFLAGS += -L $(SZIP_DIR)/lib
endif
endif
ifdef TIFF_DIR
TIFF_INCL = -I $(TIFF_DIR)/include
TIFF_LDFLAGS = -L $(TIFF_DIR)/lib
endif

PREFIX  = $(PWD)
SRC_DIR = $(PREFIX)
BIN_DIR = $(PREFIX)
PREC_FLAG = -DDOUBLEP    ## double precision is default

ALL_FLAGS = $(CXX_FLAGS) $(PREC_FLAG) $(PETSC_FLAG)

all: andromeda-gnu
	@echo $(PETSC_VERSION)


ifeq ($(PETSC_VERSION), "3.6")
PETSC_FLAG += -DPETSC_36
TAO_DIR =
TAO_INCL =
TAO_LDFLAGS =
TAO_LIBS =
else
ifeq ($(PETSC_VERSION), "3.5")
PETSC_FLAG += -DPETSC_35
TAO_DIR =
TAO_INCL =
TAO_LDFLAGS =
TAO_LIBS =
else
ifdef TAO_DIR
TAO_INCL = -I $(TAO_DIR)/include
TAO_LDFLAGS = -L $(TAO_DIR)/lib
else
TAO_DIR = /usr/local/tao-2.2.0
TAO_INCL = -I $(TAO_DIR)/include
TAO_LDFLAGS = -L $(TAO_DIR)/lib
endif
TAO_LIBS = -ltao
endif
endif

ifdef HIPGISAXS_DIR
HIPGISAXS_INCL = -I $(HIPGISAXS_DIR)/include
HIPGISAXS_LDFLAGS = -L $(HIPGISAXS_DIR)/lib
ifdef CUDA_DIR
HIPGISAXS_INCL += -I $(CUDA_DIR)/include
HIPGISAXS_LDFLAGS += -L $(CUDA_DIR)/lib64
endif
endif
HIPGISAXS_LIBS = -lhipgisaxs -lcudart

TIFF_LIBS = -ltiff

ifdef BOOST_DIR
BOOST_INCL = -I $(BOOST_DIR)/include
BOOST_LDFLAGS = -L $(BOOST_DIR)/lib
endif
BOOST_LIBS = -lboost_system -lboost_filesystem -lboost_timer -lboost_chrono


## all objs and libs                                                                                                                           
ALL_INCL= $(TAO_INCL) $(PETSC_INCL) $(HDF5_INCL) $(TIFF_INCL) $(BOOST_INCL) $(HIPGISAXS_INCL)
ALL_LDFLAGS = $(MPI_LDFLAGS) $(TAO_LDFLAGS) $(GCC_LDFLAGS) $(PETSC_LDFLAGS) \
		          $(HIPGISAXS_LDFLAGS) $(HDF5_LDFLAGS) $(TIFF_LDFLAGS) $(BOOST_LDFLAGS)
ALL_LIBS = $(MPI_LIBS) $(TAO_LIBS) $(PETSC_LIBS) $(HIPGISAXS_LIBS) $(HDF5_LIBS) $(TIFF_LIBS) $(BOOST_LIBS)

OBJECTS = ImageData.o analysis_algorithm.o objective_func.o hipgisaxs_ana.o objective_func_hipgisaxs.o

POUND_OBJS = hipgisaxs_fit_pounders.o
POUND_HIP_OBJS = $(POUND_OBJS) hipgisaxs_fit_pounders_main.o

LMVM_OBJS = hipgisaxs_fit_lmvm.o
LMVM_HIP_OBJS = $(LMVM_OBJS) hipgisaxs_fit_lmvm_main.o

F1_OBJS = objective_func_poly_one.o
F1_HIP_OBJS = $(F1_OBJS) objective_func_poly_one_main.o

PSO_OBJS = hipgisaxs_fit_pso.o hipgisaxs_fit_pso_particle.o
PSO_HIP_OBJS = $(PSO_OBJS) hipgisaxs_fit_pso_main.o

BF_OBJS = hipgisaxs_fit_bruteforce.o
BF_HIP_OBJS = $(BF_OBJS) hipgisaxs_fit_bruteforce_main.o

HIP_OBJS = hipgisaxs_fit_main.o $(POUND_OBJS) $(PSO_OBJS) $(BF_OBJS) $(LMVM_OBJS)

TEST_HIP_OBJS = hipgisaxs_fit_test.o $(POUND_OBJS) $(PSO_OBJS) $(BF_OBJS) $(LMVM_OBJS)

BINARY = analyze

all: hipgisaxs

hipgisaxs:
	(make $(BINARY) \
		"OBJECTS = $(OBJECTS) $(HIP_OBJS)")

test:
	(make $(BINARY) \
		"OBJECTS = $(OBJECTS) $(TEST_HIP_OBJS)")

pounders:
	(make $(BINARY) \
		"OBJECTS = $(OBJECTS) $(POUND_HIP_OBJS)")

f1:
	(make $(BINARY) \
		"OBJECTS = $(OBJECTS) $(F1_HIP_OBJS)")

pso:
	(make $(BINARY) \
		"OBJECTS = $(OBJECTS) $(PSO_HIP_OBJS)")

bf:
	(make $(BINARY) \
		"OBJECTS = $(OBJECTS) $(BF_HIP_OBJS)")

lmvm:
	(make $(BINARY) \
		"OBJECTS = $(OBJECTS) $(LMVM_HIP_OBJS)")

## c++ compilation                                                                                                                             
%.o: $(SRC_DIR)/%.cpp
	$(CXX) -c $< -o $@  $(ALL_FLAGS) $(ALL_INCL) 

$(BINARY): $(OBJECTS)
	$(CXX) $(ALL_FLAGS) -o $(BIN_DIR)/$@ $^ $(ALL_LDFLAGS) $(ALL_LIBS)

clean:
	rm -f *.o $(BINARY) $(HIPGISAXS_DIR)/bin/$(BINARY)

install: $(BINARY)
	@install -v $(BINARY) $(HIPGISAXS_DIR)/bin

error: errormsg

errormsg:
	@echo ""
	@echo "usage: $(MAKE) <target> [options]"
	@echo ""
	@echo "target: some examples are"
	@echo "          andromeda-intel"
	@echo "          andromeda-gnu"
	@echo "          bragg-gnu"
	@echo "          nersc-intel"
	@echo "          nersc-gnu"
	@echo "          titan-gnu"
	@echo ""
	@echo "options: none yet."
	@echo ""
	@echo "following environment variables may also be defined:"
	@echo "          HIPGISAXS_DIR"
	@echo "          TIFF_DIR"
	@echo "          BOOST_DIR"
	@echo "          PETSC_VERSION"
	@echo "          PETSC_DIR"
	@echo "          MPI_DIR"
	@echo "          CUDA_DIR"
	@echo "          HDF5_DIR"
	@echo "          SZIP_DIR"
	@echo ""

