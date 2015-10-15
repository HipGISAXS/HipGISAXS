
#PETSC_VERSION = "old"
PETSC_VERSION = "new"

#MPI_DIR = /usr/local/openmpi
#PETSC_DIR = /usr/local/petsc-3.4.3-gnu
HIPGISAXS_DIR = ../..
#HDF5_DIR = /usr/local/hdf5-1.8.9
TIFF_DIR = $(HOME)/opt/tiff-4.0.2
#BOOST_DIR = /usr/local/boost
#SZIP_DIR = /usr/local/szip-2.1
CUDA_DIR = $(CRAY_CUDATOOLKIT_DIR)
#NUM_BOOST_DIR = $(HOME)/opt/boost-gil-numeric

PREFIX = $(PWD)
INCL_DIR =$(PREFIX)
LIB_DIR = $(PREFIX)
SRC_DIR = $(PREFIX)
BIN_DIR = $(PREFIX)
OBJ_DIR = $(PREFIX)

CXX = CC
CXX_FLAGS = -std=c++11 -fopenmp -DUSE_MPI
CXX_FLAGS += -O3 #-g
PREC_FLAG = -DDOUBLEP

ifeq ($(PETSC_VERSION), "new")
  CXX_FLAGS += -DPETSC_35
  TAO_DIR =
  TAO_INCL =
  TAO_LDFLAGS =
  TAO_LIBS =
else
  TAO_DIR = $(HOME)/opt/tao-2.2.0/lib64/python2.6/site-packages/tao
  TAO_INCL = -I $(TAO_DIR)/include
  TAO_LDFLAGS = -L $(TAO_DIR)/lib
  TAO_LIBS = -ltao
endif

#MPI_LDFLAGS = -L $(MPI_DIR)/lib
#MPI_LIBS = -lmpi -lmpi_cxx


#PETSC_INCL = -I $(PETSC_DIR)/include
PETSC_LDFLAGS = -L $(PETSC_DIR)/lib #-Wl,-rpath,$(PETSC_DIR)/lib #-L /usr/syscom/opt/torque/4.1.4/lib
#PETSC_LIBS = -lpetsc -lX11 -lpthread -lf2clapack -lf2cblas -lstdc++ -ldl -lrdmacm
#PETSC_LIBS += -lpetsc #-lf2clapack -lf2cblas -lsuperlu_4.3
#PETSC_LIBS += -lX11 -lrt -lnsl -lutil -lgcc_s -u dlacon_
#PETSC_LIBS += -libverbs
#PETSC_LIBS += -lrt -lnsl -lutil -lgcc_s #-ltorque
PETSC_LIBS = -lcraypetsc_gnu_real

HIPGISAXS_INCL = -I $(HIPGISAXS_DIR)/include -I $(CUDA_DIR)/include
HIPGISAXS_LDFLAGS = -L $(HIPGISAXS_DIR)/lib -L $(CUDA_DIR)/lib64
HIPGISAXS_LIBS = -lhipgisaxs -lcudart

#HDF5_INCL = -I $(HDF5_DIR)/include
#HDF5_LDFLAGS = -L $(HDF5_DIR)/lib -L $(SZIP_DIR)/lib
#HDF5_LIBS = -lhdf5 -lz #-lsz

TIFF_INCL = -I $(TIFF_DIR)/include
TIFF_LDFLAGS = -L $(TIFF_DIR)/lib
TIFF_LIBS = -ltiff

BOOST_INCL = -I $(BOOST_DIR)/include #-I $(NUM_BOOST_DIR)/include
BOOST_LDFLAGS = -L $(BOOST_DIR)/lib
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

test:
	(make $(BINARY) \
		"OBJECTS = $(OBJECTS) $(TEST_HIP_OBJS)")

hipgisaxs:
	(make $(BINARY) \
		"OBJECTS = $(OBJECTS) $(HIP_OBJS)")

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
	$(CXX) -c $< -o $@  $(CXX_FLAGS) $(PREC_FLAG) $(ALL_INCL) 

$(BINARY): $(OBJECTS)
	$(CXX)  $(CXX_FLAGS) $(PREC_FLAG) -o $(BIN_DIR)/$@ $^ $(ALL_LDFLAGS) $(ALL_LIBS)

clean:
	rm -f *.o $(BINARY)

install: $(BINARY)
	@install -v $(BINARY) $(HIPGISAXS_DIR)/bin
