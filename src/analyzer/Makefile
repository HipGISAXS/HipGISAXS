
MPI_DIR = /usr/local/openmpi
TAO_DIR = /usr/local/petsc-3.4.3-gnu
PETSC_DIR = /usr/local/petsc-3.4.3-gnu
HIPGISAXS_DIR = /home/asarje/hipgisaxs.git
HDF5_DIR = /usr/local/hdf5-1.8.9
TIFF_DIR = /usr/local/tiff-4.0.2
BOOST_DIR = /usr/local/boost
SZIP_DIR = /usr/local/szip-2.1
CUDA_DIR = /usr/local/cuda


PREFIX = $(PWD)
BINARY = analyze
INCL_DIR =$(PREFIX)
LIB_DIR = $(PREFIX)
SRC_DIR = $(PREFIX)
BIN_DIR = $(PREFIX)
OBJ_DIR = $(PREFIX)
CXX = mpicxx
CXX_FLAGS = -std=c++0x -fopenmp
CXX_OPT_FLAGS = -g
PREC_FLAG =

########EXTERNAL DEPENDENCIES###########
MPI_LDFLAGS = -L $(MPI_DIR)/lib
MPI_LIBS = -lmpi -lmpi_cxx

TAO_INCL = -I $(TAO_DIR)/include
TAO_LDFLAGS = -L $(TAO_DIR)/lib
TAO_LIBS = -ltao

#GCC_LDFLAGS = -L $(GCC_DIR)/lib/gcc/x86_64-unknown-linux-gnu/4.7.0                                                                                                 
#GCC_LDFLAGS += -L $(GCC_DIR)/lib/gcc -L $(GCC_DIR)/lib64 -L $(GCC_DIR)/lib                                                                                         

PETSC_INCL = -I $(PETSC_DIR)/include
PETSC_LDFLAGS = -L $(PETSC_DIR)/lib -Wl,-rpath,$(PETSC_DIR)/lib #-L /usr/syscom/opt/torque/4.1.4/lib                                         
#PETSC_LIBS = -lpetsc -lX11 -lpthread -lf2clapack -lf2cblas -lstdc++ -ldl -lrdmacm                                            
PETSC_LIBS += -lpetsc -lf2clapack -lf2cblas -lsuperlu_4.3
PETSC_LIBS += -lX11 -lrt -lnsl -lutil -lgcc_s -u dlacon_
#PETSC_LIBS += -libverbs                                                                                                                                            
#PETSC_LIBS += -lrt -lnsl -lutil -lgcc_s #-ltorque                                                                                                                  

HIPGISAXS_INCL = -I $(HIPGISAXS_DIR)/include -I $(CUDA_DIR)/include
HIPGISAXS_LDFLAGS = -L $(HIPGISAXS_DIR)/lib -L $(CUDA_DIR)/lib64
HIPGISAXS_LIBS = -lhipgisaxs -lcudart

HDF5_INCL = -I $(HDF5_DIR)/include
HDF5_LDFLAGS = -L $(HDF5_DIR)/lib -L $(SZIP_DIR)/lib
HDF5_LIBS = -lhdf5 -lz -lsz

TIFF_INCL = -I $(TIFF_DIR)/include
TIFF_LDFLAGS = -L $(TIFF_DIR)/lib
TIFF_LIBS = -ltiff

BOOST_INCL = -I $(BOOST_DIR)/include -I $(NUM_BOOST_DIR)/include
BOOST_LDFLAGS = -L $(BOOST_DIR)/lib
BOOST_LIBS = -lboost_system -lboost_filesystem -lboost_timer -lboost_chrono

#ALL_INCL = $(TAO_INCL) $(PETSC_INCL) $(HIPGISAXS_INCL) $(HDF5_INCL) $(TIFF_INCL) $(BOOST_INCL)
ALL_LDFLAGS = $(MPI_LDFLAGS) $(TAO_LDFLAGS) $(GCC_LDFLAGS) $(PETSC_LDFLAGS) \
                          $(HIPGISAXS_LDFLAGS) $(HDF5_LDFLAGS) $(TIFF_LDFLAGS) $(BOOST_LDFLAGS)
ALL_LIBS = $(MPI_LIBS) $(TAO_LIBS) $(PETSC_LIBS) $(HIPGISAXS_LIBS) $(HDF5_LIBS) $(TIFF_LIBS) $(BOOST_LIBS)

#######################################

## all objs and libs                                                                                                                           
ALL_INCL= -Iinclude $(TAO_INCL) $(PETSC_INCL) $(HDF5_INCL) $(TIFF_INCL) $(BOOST_INCL) $(HIPGISAXS_INCL)
ALL_OBJ = Data.o  ImageData.o ImageDataProc.o AnaAlgorithm.o  FitPOUNDERSAlgo.o Workflow.o Dist.o  StructVar.o AnaInput.o \
	AnaOutput.o Analysis.o Analyzer.o ObjFct.o SimModel.o SimHipGISAXS.o SimTest.o

#ALL_LIB = libData.a libImageData.o libImageDataProc.o libAnaAlgorithm.a libFitPOUNDERSAlgo.a libWorkflow.a libDist.a libStructVar.a libAnaInput.a \
	libAnaOutput.a libAnalysis.a libAnalyzer.a libObjFct.a libSimModel.a libSimHipGISAXS.a

APP_OBJS = $(patsubst %, $(OBJ_DIR)/%, $(ALL_OBJ))
APP_LIBS = $(ALL_LIBS) # $(patsubst lib%.a, -l%, $(ALL_LIB))                                                                                               

## c++ compilation                                                                                                                             
_DEPS_CXX = $(HIPGISAXS_DIR)/include/analyzer/%.hpp
DEPS_CXX = $(patsubst %,$(SRC_DIR)/%,$(_DEPS_CXX))

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp     ###$(DEPS_CXX)                                                                                            
	$(CXX) -c $< -o $@  $(CXX_FLAGS) $(PREC_FLAG) $(ALL_INCL) 


$(BINARY): $(APP_OBJS)
	$(CXX)  $(CXX_FLAGS) -o $(BIN_DIR)/$@ $^ -Llib $(ALL_LDFLAGS) $(APP_LIBS) $(ALL_INCL) $(SRC_DIR)/ana_main.cpp  

clean:
	rm -f *.o $(BINARY)
