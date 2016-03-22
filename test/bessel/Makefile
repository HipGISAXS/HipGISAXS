HIPGISAXS_DIR = ../..
BOOST_DIR = /usr/local/boost
CUDA_DIR = /usr/local/cuda

PREFIX = $(PWD)
INCL_DIR = $(PREFIX)
LIB_DIR = $(PREFIX)
SRC_DIR = $(PREFIX)
BIN_DIR = $(PREFIX)
OBJ_DIR = $(PREFIX)

PREC_FLAG = -DDOUBLEP -DINTEL_AVX -DFF_CPU_OPT -DFF_CPU_OPT_AVX

CXX = icpc
CXX_FLAGS = -std=c++11 -openmp
#CXX_FLAGS += -O3
CXX_FLAGS += -g

HIPGISAXS_INCL = -I $(HIPGISAXS_DIR)/include
HIPGISAXS_LDFLAGS = -L $(HIPGISAXS_DIR)/lib
HIPGISAXS_LIBS = -lhipgisaxs

BOOST_INCL = -I $(BOOST_DIR)/include
BOOST_LDFLAGS = -L $(BOOST_DIR)/lib
BOOST_LIBS = -lboost_system -lboost_filesystem -lboost_timer -lboost_chrono


## all objs and libs                                                                                                                           
ALL_INCL= $(BOOST_INCL) $(HIPGISAXS_INCL)
ALL_LDFLAGS = $(HIPGISAXS_LDFLAGS) $(BOOST_LDFLAGS)
ALL_LIBS = $(HIPGISAXS_LIBS) $(BOOST_LIBS)

CBESSELJ = cbesselj
CBESSELJ_OBJ = cbesselj.o

UNIQUIFY = uniquify
UNIQUIFY_OBJ = uniquify.o

BINARIES = $(CBESSELJ) $(UNIQUIFY)

all: $(BINARIES)

$(CBESSELJ): $(CBESSELJ_OBJ)
	$(CXX) $(CXX_FLAGS) $(PREC_FLAG) -o $@ $^ $(ALL_LDFLAGS) $(ALL_LIBS)

$(UNIQUIFY): $(UNIQUIFY_OBJ)
	$(CXX) $(CXX_FLAGS) $(PREC_FLAG) -o $@ $^ $(ALL_LDFLAGS) $(ALL_LIBS)

## c++ compilation                                                                                                                             
%.o: %.cpp
	$(CXX) -c $< -o $@  $(CXX_FLAGS) $(PREC_FLAG) $(ALL_INCL) 

clean:
	rm -f *.o $(BINARIES)
