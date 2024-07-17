CXX := icpx
GCC_TOOLCHAIN := /cvmfs/cms-ib.cern.ch/sw/x86_64/nweek-02845/el8_amd64_gcc12/external/gcc/12.3.1-40d504be6370b5a30e3947a6e575ca28
TARGET := x86_64-redhat-linux-gnu
CXXFLAGS := -fsycl --std=c++17 -fsycl-targets=spir64_x86_64 --target=$(TARGET) --gcc-toolchain=$(GCC_TOOLCHAIN) -Wno-unknown-attributes

# Include paths from external tools
CATCH2 := $(shell scram tool tag catch2 INCLUDE)
ROOT := $(shell scram tool tag root_interface ROOT_INTERFACE_BASE)/include
TBB := $(shell scram tool tag tbb INCLUDE)
BOOST := $(shell scram tool tag boost BOOST_BASE)/include
FMT := $(shell scram tool tag fmt INCLUDE)

INCLUDES := -I$(CMSSW_BASE)/src -I$(CMSSW_RELEASE_BASE)/src -I$(CATCH2) -I$(ROOT) -I$(TBB) -I$(BOOST) -I$(FMT)

# Source files
SRC_DIR := $(CMSSW_BASE)/src/HeterogeneousTest/SYCLDevice/src
PLUGINS_DIR := $(CMSSW_BASE)/src/HeterogeneousTest/SYCLDevice/plugins

# Object files
OBJ_DIR = obj

# Create the object directory if it doesn't exist
$(shell mkdir -p $(OBJ_DIR))

# Object files
OBJ_FILES := $(OBJ_DIR)/SYCLTestDeviceAdditionAlgo.sycl.o \
             $(OBJ_DIR)/SYCLTestDeviceAdditionModule.sycl.o \
             $(OBJ_DIR)/DeviceAddition.sycl.o

SHARED_LIB := SYCLTestDeviceAdditionModule.so

# Targets
all: $(SHARED_LIB)

# Compile source files into object files
$(OBJ_DIR)/%.sycl.o: $(SRC_DIR)/%.sycl.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -fPIC $< -o $@

$(OBJ_DIR)/%.sycl.o: $(PLUGINS_DIR)/%.sycl.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -fPIC $< -o $@

$(SHARED_LIB): $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -shared -o $@

clean:
	rm -f $(OBJ_FILES) $(SHARED_LIB)

.PHONY: all clean
