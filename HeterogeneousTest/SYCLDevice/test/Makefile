# Compiler and flags
CXX := icpx
GCC_TOOLCHAIN := /cvmfs/cms-ib.cern.ch/sw/x86_64/nweek-02845/el8_amd64_gcc12/external/gcc/12.3.1-40d504be6370b5a30e3947a6e575ca28
TARGET := x86_64-redhat-linux-gnu
CXXFLAGS := -fsycl --std=c++17 -fsycl-targets=spir64_x86_64 --target=$(TARGET) --gcc-toolchain=$(GCC_TOOLCHAIN)
CATCH2 := /cvmfs/cms.cern.ch/el8_amd64_gcc12/external/catch2/2.13.6-17102db92de47c6a473c6e67627c548a/include
INCLUDES := -I$(CMSSW_BASE)/src -I$(CATCH2)

# Source files
SRC_DIR := $(CMSSW_BASE)/src/HeterogeneousTest/SYCLDevice/src
TEST_DIR := $(CMSSW_BASE)/src/HeterogeneousTest/SYCLDevice/test
UTIL_DIR := $(CMSSW_BASE)/src/HeterogeneousCore/SYCLUtilities/src

SRCS := $(SRC_DIR)/DeviceAddition.cc \
       $(TEST_DIR)/testDeviceAddition.cc \
       $(UTIL_DIR)/requireDevices.cc

# Object files
OBJ_DIR = obj

# Create the object directory if it doesn't exist
$(shell mkdir -p $(OBJ_DIR))

OBJ_FILES = $(OBJ_DIR)/DeviceAddition.cc.o \
            $(OBJ_DIR)/testDeviceAddition.cc.o \
            $(OBJ_DIR)/requireDevices.cc.o

# Output binary
EXE := SYCLtest

# Default target
all: $(EXE)

# Compile source files into object files
$(OBJ_DIR)/%.cc.o: $(SRC_DIR)/%.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.cc.o: $(TEST_DIR)/%.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.cc.o: $(UTIL_DIR)/%.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Link object files into the final binary
$(EXE): $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

# Clean target
clean:
	rm -f $(OBJ_FILES) $(EXE)

# Phony targets
.PHONY: all clean
