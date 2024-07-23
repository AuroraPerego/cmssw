# Compiler and flags
CXX := $(shell scram tool tag icx-cxxcompiler CXX)
GCC_TOOLCHAIN := $(shell scram tool tag gcc-cxxcompiler GCC_CXXCOMPILER_BASE)
TARGET := x86_64-redhat-linux-gnu
CXXFLAGS := -fsycl --std=c++17 -fsycl-targets=spir64_x86_64,intel_gpu_acm_g10,intel_gpu_pvc --target=$(TARGET) --gcc-toolchain=$(GCC_TOOLCHAIN)

# Include paths from external tools
CATCH2 := $(shell scram tool tag catch2 INCLUDE)

INCLUDES := -I$(CMSSW_BASE)/src -I$(CATCH2)

# Source files
SRC_DIR := $(CMSSW_BASE)/src/HeterogeneousTest/SYCLDevice/src
TEST_DIR := $(CMSSW_BASE)/src/HeterogeneousTest/SYCLDevice/test
UTIL_DIR := $(CMSSW_BASE)/src/HeterogeneousCore/SYCLUtilities/src

SRCS := $(SRC_DIR)/DeviceAddition.sycl.cc \
       $(TEST_DIR)/testDeviceAddition.sycl.cc \
       $(UTIL_DIR)/requireDevices.sycl.cc

# Object files
OBJ_DIR = obj

# Create the object directory if it doesn't exist
$(shell mkdir -p $(OBJ_DIR))

OBJ_FILES = $(OBJ_DIR)/DeviceAddition.sycl.cc.o \
            $(OBJ_DIR)/testDeviceAddition.sycl.cc.o \
            $(OBJ_DIR)/requireDevices.sycl.cc.o

# Output binary
EXE := SYCLtest

# Default target
all: $(EXE)

# Compile source files into object files
$(OBJ_DIR)/%.sycl.cc.o: $(SRC_DIR)/%.sycl.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.sycl.cc.o: $(TEST_DIR)/%.sycl.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.sycl.cc.o: $(UTIL_DIR)/%.sycl.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Link object files into the final binary
$(EXE): $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@

# Clean target
clean:
	rm -rf $(OBJ_DIR) $(EXE)

# Phony targets
.PHONY: all clean
