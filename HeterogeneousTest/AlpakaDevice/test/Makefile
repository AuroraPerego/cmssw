# Compiler and flags
CXX := icpx
GCC_TOOLCHAIN := /cvmfs/cms-ib.cern.ch/sw/x86_64/nweek-02845/el8_amd64_gcc12/external/gcc/12.3.1-40d504be6370b5a30e3947a6e575ca28
TARGET := x86_64-redhat-linux-gnu
CXXFLAGS := -fsycl --std=c++17 -fsycl-targets=spir64_x86_64 --target=$(TARGET) --gcc-toolchain=$(GCC_TOOLCHAIN)
CATCH2 := /cvmfs/cms.cern.ch/el8_amd64_gcc12/external/catch2/2.13.6-17102db92de47c6a473c6e67627c548a/include
INCLUDES := -I$(CMSSW_BASE)/src -I$(CMSSW_RELEASE_BASE)/src -I$(CATCH2)

# clone alpaka from https://github.com/AuroraPerego/alpaka/tree/SYCLqueue
ALPAKA_BASE := /data/user/aperego/alpaka
ALPAKA := -I$(ALPAKA_BASE)/include
# should  become $(shell scram tool tag alpaka INCLUDE)

BOOST := -I$(shell scram tool tag boost BOOST_BASE)/include
ALPAKA_FLAGS := -DALPAKA_ACC_SYCL_ENABLED -DALPAKA_SYCL_ONEAPI_CPU -DALPAKA_DISABLE_VENDOR_RNG

# Source files
TEST_DIR := $(CMSSW_BASE)/src/HeterogeneousTest/AlpakaDevice/test/alpaka

SRC := $(TEST_DIR)/testDeviceAddition.dev.cc

# Output binary
EXE := AlpakaTest

# Default target
all: $(EXE)

# Link object files into the final binary
$(EXE): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(ALPAKA) $(BOOST) $(ALPAKA_FLAGS) -Wno-deprecated-declarations -Wno-unknown-attributes  $< -o $@

# Clean target
clean:
	rm -f $(EXE)

# Phony targets
.PHONY: all clean
