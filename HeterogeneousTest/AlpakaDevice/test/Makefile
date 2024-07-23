# Compiler and flags
CXX := $(shell scram tool tag icx-cxxcompiler CXX)
GCC_TOOLCHAIN := $(shell scram tool tag gcc-cxxcompiler GCC_CXXCOMPILER_BASE)
TARGET := x86_64-redhat-linux-gnu
CXXFLAGS := -fsycl --std=c++17 -fsycl-targets=spir64_x86_64 --target=$(TARGET) --gcc-toolchain=$(GCC_TOOLCHAIN)

# Include paths from external tools
# clone alpaka from https://github.com/AuroraPerego/alpaka/tree/SYCLqueue
ALPAKA_BASE := /data/user/aperego/alpaka
ALPAKA_FLAGS := -DALPAKA_ACC_SYCL_ENABLED -DALPAKA_SYCL_ONEAPI_CPU -DALPAKA_DISABLE_VENDOR_RNG
ALPAKA := -I$(ALPAKA_BASE)/include $(ALPAKA_FLAGS) -frelaxed-template-template-args
# -frelaxed-template-template-args needed with alpaka 1.2.0
# should  become $(shell scram tool tag alpaka INCLUDE)

CATCH2 := $(shell scram tool tag catch2 INCLUDE)
BOOST := $(shell scram tool tag boost BOOST_BASE)/include

INCLUDES := -I$(CMSSW_BASE)/src -I$(CMSSW_RELEASE_BASE)/src -I$(CATCH2)-I$(BOOST) -I$(ALPAKA)

# Source files
TEST_DIR := $(CMSSW_BASE)/src/HeterogeneousTest/AlpakaDevice/test/alpaka

SRC := $(TEST_DIR)/testDeviceAddition.dev.cc

# Output binary
EXE := AlpakaTest

# Default target
all: $(EXE)

# Link object files into the final binary
$(EXE): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Wno-deprecated-declarations -Wno-unknown-attributes  $< -o $@

# Clean target
clean:
	rm -f $(EXE)

# Phony targets
.PHONY: all clean
