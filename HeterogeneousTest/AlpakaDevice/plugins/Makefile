# Compiler and flags
CXX := $(shell scram tool tag gcc-cxxcompiler CXX)
SYCL_INCLUDE := $(shell scram tool tag icx-cxxcompiler ICX_CXXCOMPILER_BASE)
CXXFLAGS := --std=c++17 -O3 -I$(SYCL_INCLUDE)/include -I$(SYCL_INCLUDE)/include/sycl -Wno-attributes

SYCL_CXX := $(shell scram tool tag icx-cxxcompiler CXX)
GCC_TOOLCHAIN := $(shell scram tool tag gcc-cxxcompiler GCC_CXXCOMPILER_BASE)
TARGET := x86_64-redhat-linux-gnu
SYCL_CXXFLAGS := -fsycl --std=c++17 -fsycl-targets=spir64_x86_64,intel_gpu_acm_g10,intel_gpu_pvc --target=$(TARGET) --gcc-toolchain=$(GCC_TOOLCHAIN)

COMMON_FLAGS := -Wno-deprecated-declarations -c -fPIC
LINK_FLAGS := -Wno-deprecated-declarations -frelaxed-template-template-args -shared

# Include paths from external tools
#
# clone alpaka from https://github.com/AuroraPerego/alpaka/tree/SYCLqueue
# should  become $(shell scram tool tag alpaka INCLUDE)
ALPAKA_BASE := /data/user/aperego/alpaka
ALPAKA_CPU_FLAGS := -DALPAKA_ACC_SYCL_ENABLED -DALPAKA_SYCL_ONEAPI_CPU -DALPAKA_DISABLE_VENDOR_RNG
ALPAKA_GPU_FLAGS := -DALPAKA_ACC_SYCL_ENABLED -DALPAKA_SYCL_ONEAPI_GPU -DALPAKA_DISABLE_VENDOR_RNG
ALPAKA := $(ALPAKA_BASE)/include

CATCH2 := $(shell scram tool tag catch2 INCLUDE)
ROOT := $(shell scram tool tag root_interface ROOT_INTERFACE_BASE)/include
TBB := $(shell scram tool tag tbb INCLUDE)
BOOST := $(shell scram tool tag boost BOOST_BASE)/include
FMT := $(shell scram tool tag fmt INCLUDE)
BOOST := $(shell scram tool tag boost BOOST_BASE)/include

INCLUDES := -I$(CMSSW_BASE)/src -I$(CMSSW_RELEASE_BASE)/src -I$(CATCH2) -I$(ROOT) -I$(TBB) -I$(BOOST) -I$(FMT) -I$(ALPAKA)

# Source files
PLUGINS_DIR := $(CMSSW_BASE)/src/HeterogeneousTest/AlpakaDevice/plugins/alpaka

SRC := $(PLUGINS_DIR)/AlpakaTestDeviceAdditionAlgo.dev.cc \
			 $(PLUGINS_DIR)/AlpakaTestDeviceAdditionModule.cc

# Object files
OBJ_DIR = obj

# Create the object directory if it doesn't exist
$(shell mkdir -p $(OBJ_DIR))

# Object files
OBJ_FILES := $(OBJ_DIR)/AlpakaTestDeviceAdditionAlgo.dev.syclcpu.o \
             $(OBJ_DIR)/AlpakaTestDeviceAdditionModule.syclcpu.o \
             $(OBJ_DIR)/AlpakaTestDeviceAdditionAlgo.dev.syclgpu.o \
             $(OBJ_DIR)/AlpakaTestDeviceAdditionModule.syclgpu.o

# Output library
SHARED_LIB := AlpakaTestDeviceAdditionModule.so

# Default target
all: $(SHARED_LIB)

$(OBJ_DIR)/%.dev.syclcpu.o: $(PLUGINS_DIR)/%.dev.cc
	$(SYCL_CXX) $(SYCL_CXXFLAGS) $(INCLUDES) $(ALPAKA_CPU_FLAGS) $(COMMON_FLAGS) -frelaxed-template-template-args $< -o $@
$(OBJ_DIR)/%.dev.syclgpu.o: $(PLUGINS_DIR)/%.dev.cc
	$(SYCL_CXX) $(SYCL_CXXFLAGS) $(INCLUDES) $(ALPAKA_GPU_FLAGS) $(COMMON_FLAGS) -frelaxed-template-template-args $< -o $@
# -frelaxed-template-template-args needed with alpaka 1.2.0

$(OBJ_DIR)/%.syclcpu.o: $(PLUGINS_DIR)/%.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(ALPAKA_CPU_FLAGS) $(COMMON_FLAGS) $< -o $@
$(OBJ_DIR)/%.syclgpu.o: $(PLUGINS_DIR)/%.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(ALPAKA_GPU_FLAGS) $(COMMON_FLAGS) $< -o $@

# Link object files into the final binary
$(SHARED_LIB): $(OBJ_FILES)
	$(SYCL_CXX) $(SYCL_CXXFLAGS) $(INCLUDES) $(LINK_FLAGS) $^ -o $@

# Clean target
clean:
	rm -rf $(OBJ_DIR) $(SHARED_LIB)

# Phony targets
.PHONY: all clean
