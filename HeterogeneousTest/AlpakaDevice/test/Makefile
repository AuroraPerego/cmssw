# Compiler and flags
CXX := icpx
GCC_TOOLCHAIN := /cvmfs/cms-ib.cern.ch/sw/x86_64/nweek-02845/el8_amd64_gcc12/external/gcc/12.3.1-40d504be6370b5a30e3947a6e575ca28
TARGET := x86_64-redhat-linux-gnu
CXXFLAGS := -fsycl --std=c++17 -fsycl-targets=spir64_x86_64 --target=$(TARGET) --gcc-toolchain=$(GCC_TOOLCHAIN)
CATCH2 := /cvmfs/cms.cern.ch/el8_amd64_gcc12/external/catch2/2.13.6-17102db92de47c6a473c6e67627c548a/include
INCLUDES := -I$(CMSSW_BASE)/src -I$(CMSSW_RELEASE_BASE)/src -I$(CATCH2)
ALPAKA := -I/data/user/aperego/alpaka/include
#$(scram tool tag alpaka INCLUDE)
BOOST := -I$(shell scram tool tag boost BOOST_BASE)/include
ALPAKA_FLAGS := -DALPAKA_ACC_SYCL_ENABLED -DALPAKA_SYCL_ONEAPI_CPU -DALPAKA_DISABLE_VENDOR_RNG


# icpx \
#   -fsycl \
# 	--gcc-toolchain=/cvmfs/cms-ib.cern.ch/sw/x86_64/nweek-02845/el8_amd64_gcc12/external/gcc/12.3.1-40d504be6370b5a30e3947a6e575ca28 \
# 	--target=x86_64-redhat-linux-gnu \
# 	--std=c++17 \
# 	-fsycl-targets=spir64_x86_64 \
# 	-I/cvmfs/cms.cern.ch/el8_amd64_gcc12/external/catch2/2.13.6-17102db92de47c6a473c6e67627c548a/include \
# 	-DALPAKA_ACC_SYCL_ENABLED \
# 	-DALPAKA_SYCL_ONEAPI_CPU \
# 	-I/cvmfs/cms-ib.cern.ch/sw/x86_64/nweek-02845/el8_amd64_gcc12/external/boost/1.80.0-a1544032d9d65904ac2112b6d35bba55/include \
# 	-I/data/user/aperego/alpaka/include  HeterogeneousTest/AlpakaDevice/test/alpaka/testDeviceAddition.dev.cc  \
# 	-o AlpakaTest \
# 	-DALPAKA_DISABLE_VENDOR_RNG \
# 	-I/data/user/aperego/cmssw/CMSSW_14_1_X_2024-07-10-2300/src \
# 	-I/cvmfs/cms-ib.cern.ch/sw/x86_64/week1/el8_amd64_gcc12/cms/cmssw/CMSSW_14_1_X_2024-07-10-2300/src \
# 	-Wno-deprecated-declarations \
# 	-Wno-unknown-attributes

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
