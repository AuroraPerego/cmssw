#ifndef HeterogeneousTest_SYCLDevice_interface_DeviceAddition_h
#define HeterogeneousTest_SYCLDevice_interface_DeviceAddition_h

#include <cstddef>

#include <sycl/sycl.hpp>

namespace cms::sycltest {

  SYCL_EXTERNAL void add_vectors_f(const float* __restrict__ in1,
                                   const float* __restrict__ in2,
                                   float* __restrict__ out,
                                   size_t size);

  SYCL_EXTERNAL void add_vectors_d(const double* __restrict__ in1,
                                   const double* __restrict__ in2,
                                   double* __restrict__ out,
                                   size_t size);

}  // namespace cms::sycltest

#endif  // HeterogeneousTest_SYCLDevice_interface_DeviceAddition_h
