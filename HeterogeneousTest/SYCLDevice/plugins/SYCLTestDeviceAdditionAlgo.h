#ifndef HeterogeneousTest_SYCLDevice_plugins_CUDATestDeviceAdditionAlgo_h
#define HeterogeneousTest_SYCLDevice_plugins_CUDATestDeviceAdditionAlgo_h

#include <cstddef>

namespace HeterogeneousTestSYCLDevicePlugins {

  void wrapper_add_vectors_f(const float* __restrict__ in1,
                             const float* __restrict__ in2,
                             float* __restrict__ out,
                             size_t size,
                             sycl::queue& queue);

}  // namespace HeterogeneousTestSYCLDevicePlugins

#endif  // HeterogeneousTest_SYCLDevice_plugins_CUDATestDeviceAdditionAlgo_h
