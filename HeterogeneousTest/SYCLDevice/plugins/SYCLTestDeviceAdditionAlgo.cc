#include <cstddef>

#include <sycl/sycl.hpp>

#include "HeterogeneousTest/SYCLDevice/interface/DeviceAddition.h"

#include "SYCLTestDeviceAdditionAlgo.h"

namespace HeterogeneousTestSYCLDevicePlugins {

  void kernel_add_vectors_f(const float* __restrict__ in1,
                            const float* __restrict__ in2,
                            float* __restrict__ out,
                            size_t size) {
    cms::sycltest::add_vectors_f(in1, in2, out, size);
  }

  void wrapper_add_vectors_f(const float* __restrict__ in1,
                             const float* __restrict__ in2,
                             float* __restrict__ out,
                             size_t size,
                             sycl::queue& queue) {
    queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(sycl::nd_range<1>{32 * 32, 32},
                       [=](sycl::nd_item<1> idx) { kernel_add_vectors_f(in1, in2, out, size); });
    });
  }

}  // namespace HeterogeneousTestSYCLDevicePlugins
