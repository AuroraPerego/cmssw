#include <cstddef>
#include <cstdint>

#include <sycl/sycl.hpp>

#include "HeterogeneousTest/SYCLDevice/interface/DeviceAddition.h"

namespace cms::sycltest {

  void add_vectors_f(const float* __restrict__ in1,
                     const float* __restrict__ in2,
                     float* __restrict__ out,
                     size_t size) {
    auto item = sycl::ext::oneapi::experimental::this_nd_item<1>();
    uint32_t thread = item.get_global_linear_id();
    uint32_t stride = item.get_local_range(0) * item.get_group_range(0);

    for (size_t i = thread; i < size; i += stride) {
      out[i] = in1[i] + in2[i];
    }
  }

  void add_vectors_d(const double* __restrict__ in1,
                     const double* __restrict__ in2,
                     double* __restrict__ out,
                     size_t size) {
    auto item = sycl::ext::oneapi::experimental::this_nd_item<1>();
    uint32_t thread = item.get_global_linear_id();
    uint32_t stride = item.get_local_range(0) * item.get_group_range(0);

    for (size_t i = thread; i < size; i += stride) {
      out[i] = in1[i] + in2[i];
    }
  }

}  // namespace cms::sycltest
