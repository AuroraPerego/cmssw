#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <sycl/sycl.hpp>

#include "HeterogeneousCore/SYCLUtilities/interface/requireDevices.h"
#include "HeterogeneousTest/SYCLDevice/interface/DeviceAddition.h"

void kernel_add_vectors_f(const float* __restrict__ in1,
                          const float* __restrict__ in2,
                          float* __restrict__ out,
                          size_t size) {
  cms::sycltest::add_vectors_f(in1, in2, out, size);
}

TEST_CASE("HeterogeneousTest/SYCLDevice test", "[syclTestDeviceAddition]") {
  cms::sycltest::requireDevices();

  sycl::queue queue{sycl::cpu_selector_v, sycl::property::queue::in_order()};

  // random number generator with a gaussian distribution
  std::random_device rd{};
  std::default_random_engine rand{rd()};
  std::normal_distribution<float> dist{0., 1.};

  // tolerance
  constexpr float epsilon = 0.000001;

  // buffer size
  constexpr size_t size = 1024 * 1024;

  // allocate input and output host buffers
  std::vector<float> in1_h(size);
  std::vector<float> in2_h(size);
  std::vector<float> out_h(size);

  // fill the input buffers with random data, and the output buffer with zeros
  for (size_t i = 0; i < size; ++i) {
    in1_h[i] = dist(rand);
    in2_h[i] = dist(rand);
    out_h[i] = 0.;
  }

  SECTION("Test add_vectors_f") {
    // allocate input and output buffers on the device
    float* in1_d = sycl::malloc_device<float>(size, queue);
    float* in2_d = sycl::malloc_device<float>(size, queue);
    float* out_d = sycl::malloc_device<float>(size, queue);

    // copy the input data to the device
    REQUIRE_NOTHROW(queue.memcpy(in1_d, in1_h.data(), size * sizeof(float)));
    REQUIRE_NOTHROW(queue.memcpy(in2_d, in2_h.data(), size * sizeof(float)));

    // fill the output buffer with zeros
    REQUIRE_NOTHROW(queue.memset(out_d, 0, size * sizeof(float)));

    // launch the 1-dimensional kernel for vector addition
    REQUIRE_NOTHROW(queue.submit([&](sycl::handler& cgh) {
      cgh.parallel_for(sycl::nd_range<1>{32 * 32, 32},
                       [=](sycl::nd_item<1> idx) { kernel_add_vectors_f(in1_d, in2_d, out_d, size); });
    }));

    // copy the results from the device to the host
    REQUIRE_NOTHROW(queue.memcpy(out_h.data(), out_d, size * sizeof(float)));

    // wait for all the operations to complete
    REQUIRE_NOTHROW(queue.wait());

    // check the results
    for (size_t i = 0; i < size; ++i) {
      float sum = in1_h[i] + in2_h[i];
      CHECK_THAT(out_h[i], Catch::Matchers::WithinAbs(sum, epsilon));
    }
  }
}
