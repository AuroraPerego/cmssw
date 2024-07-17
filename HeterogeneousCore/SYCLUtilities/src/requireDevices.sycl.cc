#include <cstdlib>
#include <iostream>

#include <sycl/sycl.hpp>

#include "HeterogeneousCore/SYCLUtilities/interface/requireDevices.h"

namespace cms::sycltest {
  bool testDevices() {
    std::vector<sycl::device> devices = sycl::device::get_devices();
    if (devices.size() == 0) {
      std::cerr << "No SYCL devices available, the test will be skipped."
                << "\n";
      return false;
    }
    return true;
  }

  void requireDevices() {
    if (not testDevices()) {
      exit(EXIT_SUCCESS);
    }
  }
}  // namespace cms::sycltest
