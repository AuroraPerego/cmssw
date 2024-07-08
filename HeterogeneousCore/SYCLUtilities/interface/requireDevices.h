#ifndef HeterogeneousCore_SYCLUtilities_requireDevices_h
#define HeterogeneousCore_SYCLUtilities_requireDevices_h

/**
 * These functions are meant to be called only from unit tests.
 */
namespace cms {
  namespace sycltest {
    /// In presence of SYCL devices, return true; otherwise print message and return false
    bool testDevices();

    /// Print message and exit if there are no SYCL devices
    void requireDevices();
  }  // namespace sycltest
}  // namespace cms

#endif  // HeterogeneousCore_SYCLUtilities_requireDevices_h
