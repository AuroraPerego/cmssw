#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

#include <sycl/sycl.hpp>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SYCLTestDeviceAdditionAlgo.h"

class SYCLTestDeviceAdditionModule : public edm::global::EDAnalyzer<> {
public:
  explicit SYCLTestDeviceAdditionModule(edm::ParameterSet const& config);
  ~SYCLTestDeviceAdditionModule() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void analyze(edm::StreamID, edm::Event const& event, edm::EventSetup const& setup) const override;

private:
  const uint32_t size_;
};

SYCLTestDeviceAdditionModule::SYCLTestDeviceAdditionModule(edm::ParameterSet const& config)
    : size_(config.getParameter<uint32_t>("size_")) {}

void SYCLTestDeviceAdditionModule::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<uint32_t>("size_", 1024 * 1024);
  descriptions.addWithDefaultLabel(desc);
}

void SYCLTestDeviceAdditionModule::analyze(edm::StreamID, edm::Event const& event, edm::EventSetup const& setup) const {
  sycl::queue queue{sycl::cpu_selector_v, sycl::property::queue::in_order()};

  // random number generator with a gaussian distribution
  std::random_device rd{};
  std::default_random_engine rand{rd()};
  std::normal_distribution<float> dist{0., 1.};

  // tolerance
  constexpr float epsilon = 0.000001;

  // allocate input and output host buffers
  std::vector<float> in1_h(size_);
  std::vector<float> in2_h(size_);
  std::vector<float> out_h(size_);

  // fill the input buffers with random data, and the output buffer with zeros
  for (size_t i = 0; i < size_; ++i) {
    in1_h[i] = dist(rand);
    in2_h[i] = dist(rand);
    out_h[i] = 0.;
  }

  // allocate input and output buffers on the device
  float* in1_d = sycl::malloc_device<float>(size_, queue);
  float* in2_d = sycl::malloc_device<float>(size_, queue);
  float* out_d = sycl::malloc_device<float>(size_, queue);

  // copy the input data to the device
  queue.memcpy(in1_d, in1_h.data(), size_ * sizeof(float));
  queue.memcpy(in2_d, in2_h.data(), size_ * sizeof(float));

  // fill the output buffer with zeros
  queue.memset(out_d, 0, size_ * sizeof(float));

  // launch the 1-dimensional kernel for vector addition
  HeterogeneousTestSYCLDevicePlugins::wrapper_add_vectors_f(in1_d, in2_d, out_d, size_, queue);

  // copy the results from the device to the host
  queue.memcpy(out_h.data(), out_d, size_ * sizeof(float));

  // wait for all the operations to complete
  queue.wait();

  // check the results
  for (size_t i = 0; i < size_; ++i) {
    float sum = in1_h[i] + in2_h[i];
    assert(out_h[i] < sum + epsilon);
    assert(out_h[i] > sum - epsilon);
  }

  std::cout << "All tests passed.\n";
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SYCLTestDeviceAdditionModule);
