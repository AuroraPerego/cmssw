#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>

#include <alpaka/alpaka.hpp>
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "DataFormats/HGCalReco/interface/HGCalSoAClusters.h"
#include "DataFormats/HGCalReco/interface/HGCalSoARecHitsHostCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoAClustersDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/HGCalSoARecHitsExtraDeviceCollection.h"
#include "DataFormats/HGCalReco/interface/alpaka/TracksterSoADeviceCollection.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "CLUEstering/CLUEstering.hpp"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // FIXME this file should not contain alpaka stuff, should be a .cc and the alpaka part should go in a .dev.cc
  struct KernelFillTracksters {
    template <typename TAcc, typename MapViewT, typename TsSoAViewT, typename LcViewT>
    ALPAKA_FN_ACC void operator()(
        const TAcc& acc, MapViewT dev_map, TsSoAViewT tsView, LcViewT lcView, int32_t n_tracksters) const {
      // Each work-item handles one trackster index
      for (auto tsIdx : alpaka::uniformElements(acc, n_tracksters)) {
        auto& trackster = tsView[tsIdx];
        auto keys_ptr = dev_map[tsIdx].data();
        int32_t keys_n = dev_map[tsIdx].size();

        // Initialize some SoA fields (example fields: energy, emEnergy, barycenterX/Y/Z)
        float accumEnergy = 0.f;
        float bx = 0.f, by = 0.f, bz = 0.f;

        // iterate the layer cluster keys and gather information
        for (int k = 0; k < keys_n; ++k) {
          const auto lc_idx = keys_ptr[k];  // layer cluster index

          // read LC energy / position from layerClustersView
          float lcE = lcView.energy()[lc_idx];
          float lcX = lcView.x()[lc_idx];
          float lcY = lcView.y()[lc_idx];
          float lcZ = lcView.z()[lc_idx];

          trackster.vertices()[k] = lc_idx;
          trackster.vertex_multiplicity()[k] = 1;
          const float fraction = 1.f / trackster.vertex_multiplicity()[k];

          accumEnergy += lcE * fraction;

          // compute barycenter weighted by lcE
          bx += (lcE * fraction) * lcX;
          by += (lcE * fraction) * lcY;
          bz += (lcE * fraction) * lcZ;
        }  // keys loop

        // normalize barycenter if needed (avoid division by zero)
        if (accumEnergy != 0.f) {
          bx /= accumEnergy;
          by /= accumEnergy;
          bz /= accumEnergy;
        }

        // --- write results into SoA at position tsIdx ---
        trackster.raw_energy() = accumEnergy;
        trackster.barycenter().x() = bx;
        trackster.barycenter().y() = by;
        trackster.barycenter().z() = bz;

        // FIXME
        // trackster.calculateRawPt();
        // trackster.calculateRawEmPt();

        // compute trackster time
        constexpr float c = 29.9792458;  // cm/ns
        float tracksterTime = 0.f;
        int num = 0;
        for (int k = 0; k < keys_n; ++k) {
          const auto time = lcView.time()[keys_ptr[k]];
          if (time > 0.f) {
            // calcolo delta T (assuming Test Beam setup)
            float deltaT = (lcView.z()[keys_ptr[k]] - trackster.barycenter().z()) / c;
            tracksterTime += (time - deltaT);
            num++;
          }
        }
        if (tracksterTime > 0.f) {
          trackster.time() = tracksterTime / num;
          trackster.timeError() = 0.f;
        } else {
          trackster.time() = 0.f;
          trackster.timeError() = -1.f;
        }
      }
    }
  };

  class HeterogeneousTracksterProducer : public stream::EDProducer<> {
  public:
    HeterogeneousTracksterProducer(edm::ParameterSet const& config)
        : EDProducer(config),
          deviceTokenSoAClusters_{consumes(config.getParameter<edm::InputTag>("layerClusters"))},
          legacyTrackstersToken_{produces()},
          rho_(config.getParameter<double>("rho_c")) {
      auto dc_vec = config.getParameter<std::vector<double>>("dc");
      auto dm_vec = config.getParameter<std::vector<double>>("dm");

      if (dc_vec.size() != 3 || dm_vec.size() != 3) {
        throw cms::Exception("Configuration") << "Parameters 'dc' and 'dm' must each have exactly 3 elements.";
      }

      for (size_t i = 0; i < 3; ++i) {
        dc_[i] = static_cast<float>(dc_vec[i]);
        dm_[i] = static_cast<float>(dm_vec[i]);
      }
    }
    ~HeterogeneousTracksterProducer() override = default;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("layerClusters", edm::InputTag("hgcalRecHitsLayerClustersSoA"));
      desc.add<double>("rho_c", 0.6);
      desc.add<std::vector<double>>("dc", {2., 2., 2});
      desc.add<std::vector<double>>("dm", {1.8, 1.8, 2});
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(device::Event& iEvent, device::EventSetup const& iSetup) override {
      const auto& lc = iEvent.get(deviceTokenSoAClusters_);
      auto& queue = iEvent.queue();
      auto x = const_cast<float*>(lc.view().x().data());
      auto y = const_cast<float*>(lc.view().y().data());
      auto z = const_cast<float*>(lc.view().z().data());
      auto E = const_cast<float*>(lc.view().energy().data());

      const int32_t n = static_cast<int32_t>(lc->metadata().size());
      if (n > 0) {
        auto d_clIndex =
            cms::alpakatools::make_device_buffer<int[]>(queue, n);  // temporary buffer needed by CLUEstering
        auto dp_clIndex = const_cast<int*>(d_clIndex.data());
        clue::PointsDevice<3> d_points(queue, n, x, y, z, E, dp_clIndex);

        clue::Clusterer<3> algo(queue, dc_, rho_, dm_);
        algo.make_clusters(queue, d_points);
        // create hosts points and do the copy
        clue::PointsHost<3> h_points(queue, n);
        clue::copyToHost(queue, h_points, d_points);
        //alpaka::wait(queue);

        auto time = cms::alpakatools::make_host_buffer<float[]>(n);
        alpaka::memcpy(queue,
                       cms::alpakatools::make_host_view(h_points.view().coords[0], n),
                       cms::alpakatools::make_device_view(alpaka::getDev(queue), x, n),
                       (unsigned int)n);
        alpaka::memcpy(queue,
                       cms::alpakatools::make_host_view(h_points.view().coords[1], n),
                       cms::alpakatools::make_device_view(alpaka::getDev(queue), y, n),
                       (unsigned int)n);
        alpaka::memcpy(queue,
                       cms::alpakatools::make_host_view(h_points.view().coords[2], n),
                       cms::alpakatools::make_device_view(alpaka::getDev(queue), z, n),
                       (unsigned int)n);
        alpaka::memcpy(queue,
                       cms::alpakatools::make_host_view(h_points.weights(), n),
                       cms::alpakatools::make_device_view(alpaka::getDev(queue), E, n),
                       (unsigned int)n);
        alpaka::memcpy(
            queue,
            cms::alpakatools::make_host_view(alpaka::getPtrNative(time), n),
            cms::alpakatools::make_device_view(alpaka::getDev(queue), const_cast<float*>(lc.view().time().data()), n),
            (unsigned int)n);
        alpaka::wait(queue);

        // compute trackster properties
        bool energyWeight = true;
        auto xHost = h_points.coords(0).data();
        auto yHost = h_points.coords(1).data();
        auto zHost = h_points.coords(2).data();
        auto EHost = h_points.weights();

        // DEBUG PRINT
        std::cout << "Event Number of LCs " << n << std::endl;
        std::unordered_map<float, std::vector<int>> map;

        for (int i = 0; i < lc->metadata().size(); ++i) {
          map[zHost[i]].push_back(i);
        }
        for (const auto& [Z, indices] : map) {
          std::cout << "z = " << Z << " -> Clusters : ";
          for (auto i : indices)
            std::cout << "\t( " << xHost[i] << ", " << yHost[i] << ", " << zHost[i] << ", " << EHost[i] << ", "
                      << time[i] << ")" << std::endl;
          std::cout << std::endl;
        }
        // END DEBUG PRINT

        // BEGIN TRACKSTERS SOA ATTEMPT
        // call get clusters with device points to get the trackster soa size
        const auto devTsMap = clue::get_clusters(queue, d_points);
        const auto devTsMap_v = devTsMap.view();
        std::size_t sizeTs = 0;
        alpaka::memcpy(queue, &sizeTs, devTsMap.size(), sizeof(std::size_t));
        alpaka::wait(queue);
        // do the stuff below only if sizeTs > 0
        TracksterSoADeviceCollection tracksterSoA(sizeTs, queue);
        auto tracksterSoA_v = tracksterSoA.view();
        if (sizeTs != 0) {
          // define work div
          auto work_div = cms::alpakatools::make_workdiv<Acc1D>(static_cast<uint32_t>(sizeTs), 256);
          alpaka::exec<Acc1D>(queue,
                              work_div,
                              KernelFillTracksters{},
                              devTsMap_v,
                              tracksterSoA_v,
                              lc.view(),
                              static_cast<int32_t>(sizeTs));
        }
        // END TRACKSTERS SOA ATTEMPT

        // get LCs indices in tracksters and fill the trackster collection
        const auto tsMap = clue::get_clusters(h_points);
        auto tracksters = std::vector<ticl::Trackster>(tsMap.size());
        std::cout << "Event Number of Tracksters " << tsMap.size() << std::endl;

        for (long unsigned int i = 0; i < tsMap.size(); ++i) {
          auto& trackster = tracksters[i];

          const auto [beginLC, endLC] = tsMap.equal_range(i);
          std::copy(beginLC, endLC, std::back_inserter(trackster.vertices()));
          tracksters[i].vertex_multiplicity().resize(trackster.vertices().size(), 1);

          size_t N = trackster.vertices().size();
          if (N == 0)  // useless?
            continue;

          Eigen::Vector3f point;
          point << 0., 0., 0.;
          Eigen::Vector3f barycenter;
          barycenter << 0., 0., 0.;

          auto fillPoint = [&](const float x, const float y, const float z, const float weight = 1.f) {
            point[0] = weight * x;
            point[1] = weight * y;
            point[2] = weight * z;
          };

          // Initialize this trackster with default, dummy values
          trackster.setRawEnergy(0.f);
          trackster.setRawEmEnergy(0.f);
          trackster.setRawPt(0.f);
          trackster.setRawEmPt(0.f);

          float weight = 1.f / N;

          std::vector<float> layerClusterEnergies;

          for (size_t i = 0; i < N; ++i) {
            auto lcIdx = trackster.vertices(i);
            auto fraction = 1.f / trackster.vertex_multiplicity(i);
            trackster.addToRawEnergy(EHost[lcIdx] * fraction);
            // trackster.addToRawEmEnergy(EHost[lcIdx] * fraction);

            // Compute the weighted barycenter.
            if (energyWeight)
              weight = EHost[lcIdx] * fraction;
            fillPoint(xHost[lcIdx], yHost[lcIdx], zHost[lcIdx], weight);
            for (size_t j = 0; j < 3; ++j)
              barycenter[j] += point[j];

            layerClusterEnergies.push_back(EHost[lcIdx]);
          }
          float raw_energy = trackster.raw_energy();
          float inv_raw_energy = 1.f / raw_energy;
          if (energyWeight)
            barycenter *= inv_raw_energy;
          trackster.setBarycenter(ticl::Trackster::Vector(barycenter));

          trackster.calculateRawPt();
          trackster.calculateRawEmPt();

          // compute trackster time
          constexpr float c = 29.9792458;  // cm/ns
          float tracksterTime = 0.f;
          int num = 0;
          for (size_t i = 0; i < N; ++i) {
            if (time[i] > 0.f) {
              // calcolo delta T (assuming Test Beam setup)
              float deltaT = (zHost[i] - trackster.barycenter().z()) / c;
              tracksterTime += (time[i] - deltaT);
              num++;
            }
          }
          if (tracksterTime > 0.f)
            trackster.setTimeAndError(tracksterTime / num, 0.f);
          else
            trackster.setTimeAndError(-99.f, -1.f);

          std::cout << "  LC in TS: ";
          for (const auto& lc : trackster.vertices())
            std::cout << lc << " ";
          std::cout << std::endl;
          std::cout << "  energy raw: " << trackster.raw_energy() << std::endl;
          std::cout << "  barycenter: " << trackster.barycenter().x() << ", " << trackster.barycenter().y() << ", "
                    << trackster.barycenter().z() << std::endl;
          std::cout << "  time: " << trackster.time() << std::endl;
        }

        iEvent.emplace(legacyTrackstersToken_, std::move(tracksters));
      } else {
        auto tracksters = std::vector<ticl::Trackster>();
        iEvent.emplace(legacyTrackstersToken_, std::move(tracksters));
      }
    }

  private:
    device::EDGetToken<HGCalSoAClustersDeviceCollection> const deviceTokenSoAClusters_;
    edm::EDPutTokenT<std::vector<ticl::Trackster>> const legacyTrackstersToken_;
    float rho_;
    std::array<float, 3> dc_;
    std::array<float, 3> dm_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(HeterogeneousTracksterProducer);
