#ifndef DataFormats_EgammaReco_SuperClusterSoA_H
#define DataFormats_EgammaReco_SuperClusterSoA_H

// Authors: Simone Balducci, Felice Pantaleo, Wahid Redjeb, Aurora Perego, Leonardo Beltrame

#include "DataFormats/SoATemplate/interface/SoABlocks.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/TICL/interface/AssociationMap.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterSoA.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include <xtd/xtd.h>

namespace reco {

  template <std::size_t Size, bool Boolean>
  using ClustersLayout = typename ticl::AssociationMapLayout<uint32_t, uint32_t>::template Layout<Size, Boolean>;

  GENERATE_SOA_LAYOUT(SuperClusterSoAExtra,
                      SOA_COLUMN(int, seed),
                      SOA_COLUMN(float, rawEnergy),
                      SOA_COLUMN(float, phiWidth),
                      SOA_COLUMN(float, etaWidth))

  GENERATE_SOA_BLOCKS(SuperClusterSoALayout,
                      SOA_BLOCK(position, CaloClusterSoAPosition),
                      SOA_BLOCK(energy, CaloClusterSoAEnergy),
                      SOA_BLOCK(indexes, CaloClusterSoAIndexes),
                      SOA_BLOCK(timing, CaloClusterSoATiming),
                      SOA_BLOCK(extra, SuperClusterSoAExtra),
                      SOA_BLOCK(clusters, ClustersLayout),
                      SOA_CONST_VIEW_METHODS(
                        SOA_HOST_DEVICE auto r(std::integral auto idx) const {
                          const auto x = this->position()[idx].x();
                          const auto y = this->position()[idx].y();
                          return xtd::sqrt(x * x + y * y);
                        }
                        SOA_HOST_DEVICE auto phi(std::integral auto idx) const {
                          const auto x = this->position()[idx].x();
                          const auto y = this->position()[idx].y();
                          return xtd::atan2(y, x);
                        }
                        SOA_HOST_DEVICE auto eta(std::integral auto idx) const {
                          const auto x = this->position()[idx].x();
                          const auto y = this->position()[idx].y();
                          const auto z = this->position()[idx].z();
                          const auto theta = xtd::atan2(xtd::sqrt(x * x + y * y), z);
                          return -xtd::log(xtd::tan(theta / 2));
                        }
                        // NOTE: this and the next probably don't work on GPU
                        // and they assume that seedID is the ID of the seed caloCluster
                        SOA_HOST_DEVICE auto seedCrysIEtaOrIx(std::integral auto idx) const {
                          auto detid = this->indexes()[idx].seedID();
                          int ietaorix = 0;
                          if (detid.subdetId() == EcalBarrel) {
                            EBDetId ebdetid(detid);
                            ietaorix = ebdetid.ieta();
                          } else if (detid.subdetId() == EcalEndcap) {
                            EEDetId eedetid(detid);
                            ietaorix = eedetid.ix();
                          }
                          return ietaorix;
                        }
                        SOA_HOST_DEVICE auto seedCrysIPhiOrIy(std::integral auto idx) const {
                          auto detid = this->indexes()[idx].seedID();
                          int iphioriy = 0;
                          if (detid.subdetId() == EcalBarrel) {
                            EBDetId ebdetid(detid);
                            iphioriy = ebdetid.iphi();
                          } else if (detid.subdetId() == EcalEndcap) {
                            EEDetId eedetid(detid);
                            iphioriy = eedetid.iy();
                          }
                          return iphioriy;
                        }
                        SOA_HOST_DEVICE auto computeRawEnergy(std::integral auto idx) {
                          float raw_energy = 0.f;
                          for (const auto& cl : this->clusters()[idx])
                            raw_energy += cl; //.energy();  TODO: fix this based on what is saved in the map
                          return raw_energy;
                        }
                      )
  )

  using SuperClusterSoA = SuperClusterSoALayout<>;
  using SuperClusterSoAView = SuperClusterSoA::View;
  using SuperClusterSoAConstView = SuperClusterSoA::ConstView;

}  // namespace reco

#endif
