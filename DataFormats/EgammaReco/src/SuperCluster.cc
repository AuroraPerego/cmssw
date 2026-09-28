#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include <algorithm>
using namespace reco;

SuperCluster::SuperCluster(double energy, const math::XYZPoint& position)
    : CaloCluster(energy, position),
      rawEnergy_(0),
      phiWidth_(0),
      etaWidth_(0) {}

SuperCluster::SuperCluster(double energy,
                           const math::XYZPoint& position,
                           const CaloClusterPtr& seed,
                           const CaloClusterPtrVector& clusters,
                           double phiWidth,
                           double etaWidth)
    : CaloCluster(energy, position), rawEnergy_(0) {
  phiWidth_ = phiWidth;
  etaWidth_ = etaWidth;
  seed_ = seed;

  // set references to constituent basic clusters and update list of rechits
  for (CaloClusterPtrVector::const_iterator bcit = clusters.begin(); bcit != clusters.end(); ++bcit) {
    clusters_.push_back((*bcit));

    // updated list of used hits
    const std::vector<std::pair<DetId, float> >& v1 = (*bcit)->hitsAndFractions();
    for (std::vector<std::pair<DetId, float> >::const_iterator diIt = v1.begin(); diIt != v1.end(); ++diIt) {
      hitsAndFractions_.push_back((*diIt));
    }  // loop over rechits
  }  // loop over basic clusters

  computeRawEnergy();
}

void SuperCluster::computeRawEnergy() {
  rawEnergy_ = 0.;
  for (CaloClusterPtrVector::const_iterator bcItr = clustersBegin(); bcItr != clustersEnd(); bcItr++) {
    rawEnergy_ += (*bcItr)->energy();
  }
}
