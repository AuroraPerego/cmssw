#ifndef DataFormats_EgammaReco_SuperClusterHostCollection_H
#define DataFormats_EgammaReco_SuperClusterHostCollection_H

#include "DataFormats/EgammaReco/interface/SuperClusterSoA.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include <alpaka/alpaka.hpp>

namespace reco {

  using SuperClusterHostCollection = PortableHostCollection<SuperClusterSoA>;

}  // namespace reco

#endif
