#ifndef DataFormats_EgammaReco_SuperClusterDeviceCollection_H
#define DataFormats_EgammaReco_SuperClusterDeviceCollection_H

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/EgammaReco/interface/SuperClusterSoA.h"
#include <alpaka/alpaka.hpp>

namespace ALPAKA_ACCELERATOR_NAMESPACE::reco {

  using SuperClusterDeviceCollection = PortableCollection<::reco::SuperClusterSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::reco

#endif
