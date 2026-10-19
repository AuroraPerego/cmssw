#ifndef Validation_EGamma_Common_h
#define Validation_EGamma_Common_h

namespace {
  using ClusterAssociationMap = ticl::TICLAssociationMap<ticl::mapWithSharedEnergyAndScore>;
  using TrackAssociationMap = ticl::TICLAssociationMap<ticl::mapWithFraction>; // only score for tracks
}

#endif
