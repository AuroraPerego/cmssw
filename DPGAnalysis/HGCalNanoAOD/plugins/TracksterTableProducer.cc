#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"

#include "DataFormats/TICL/interface/Trackster.h"
typedef SimpleCollectionFlatTableProducer<ticl::Trackster> TracksterCollectionTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TracksterCollectionTableProducer);
