#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"

#include "DataFormats/TICLReco/interface/TICLCandidate.h"
typedef SimpleCollectionFlatTableProducer<TICLCandidate> TICLCandidateTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TICLCandidateTableProducer);
