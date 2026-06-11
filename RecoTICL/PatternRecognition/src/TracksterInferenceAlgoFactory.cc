#include "RecoTICL/PatternRecognition/interface/TracksterInferenceAlgoFactory.h"
#include "RecoTICL/PatternRecognition/interface/TracksterInferenceByPFN.h"
#include "RecoTICL/PatternRecognition/interface/TracksterInferenceByDNN.h"
#include "RecoTICL/PatternRecognition/interface/TracksterInferenceByCNN.h"

#include "FWCore/ParameterSet/interface/ValidatedPluginFactoryMacros.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginMacros.h"

EDM_REGISTER_VALIDATED_PLUGINFACTORY(TracksterInferenceAlgoFactory, "TracksterInferenceAlgoFactory");
DEFINE_EDM_VALIDATED_PLUGIN(TracksterInferenceAlgoFactory, ticl::TracksterInferenceByPFN, "TracksterInferenceByPFN");
DEFINE_EDM_VALIDATED_PLUGIN(TracksterInferenceAlgoFactory, ticl::TracksterInferenceByDNN, "TracksterInferenceByDNN");
DEFINE_EDM_VALIDATED_PLUGIN(TracksterInferenceAlgoFactory, ticl::TracksterInferenceByCNN, "TracksterInferenceByCNN");
