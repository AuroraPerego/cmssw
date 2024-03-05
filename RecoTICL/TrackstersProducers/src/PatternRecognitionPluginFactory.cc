#include "RecoTICL/TrackstersProducers/interface/PatternRecognitionPluginFactory.h"
#include "RecoTICL/TrackstersProducers/interface/PatternRecognitionbyCA.h"
#include "RecoTICL/TrackstersProducers/interface/PatternRecognitionbyCLUE3D.h"
#include "RecoTICL/TrackstersProducers/interface/PatternRecognitionbyFastJet.h"
#include "RecoTICL/TrackstersProducers/interface/PatternRecognitionbyPassthrough.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginFactoryMacros.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginMacros.h"

EDM_REGISTER_VALIDATED_PLUGINFACTORY(PatternRecognitionFactory, "PatternRecognitionFactory");
EDM_REGISTER_VALIDATED_PLUGINFACTORY(PatternRecognitionHFNoseFactory, "PatternRecognitionHFNoseFactory");
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionFactory, ticl::PatternRecognitionbyCA<TICLLayerTiles>, "CA");
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionFactory, ticl::PatternRecognitionbyCLUE3D<TICLLayerTiles>, "CLUE3D");
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionFactory, ticl::PatternRecognitionbyFastJet<TICLLayerTiles>, "FastJet");
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionFactory,
                            ticl::PatternRecognitionbyPassthrough<TICLLayerTiles>,
                            "Passthrough");
DEFINE_EDM_VALIDATED_PLUGIN(PatternRecognitionHFNoseFactory, ticl::PatternRecognitionbyCA<TICLLayerTilesHFNose>, "CA");
