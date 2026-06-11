// Author: Felice Pantaleo - felice.pantaleo@cern.ch
// Date: 07/2024

#ifndef RecoTICL_PatternRecognition_TracksterInferenceAlgoFactory_H__
#define RecoTICL_PatternRecognition_TracksterInferenceAlgoFactory_H__

#include "FWCore/PluginManager/interface/PluginFactory.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "RecoTICL/PatternRecognition/interface/TracksterInferenceAlgoBase.h"

typedef edmplugin::PluginFactory<ticl::TracksterInferenceAlgoBase*(const edm::ParameterSet&,
                                                                   ticl::TICLONNXGlobalCache const*)>
    TracksterInferenceAlgoFactory;

#endif
