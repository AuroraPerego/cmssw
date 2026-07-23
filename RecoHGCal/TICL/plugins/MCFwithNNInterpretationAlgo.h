#ifndef RecoHGCal_TICL_MCFwithNNInterpretationAlgo_H
#define RecoHGCal_TICL_MCFwithNNInterpretationAlgo_H

#include "RecoHGCal/TICL/interface/TICLInterpretationAlgoBase.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace ticl {

  class MCFwithNNInterpretationAlgo : public TICLInterpretationAlgoBase<reco::Track> {
  public:
    explicit MCFwithNNInterpretationAlgo(const edm::ParameterSet& conf, TICLONNXGlobalCache const* cache);

    void initialize(const HGCalDDDConstants* hgcons,
                    const hgcal::RecHitTools rhtools,
                    const edm::ESHandle<MagneticField> bfieldH,
                    const edm::ESHandle<Propagator> propH) override;

    void makeCandidates(const Inputs& input,
                        edm::Handle<MtdHostCollection> inputTimingh,
                        std::vector<Trackster>& resultTracksters,
                        std::vector<int>& resultCandidate,
                        std::vector<bool>& maskedTracksters) override;

    static void fillPSetDescription(edm::ParameterSetDescription& desc);

  private:
    float drCut_;
    float tsTsScoreShift_;
    float trackTsScoreShift_;
    float tsTsScoreWeight_;
    float trackTsScoreWeight_;
    int neutralPenalty_;
    int tracksterInit_;
    int trackInit_;
    int manyPenalty_;

    cms::Ort::ONNXRuntime const* onnxSessionTracks_ = nullptr;
    cms::Ort::ONNXRuntime const* onnxSessionTracksters_ = nullptr;
    const std::vector<std::string> inputNames_;
    const std::vector<std::string> outputNames_;

    const HGCalDDDConstants* hgcons_ = nullptr;
    hgcal::RecHitTools rhtools_;
    edm::ESHandle<MagneticField> bfield_;
    edm::ESHandle<Propagator> propagator_;

    std::unique_ptr<GeomDet> firstDisk_[2];

    void buildLayers();
  };

}  // namespace ticl

#endif
