#ifndef RecoHGCal_TICL_MCFwithNNInterpretationAlgo_H
#define RecoHGCal_TICL_MCFwithNNInterpretationAlgo_H

#include "RecoHGCal/TICL/interface/TICLInterpretationAlgoBase.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
//#include "GeometryHGCalCommonData/interface/HGCalDDDConstants.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include <memory>
#include <vector>

namespace ticl {

  class MCFwithNNInterpretationAlgo : public TICLInterpretationAlgoBase<reco::Track> {
  public:
    MCFwithNNInterpretationAlgo(const edm::ParameterSet &conf, TICLONNXGlobalCache const *cache);
    ~MCFwithNNInterpretationAlgo() override = default;

    void initialize(const HGCalDDDConstants *hgcons,
                    const hgcal::RecHitTools rhtools,
                    const edm::ESHandle<MagneticField> bfieldH,
                    const edm::ESHandle<Propagator> propH) override;

    void makeCandidates(const Inputs &input,
                        edm::Handle<MtdHostCollection> inputTimingh,
                        std::vector<Trackster> &resultTracksters,
                        std::vector<int> &resultCandidate) override;

    static void fillPSetDescription(edm::ParameterSetDescription &desc);

  private:
    // Configurable parameters (the "x" array from Python)
    float drCut_;
    float tsTsScoreShift_;
    float trackTsScoreShift_;
    float tsTsScoreWeight_;
    float trackTsScoreWeight_;
    int neutralPenalty_;
    int tracksterInit_;
    int trackInit_;

    const cms::Ort::ONNXRuntime *onnxSessionTracks_;
    const cms::Ort::ONNXRuntime *onnxSessionTracksters_;

    const HGCalDDDConstants *hgcons_ = nullptr;
    hgcal::RecHitTools rhtools_;
    edm::ESHandle<MagneticField> bfield_;
    edm::ESHandle<Propagator> propagator_;

    std::unique_ptr<GeomDet> firstDisk_[2];

    void buildLayers();

    Vector propagateTrackster(const Trackster &t,
                              const unsigned idx,
                              float zVal,
                              std::array<TICLLayerTile, 2> &tracksterTiles);

    void findTrackstersInWindow(const edm::MultiSpan<Trackster> &tracksters,
                                const std::vector<std::pair<Vector, unsigned>> &seedingCollection,
                                const std::array<TICLLayerTile, 2> &tracksterTiles,
                                const std::vector<Vector> &tracksterPropPoints,
                                float delta,
                                unsigned trackstersSize,
                                std::vector<std::vector<unsigned>> &resultCollection,
                                bool useMask);

    // Scoring helpers
    float normTracks(float x, float y) const;
    float normTracksters(float x, float y) const;
    float computeScore(
        float refPt, float refEta, float refPhi, float refP, float tsEta, float tsPhi, float tsEnergy) const;
    float deltaR2(float eta1, float phi1, float eta2, float phi2) const;
  };

}  // namespace ticl

#endif
