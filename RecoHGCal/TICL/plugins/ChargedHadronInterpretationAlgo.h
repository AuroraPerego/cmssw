#ifndef RecoHGCal_TICL_ChargedHadronInterpretationAlgo_h
#define RecoHGCal_TICL_ChargedHadronInterpretationAlgo_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "RecoHGCal/TICL/interface/TICLInterpretationAlgoBase.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"

namespace ticl {

  class ChargedHadronInterpretationAlgo : public TICLInterpretationAlgoBase<reco::Track> {
  public:
    ChargedHadronInterpretationAlgo(const edm::ParameterSet &conf, TICLONNXGlobalCache const* cache);

    ~ChargedHadronInterpretationAlgo() override;

    void makeCandidates(const Inputs &input,
                        edm::Handle<MtdHostCollection> inputTiming_h,
                        std::vector<Trackster> &resultTracksters,
                        std::vector<int> &resultCandidate,
                        std::vector<bool> &maskedTracksters) override;

    // Arbitration mode: one charged-hadron hypothesis per track-linked (merged)
    // trackster, reusing the geometric association of makeCandidates. Neutral
    // leftovers are not emitted; the producer derives neutrals from unclaimed
    // tracksters after arbitration.
    void makeOpinions(const Inputs &input,
                      edm::Handle<MtdHostCollection> inputTiming_h,
                      std::vector<Trackster> &hypothesisTracksters,
                      std::vector<Hypothesis> &hypotheses) override;

    void initialize(const HGCalDDDConstants *hgcons,
                    const hgcal::RecHitTools rhtools,
                    const edm::ESHandle<MagneticField> bfieldH,
                    const edm::ESHandle<Propagator> propH) override;

    static void fillPSetDescription(edm::ParameterSetDescription &iDesc);

  private:
    void buildLayers();

    // MCF+NN parameters
    double drCut_;
    double tsTsScoreShift_;
    double trackTsScoreShift_;
    double tsTsScoreWeight_;
    double trackTsScoreWeight_;

    cms::Ort::ONNXRuntime const* onnxSessionTracks_ = nullptr;
    cms::Ort::ONNXRuntime const* onnxSessionTracksters_ = nullptr;
    const std::vector<std::string> inputNames_;
    const std::vector<std::string> outputNames_;

    const HGCalDDDConstants *hgcons_;

    std::unique_ptr<GeomDet> firstDisk_[2];

    hgcal::RecHitTools rhtools_;

    edm::ESHandle<MagneticField> bfield_;
    edm::ESHandle<Propagator> propagator_;
  };

}  // namespace ticl

#endif
