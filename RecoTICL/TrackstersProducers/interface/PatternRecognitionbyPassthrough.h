// Author: Felice Pantaleo - felice.pantaleo@cern.ch
// Date: 05/2024

#ifndef __RecoTICL_TrackstersProducers_PatternRecognitionbyPassthrough_H__
#define __RecoTICL_TrackstersProducers_PatternRecognitionbyPassthrough_H__
#include <memory>  // unique_ptr
#include "RecoTICL/TrackstersProducers/interface/PatternRecognitionAlgoBase.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

namespace ticl {
  template <typename TILES>
  class PatternRecognitionbyPassthrough final : public PatternRecognitionAlgoBaseT<TILES> {
  public:
    PatternRecognitionbyPassthrough(const edm::ParameterSet& conf, edm::ConsumesCollector);
    ~PatternRecognitionbyPassthrough() override = default;

    void makeTracksters(const typename PatternRecognitionAlgoBaseT<TILES>::Inputs& input,
                        std::vector<Trackster>& result,
                        std::unordered_map<int, std::vector<int>>& seedToTracksterAssociation) override;

    void filter(std::vector<Trackster>& output,
                const std::vector<Trackster>& inTracksters,
                const typename PatternRecognitionAlgoBaseT<TILES>::Inputs& input,
                std::unordered_map<int, std::vector<int>>& seedToTracksterAssociation) override;

    static void fillPSetDescription(edm::ParameterSetDescription& iDesc);

  private:
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
    hgcal::RecHitTools rhtools_;
  };

}  // namespace ticl

#endif  // __RecoTICL_TrackstersProducers_PatternRecognitionbyPassthrough_H__
