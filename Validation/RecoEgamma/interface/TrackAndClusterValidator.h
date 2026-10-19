#ifndef Validation_EGamma_TrackAndClusterValidator_h
#define Validation_EGamma_TrackAndClusterValidator_h

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"
#include "PhysicsTools/TruthInfo/interface/BranchHitAssociator.h"
#include "PhysicsTools/TruthInfo/interface/RecoHitAdapters.h"
#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/TruthInfo/interface/LogicalGraphHitIndex.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "Validation/RecoEgamma/interface/Common.h"

namespace egamma {

  struct TrackAndClusterValidatorHistograms {
    dqm::reco::MonitorElement *den_SuperCluster_eta, *num_SuperCluster_eta, *num_SuperClusterDuplicate_eta;
    dqm::reco::MonitorElement *den_SuperCluster_ene,  *num_SuperCluster_ene,  *num_SuperClusterDuplicate_ene;
    dqm::reco::MonitorElement *den_SuperClusterFake_eta, *num_SuperClusterFake_eta;
    dqm::reco::MonitorElement *den_SuperClusterFake_ene,  *num_SuperClusterFake_ene;

    dqm::reco::MonitorElement *den_GsfTrack_eta, *num_GsfTrack_eta;
    dqm::reco::MonitorElement *den_GsfTrack_pt,  *num_GsfTrack_pt;
    dqm::reco::MonitorElement *den_GsfTrackFake_eta, *num_GsfTrackFake_eta;
    dqm::reco::MonitorElement *den_GsfTrackFake_pt,  *num_GsfTrackFake_pt;
    dqm::reco::MonitorElement *track_purity, *track_response_pt;
  };

  class TrackAndClusterValidator {
  public:
    typedef dqm::legacy::DQMStore DQMStore;
    typedef dqm::legacy::MonitorElement MonitorElement;

    TrackAndClusterValidator(edm::EDGetTokenT<std::vector<reco::GsfTrack>> tracks,
                  std::vector<edm::EDGetTokenT<reco::SuperClusterCollection>> superClusters,
                  const float mp, const float ptmin, const float maxeta, const float shEne, const float score, const float trackScore);
    ~TrackAndClusterValidator();

    using Histograms = TrackAndClusterValidatorHistograms;

    void bookHistos(DQMStore::IBooker& ibook, Histograms& histograms, std::string baseDir) const;

    void fillHistos(const edm::Event& event,
                    const Histograms& histograms,
                    const truth::Graph& graph,
                    const truth::BranchHitAssociator& caloAssociator,
                    const std::vector<unsigned int>& selected,
                    const std::vector<bool>& assignable,
                    const TrackAssociationMap& branchToGsf,
                    const TrackAssociationMap& gsfToBranch
    ) const;

  private:
    edm::EDGetTokenT<std::vector<reco::GsfTrack>> trackToken_;
    std::vector<edm::EDGetTokenT<reco::SuperClusterCollection>> superclusterTokens_;
    const float minPurity_;
    const float electronPtMin_;
    const float electronMaxAbsEta_;
    const float shEnergyCut_;
    const float clScoreCut_;
    const float trackScoreCut_;
  };
}
#endif


