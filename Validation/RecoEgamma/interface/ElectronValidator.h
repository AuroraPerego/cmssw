#ifndef Validation_EGamma_ElectronValidator_h
#define Validation_EGamma_ElectronValidator_h

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"
#include "PhysicsTools/TruthInfo/interface/BranchHitAssociator.h"
#include "PhysicsTools/TruthInfo/interface/RecoHitAdapters.h"
#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/TruthInfo/interface/LogicalGraphHitIndex.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "Validation/RecoEgamma/interface/Common.h"

namespace egamma {

  struct ElectronValidatorHistograms {
    std::vector<dqm::reco::MonitorElement*> num_ele_eta;
    std::vector<dqm::reco::MonitorElement*> num_ele_ene;
    std::vector<dqm::reco::MonitorElement*> num_eleDuplicate_eta;
    std::vector<dqm::reco::MonitorElement*> num_eleDuplicate_ene;
    std::vector<dqm::reco::MonitorElement*> den_ele_eta;
    std::vector<dqm::reco::MonitorElement*> den_ele_ene;
    std::vector<dqm::reco::MonitorElement*> num_ele_fake_eta;
    std::vector<dqm::reco::MonitorElement*> den_ele_fake_eta;
    std::vector<dqm::reco::MonitorElement*> num_ele_fake_ene;
    std::vector<dqm::reco::MonitorElement*> den_ele_fake_ene;
    std::vector<dqm::reco::MonitorElement*> num_track_pt;
    std::vector<dqm::reco::MonitorElement*> num_track_eta;
    std::vector<dqm::reco::MonitorElement*> den_track_pt;
    std::vector<dqm::reco::MonitorElement*> den_track_eta;
    std::vector<dqm::reco::MonitorElement*> num_track_fake_pt;
    std::vector<dqm::reco::MonitorElement*> num_track_fake_eta;
    std::vector<dqm::reco::MonitorElement*> den_track_fake_pt;
    std::vector<dqm::reco::MonitorElement*> den_track_fake_eta;
    std::vector<dqm::reco::MonitorElement*> ele_mismatch_eta;
    std::vector<dqm::reco::MonitorElement*> ele_mismatch_ene;
  };

  class ElectronValidator {
  public:
    typedef dqm::legacy::DQMStore DQMStore;
    typedef dqm::legacy::MonitorElement MonitorElement;

    ElectronValidator(std::vector<edm::EDGetTokenT<std::vector<reco::GsfElectronCore>>> electronCores,
                  std::vector<edm::EDGetTokenT<std::vector<reco::GsfElectron>>> electrons,
                  std::vector<std::string> labels, const float mp,
                  const float ptmin, const float maxeta, const float shEne, const float score, const float trackScore);
    ~ElectronValidator();

    using Histograms = ElectronValidatorHistograms;

    void bookHistos(DQMStore::IBooker& ibook, Histograms& histograms, std::string baseDir) const;

    void fillHistos(const edm::Event& event,
                    const Histograms& histograms,
                    const truth::Graph& graph,
                    const truth::BranchHitAssociator& caloAssociator,
                    const std::vector<unsigned int>& selected,
                    const std::vector<bool>& assignable,
                    const TrackAssociationMap& branchToGsf,
                    const TrackAssociationMap& gsfToBranch,
                    const reco::GsfTrack* first
    ) const;

  private:
    std::vector<edm::EDGetTokenT<std::vector<reco::GsfElectronCore>>> electronCoreTokens_;
    std::vector<edm::EDGetTokenT<std::vector<reco::GsfElectron>>> electronTokens_;
    const uint16_t nCollections_;
    std::vector<std::string> labels_;
    const float minPurity_;
    const float electronPtMin_;
    const float electronMaxAbsEta_;
    const float shEnergyCut_;
    const float clScoreCut_;
    const float trackScoreCut_;
  };
}
#endif

