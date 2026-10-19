#ifndef Validation_EGamma_SeedValidator_h
#define Validation_EGamma_SeedValidator_h

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"
#include "PhysicsTools/TruthInfo/interface/BranchHitAssociator.h"
#include "PhysicsTools/TruthInfo/interface/RecoHitAdapters.h"
#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/TruthInfo/interface/LogicalGraphHitIndex.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "Validation/RecoEgamma/interface/Common.h"

namespace egamma {
  struct SeedValidatorHistograms {
    dqm::reco::MonitorElement *den_EcalSeed_eta, *num_EcalSeed_eta, *num_EcalDuplicateSeed_eta;
    dqm::reco::MonitorElement *den_EcalSeed_pt, *num_EcalSeed_pt, *num_EcalDuplicateSeed_pt;

    dqm::reco::MonitorElement *den_TrkSeed_eta, *num_TrkSeed_eta, *num_TrkDuplicateSeed_eta;
    dqm::reco::MonitorElement *den_TrkSeed_pt, *num_TrkSeed_pt, *num_TrkDuplicateSeed_pt;

    dqm::reco::MonitorElement *den_BothSeed_eta, *num_BothSeed_eta;
    dqm::reco::MonitorElement *den_BothSeed_pt, *num_BothSeed_pt;

    dqm::reco::MonitorElement *den_MergedSeed_eta, *num_MergedSeed_eta, *num_MergedDuplicateSeed_eta;
    dqm::reco::MonitorElement *den_MergedSeed_pt, *num_MergedSeed_pt, *num_MergedDuplicateSeed_pt;
  };

  class SeedValidator {
  public:
    typedef dqm::legacy::DQMStore DQMStore;
    typedef dqm::legacy::MonitorElement MonitorElement;

    SeedValidator(edm::EDGetTokenT<std::vector<reco::ElectronSeed>> ecalSeeds,
                  edm::EDGetTokenT<std::vector<reco::ElectronSeed>> trackerSeeds,
                  edm::EDGetTokenT<std::vector<reco::ElectronSeed>> mergeSeeds,
                  const float ptmin, const float maxeta);
    ~SeedValidator();

    using Histograms = SeedValidatorHistograms;

    void bookHistos(DQMStore::IBooker &ibook, Histograms &histograms, std::string baseDir) const;

    void fillHistos(const edm::Event &event,
                    const Histograms &histograms,
                    const truth::Graph& graph,
                    const truth::BranchHitAssociator& caloAssociator,
                    const truth::BranchHitAssociator& trackAssociator,
                    const std::vector<bool>& assignable
    ) const;

  private:
    edm::EDGetTokenT<std::vector<reco::ElectronSeed>> ecalSeedToken_;
    edm::EDGetTokenT<std::vector<reco::ElectronSeed>> trackerSeedToken_;
    edm::EDGetTokenT<std::vector<reco::ElectronSeed>> mergeSeedToken_;
    const float electronPtMin_;
    const float electronMaxAbsEta_;
  };
}  // namespace egamma
#endif

