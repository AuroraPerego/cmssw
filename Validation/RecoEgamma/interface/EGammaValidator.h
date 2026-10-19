#ifndef Validation_EGamma_EGammaValidator_h
#define Validation_EGamma_EGammaValidator_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"

#include "PhysicsTools/TruthInfo/interface/BranchHitAssociator.h"
#include "PhysicsTools/TruthInfo/interface/RecoHitAdapters.h"
#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/TruthInfo/interface/LogicalGraphHitIndex.h"
#include "SimDataFormats/TruthInfo/interface/Particle.h"

#include "SimDataFormats/Associations/interface/TICLAssociationMap.h"

#include "DQMServices/Core/interface/DQMGlobalEDAnalyzer.h"

#include "Validation/RecoEgamma/interface/SeedValidator.h"
#include "Validation/RecoEgamma/interface/ElectronValidator.h"
#include "Validation/RecoEgamma/interface/TrackAndClusterValidator.h"

#include <vector>

struct EGammaHistograms {
  egamma::SeedValidatorHistograms seedHistos;
  egamma::ElectronValidatorHistograms eleHistos;
  egamma::TrackAndClusterValidatorHistograms trackAndClusHistos;
};

class EGammaValidator : public DQMGlobalEDAnalyzer<EGammaHistograms> {
public:
  using Histograms = EGammaHistograms;

  /// Constructor
  EGammaValidator(const edm::ParameterSet& pset);

  /// Destructor
  ~EGammaValidator() override;

  /// Method called once per event
  void dqmAnalyze(const edm::Event&, const edm::EventSetup&, const Histograms&) const override;
  /// Method called to book the DQM histograms
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&, Histograms&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  const edm::EDGetTokenT<truth::Graph> graphToken_;
  const edm::EDGetTokenT<truth::LogicalGraphHitIndex> hitIndexToken_;
  const edm::EDGetTokenT<std::vector<unsigned int>> selectedRootsToken_;
  const edm::EDGetTokenT<std::vector<unsigned int>> assignableRootsToken_;
  const edm::EDGetTokenT<reco::GsfTrackCollection> gsfTracksToken_;
  std::vector<edm::EDGetTokenT<reco::SuperClusterCollection>> superClustersToken_;
  std::vector<edm::EDGetTokenT<std::vector<reco::GsfElectronCore>>> electronCoresToken_;
  std::vector<edm::EDGetTokenT<std::vector<reco::GsfElectron>>> electronsToken_;
  const edm::EDGetTokenT<std::vector<reco::ElectronSeed>> ecalSeedsToken_;
  const edm::EDGetTokenT<std::vector<reco::ElectronSeed>> trackerSeedsToken_;
  const edm::EDGetTokenT<std::vector<reco::ElectronSeed>> mergedSeedsToken_;
  std::unique_ptr<egamma::SeedValidator> seedValidator_;
  std::unique_ptr<egamma::ElectronValidator> eleValidator_;
  std::unique_ptr<egamma::TrackAndClusterValidator> trkAndClusValidator_;
  const float minPurity_, electronPtMin_, electronMaxAbsEta_, shEnergyCut_, clScoreCut_, trackScoreCut_;
  std::string dirName_;
};

#endif

