#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
//#include "SimCalorimetry/HGCalAssociatorProducers/interface/AssociatorTools.h"
#include "Validation/RecoEgamma/interface/EGammaValidator.h"

EGammaValidator::EGammaValidator(const edm::ParameterSet& cfg)
  : graphToken_(consumes(cfg.getParameter<edm::InputTag>("graph"))),
    hitIndexToken_(consumes(cfg.getParameter<edm::InputTag>("hitIndex"))),
    selectedRootsToken_(consumes(cfg.getParameter<edm::InputTag>("selectedRoots"))),
    assignableRootsToken_(consumes(cfg.getParameter<edm::InputTag>("assignableRoots"))),
    gsfTracksToken_(consumes(cfg.getParameter<edm::InputTag>("gsfTracks"))),
    ecalSeedsToken_(consumes<std::vector<reco::ElectronSeed>>(cfg.getParameter<edm::InputTag>("ecalBasedSeeds"))),
    trackerSeedsToken_(consumes<std::vector<reco::ElectronSeed>>(cfg.getParameter<edm::InputTag>("tkBasedSeeds"))),
    mergedSeedsToken_(consumes<std::vector<reco::ElectronSeed>>(cfg.getParameter<edm::InputTag>("mergedSeeds"))),
    minPurity_(cfg.getParameter<float>("minPurity")),
    electronPtMin_(cfg.getParameter<float>("electronPtMin")),
    electronMaxAbsEta_(cfg.getParameter<float>("electronMaxAbsEta")),
    shEnergyCut_(cfg.getParameter<float>("shEnergyCut")),
    clScoreCut_(cfg.getParameter<float>("clScoreCut")),
    trackScoreCut_(cfg.getParameter<float>("trackScoreCut")),
    dirName_(cfg.getParameter<std::string>("dirName")) {

  std::vector<std::string> labels;
  for (auto const &tag : cfg.getParameter<std::vector<edm::InputTag>>("electronCores")) {
    electronCoresToken_.emplace_back(consumes<std::vector<reco::GsfElectronCore>>(tag));
    labels.emplace_back(tag.label());
  }

  for (auto const &tag : cfg.getParameter<std::vector<edm::InputTag>>("electrons")) {
    electronsToken_.emplace_back(consumes<std::vector<reco::GsfElectron>>(tag));
    labels.emplace_back(tag.label());
  }

  for (auto const &tag : cfg.getParameter<std::vector<edm::InputTag>>("superClusters")) {
    superClustersToken_.emplace_back(consumes<reco::SuperClusterCollection>(tag));
  }

  seedValidator_ = std::make_unique<egamma::SeedValidator>(ecalSeedsToken_, trackerSeedsToken_, mergedSeedsToken_, electronPtMin_, electronMaxAbsEta_);
  eleValidator_ = std::make_unique<egamma::ElectronValidator>(electronCoresToken_, electronsToken_, labels, minPurity_, electronPtMin_, electronMaxAbsEta_, shEnergyCut_, clScoreCut_, trackScoreCut_);
  trkAndClusValidator_ = std::make_unique<egamma::TrackAndClusterValidator>(gsfTracksToken_, superClustersToken_, minPurity_, electronPtMin_, electronMaxAbsEta_, shEnergyCut_, clScoreCut_, trackScoreCut_);

}

EGammaValidator::~EGammaValidator() {}

void EGammaValidator::bookHistograms(DQMStore::IBooker& ibook,
                                    edm::Run const&,
                                    edm::EventSetup const& setup,
                                    Histograms& histograms) const {
  seedValidator_->bookHistos(ibook, histograms.seedHistos, dirName_+"electronSeeds");
  eleValidator_->bookHistos(ibook, histograms.eleHistos, dirName_+"gsfElectrons");
  trkAndClusValidator_->bookHistos(ibook, histograms.trackAndClusHistos, dirName_+"tracksAndClusters");
}


void EGammaValidator::dqmAnalyze(const edm::Event& event,
                                const edm::EventSetup& setup,
                                const Histograms& histograms) const {
  auto const& graph = event.get(graphToken_);
  auto const& selected = event.get(selectedRootsToken_);
  std::vector<bool> assignable(graph.nParticles(), false);
  for (auto id : event.get(assignableRootsToken_))
    if (id < assignable.size())
      assignable[id] = true;

  const auto & hitIndex = event.get(hitIndexToken_);
  // One associator per event, on the tracker channel, counting shared hits.
  const truth::BranchHitAssociator trackAssociator(hitIndex,
                                              selected,
                                              truth::BranchHitAssociator::Metric::SharedHits,
                                              truth::HitChannel::Tracker,
                                              /*emptyRootsMeansAll=*/false);

  const truth::BranchHitAssociator caloAssociator(hitIndex,
                                                  selected,
                                                  truth::BranchHitAssociator::Metric::SharedEnergy,
                                                  truth::HitChannel::Calo,
                                                  /*emptyRootsMeansAll=*/false);

  std::vector<uint32_t> rootRow(graph.nParticles(), truth::BranchMatch::kInvalidRoot);
  for (std::size_t i = 0; i < selected.size(); ++i)
    rootRow[selected[i]] = static_cast<uint32_t>(i);

  // reco
  auto const& gsfTracks = event.get(gsfTracksToken_);

  TrackAssociationMap branchToGsf(static_cast<unsigned int>(selected.size()));
  TrackAssociationMap gsfToBranch(static_cast<unsigned int>(gsfTracks.size()));
  for (std::size_t t = 0; t < gsfTracks.size(); ++t) {
    const auto hits = truth::recoHits(gsfTracks[t]);
    for (auto const& match : trackAssociator.bestBranches(std::span<const truth::RecoHit>(hits))) {
      if (match.rootParticleId >= assignable.size() || !assignable[match.rootParticleId])
        continue;
      const uint32_t row = rootRow[match.rootParticleId];
      branchToGsf.insert(row, t, match.reverseScore);
      gsfToBranch.insert(t, row, match.score);
    }
  }
  branchToGsf.sort(false);
  gsfToBranch.sort(false);

  std::vector<edm::Handle<reco::SuperClusterCollection>> superClusters_h(superClustersToken_.size());
  edm::MultiSpan<reco::SuperCluster> superClusters;
  for (unsigned int i = 0; i < superClustersToken_.size(); ++i) {
    event.getByToken(superClustersToken_[i], superClusters_h[i]);
    //Fill MultiSpan
    superClusters.add(*superClusters_h[i]);
  }

  seedValidator_->fillHistos(event, histograms.seedHistos, graph, caloAssociator, trackAssociator, assignable);
  eleValidator_->fillHistos(event, histograms.eleHistos, graph, caloAssociator, selected, assignable, branchToGsf, gsfToBranch, &gsfTracks[0]);
  trkAndClusValidator_->fillHistos(event, histograms.trackAndClusHistos, graph, caloAssociator, selected, assignable, branchToGsf, gsfToBranch);
}


void EGammaValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // hgcalValidator
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("graph", edm::InputTag("truthLogicalGraphProducer"));
  desc.add<edm::InputTag>("hitIndex", edm::InputTag("truthLogicalGraphHitIndexProducer"));
  desc.add<edm::InputTag>("selectedRoots", edm::InputTag("truthBranchTargets", "selectedRoots"))
      ->setComment("the candidate particles");
  desc.add<edm::InputTag>("assignableRoots", edm::InputTag("truthBranchTargets", "assignableRoots"))
      ->setComment("the candidates a reco object may be assigned to: no parton, boson or beam particle");
  desc.add<edm::InputTag>("gsfTracks", edm::InputTag("electronGsfTracks"));
  desc.add<std::vector<edm::InputTag>>("electronCores", {edm::InputTag("ecalDrivenGsfElectronCores")});
  desc.add<std::vector<edm::InputTag>>("electrons", {edm::InputTag("ecalDrivenGsfElectrons")});
  desc.add<std::vector<edm::InputTag>>("superClusters", {edm::InputTag("particleFlowSuperClusterECAL", "particleFlowSuperClusterECALBarrel")});
  desc.add<edm::InputTag>("ecalBasedSeeds", edm::InputTag("ecalDrivenElectronSeeds"));
  desc.add<edm::InputTag>("tkBasedSeeds", edm::InputTag("trackerDrivenElectronSeeds","SeedsForGsf"));
  desc.add<edm::InputTag>("mergedSeeds", edm::InputTag("electronMergedSeeds"));
  desc.add<float>("minPurity", 0.75);
  desc.add<float>("electronPtMin", 5.);
  desc.add<float>("electronMaxAbsEta", 3.);
  desc.add<float>("shEnergyCut", 0.5);
  desc.add<float>("clScoreCut", 0.5);
  desc.add<float>("trackScoreCut", 0.9);
  desc.add<std::string>("dirName", "Egamma/truthGraph/");
  descriptions.add("egammaValidator", desc);
}

DEFINE_FWK_MODULE(EGammaValidator);
