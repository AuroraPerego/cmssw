#include "Validation/RecoEgamma/interface/SeedValidator.h"

namespace egamma {

  SeedValidator::SeedValidator(edm::EDGetTokenT<std::vector<reco::ElectronSeed>> ecalSeeds,
                               edm::EDGetTokenT<std::vector<reco::ElectronSeed>> trackerSeeds,
                               edm::EDGetTokenT<std::vector<reco::ElectronSeed>> mergeSeeds,
                               const float ptmin, const float maxeta)
      : ecalSeedToken_(ecalSeeds), trackerSeedToken_(trackerSeeds), mergeSeedToken_(mergeSeeds), electronPtMin_(ptmin), electronMaxAbsEta_(maxeta) {}

  SeedValidator::~SeedValidator() {};

  void SeedValidator::bookHistos(DQMStore::IBooker& ibook, Histograms& histograms, std::string baseDir) const {
    ibook.cd();
    ibook.setCurrentFolder(baseDir);

    auto bookEta = [&](const std::string& numDem) {
      return ibook.book1D(numDem + "_vs_eta", numDem + " eta;#eta (GeV)", 50, -3, 3);
    };

    auto bookPt = [&](const std::string& numDem) {
      return ibook.book1D(numDem + "_vs_pt", numDem + " pt;#p_T (GeV)", 50, 1, 100);
    };

    histograms.den_EcalSeed_eta = bookEta("denEcalSeed");
    histograms.num_EcalSeed_eta = bookEta("numEcalSeed");
    histograms.num_EcalDuplicateSeed_eta = bookEta("numEcalSeedDup");
    histograms.den_EcalSeed_pt = bookPt("denEcalSeed");
    histograms.num_EcalSeed_pt = bookPt("numEcalSeed");
    histograms.num_EcalDuplicateSeed_pt = bookPt("numEcalSeedDup");

    histograms.den_TrkSeed_eta = bookEta("denTrkSeed");
    histograms.num_TrkSeed_eta = bookEta("numTrkSeed");
    histograms.num_TrkDuplicateSeed_eta = bookEta("numTrkSeedDup");
    histograms.den_TrkSeed_pt = bookPt("denTrkSeed");
    histograms.num_TrkSeed_pt = bookPt("numTrkSeed");
    histograms.num_TrkDuplicateSeed_pt = bookPt("numTrkSeedDup");

    histograms.den_BothSeed_eta = bookEta("denBothSeed");
    histograms.num_BothSeed_eta = bookEta("numBothSeed");
    histograms.den_BothSeed_pt = bookPt("denBothSeed");
    histograms.num_BothSeed_pt = bookPt("numBothSeed");

    histograms.den_MergedSeed_eta = bookEta("denMergedSeed");
    histograms.num_MergedSeed_eta = bookEta("numMergedSeed");
    histograms.num_MergedDuplicateSeed_eta = bookEta("numMergedSeedDup");
    histograms.den_MergedSeed_pt = bookPt("denMergedSeed");
    histograms.num_MergedSeed_pt = bookPt("numMergedSeed");
    histograms.num_MergedDuplicateSeed_pt = bookPt("numMergedSeedDup");
  };

  void SeedValidator::fillHistos(const edm::Event& event,
                                 const Histograms& histograms,
                                 const truth::Graph& graph,
                                 const truth::BranchHitAssociator& caloAssociator,
                                 const truth::BranchHitAssociator& trackAssociator,
                                 const std::vector<bool>& assignable
  ) const {

    auto const& ecalSeeds = event.get(ecalSeedToken_);
    auto const& trackerSeeds = event.get(trackerSeedToken_);
    auto const& mergedSeeds = event.get(mergeSeedToken_);

    constexpr uint32_t kNone = truth::BranchMatch::kInvalidRoot;

    // Best assignable root for a set of hits, or kNone if nothing passes the score cut.
    auto bestRoot = [&](truth::BranchHitAssociator const& assoc,
                        std::vector<truth::RecoHit> const& hits) -> uint32_t {
      if (hits.empty()) return kNone;
      uint32_t bestId = kNone;
      float bestScore = 1.;
      for (auto const& m : assoc.bestBranches(std::span<const truth::RecoHit>(hits))) {
        if (m.rootParticleId >= assignable.size() || !assignable[m.rootParticleId]) continue;
        if (m.score < bestScore) { bestScore = m.score; bestId = m.rootParticleId; }
      }
      return bestId;
    };

    auto caloRoot = [&](reco::ElectronSeed const& s) -> uint32_t {
      if (s.caloCluster().isNull()) return kNone;
      return bestRoot(caloAssociator, truth::recoHits(*s.caloCluster()));
    };
    auto trackRoot = [&](reco::ElectronSeed const& s) -> uint32_t {
      if (s.ctfTrack().isNull()) return kNone;
      return bestRoot(trackAssociator, truth::recoHits(*s.ctfTrack()));
    };

    // Number of seeds per truth root, one vector per collection.
    std::vector<unsigned> nEcal(graph.nParticles(), 0), nTrk(graph.nParticles(), 0),
        nMerged(graph.nParticles(), 0);

    for (auto const& s : ecalSeeds)
      if (auto r = caloRoot(s); r != kNone) ++nEcal[r];

    for (auto const& s : trackerSeeds)
      if (auto r = trackRoot(s); r != kNone) ++nTrk[r];

    auto mismatch = 0;
    for (auto const& s : mergedSeeds) {
      const uint32_t rc = caloRoot(s);
      const uint32_t rt = trackRoot(s);
      uint32_t r = kNone;
      if (rc != kNone && rt != kNone) {
        if (rc == rt) {
          r = rc;          // ecal and tracker parts agree on the same electron
        } else {
          ++mismatch;
        }
      } else {
        r = (rc != kNone) ? rc : rt;   // only one driver present
      }
      if (r != kNone) ++nMerged[r];
    }
    //mergedMismatch_->Fill(mismatch); TODO fix this plot

    auto fillSeed = [](dqm::reco::MonitorElement* den,
                       dqm::reco::MonitorElement* num,
                       dqm::reco::MonitorElement* numDup,
                       unsigned n, double x) {
      den->Fill(x);
      if (n >= 1) num->Fill(x);
      if (numDup && n > 1) numDup->Fill(x);
    };

    for (uint32_t id = 0; id < graph.nParticles(); ++id) {
      const truth::Particle electron(&graph, id);
      if (!assignable[id] || electron.id() != electron.lastCopy().id() ||
          std::abs(electron.pdgId()) != 11 || !electron.hasGen() || electron.status() != 1)
        continue;
      auto const& p4 = electron.momentum();
      if (p4.pt() < electronPtMin_ || std::abs(p4.eta()) > electronMaxAbsEta_)
        continue;
      const double eta = p4.eta();
      const double pt = p4.pt();

      const unsigned nBoth = (nEcal[id] >= 1 && nTrk[id] >= 1) ? 1 : 0;

      fillSeed(histograms.den_EcalSeed_eta, histograms.num_EcalSeed_eta,
               histograms.num_EcalDuplicateSeed_eta, nEcal[id], eta);
      fillSeed(histograms.den_EcalSeed_pt, histograms.num_EcalSeed_pt,
               histograms.num_EcalDuplicateSeed_pt, nEcal[id], pt);

      fillSeed(histograms.den_TrkSeed_eta, histograms.num_TrkSeed_eta,
               histograms.num_TrkDuplicateSeed_eta, nTrk[id], eta);
      fillSeed(histograms.den_TrkSeed_pt, histograms.num_TrkSeed_pt,
               histograms.num_TrkDuplicateSeed_pt, nTrk[id], pt);

      fillSeed(histograms.den_BothSeed_eta, histograms.num_BothSeed_eta, nullptr, nBoth, eta);
      fillSeed(histograms.den_BothSeed_pt, histograms.num_BothSeed_pt, nullptr, nBoth, pt);

      fillSeed(histograms.den_MergedSeed_eta, histograms.num_MergedSeed_eta,
               histograms.num_MergedDuplicateSeed_eta, nMerged[id], eta);
      fillSeed(histograms.den_MergedSeed_pt, histograms.num_MergedSeed_pt,
               histograms.num_MergedDuplicateSeed_pt, nMerged[id], pt);

      //nDuplicates_->Fill(nEcal[id]);   // unchanged: still a standalone histogram
    }
  };
}

