#include "Validation/RecoEgamma/interface/ElectronValidator.h"

namespace egamma {

  ElectronValidator::ElectronValidator(std::vector<edm::EDGetTokenT<std::vector<reco::GsfElectronCore>>> electronCores,
                  std::vector<edm::EDGetTokenT<std::vector<reco::GsfElectron>>> electrons,
                  std::vector<std::string> labels, const float mp,
      const float ptmin, const float maxeta, const float shEne, const float score, const float trackScore)
             : electronCoreTokens_(electronCores),
               electronTokens_(electrons),
               nCollections_(electronCores.size() + electrons.size()),
               labels_(labels), minPurity_(mp), electronPtMin_(ptmin), electronMaxAbsEta_(maxeta), shEnergyCut_(shEne), clScoreCut_(score), trackScoreCut_(trackScore) {}

  ElectronValidator::~ElectronValidator() {};

  void ElectronValidator::bookHistos(DQMStore::IBooker& ibook, Histograms& histograms, std::string baseDir) const {
    ibook.cd();
    ibook.setCurrentFolder(baseDir);

    // loop and fill up to nCollections
    histograms.num_ele_eta.reserve(nCollections_);
    histograms.num_ele_ene.reserve(nCollections_);
    histograms.num_eleDuplicate_eta.reserve(nCollections_);
    histograms.num_eleDuplicate_ene.reserve(nCollections_);
    histograms.den_ele_eta.reserve(nCollections_);
    histograms.den_ele_ene.reserve(nCollections_);
    histograms.num_ele_fake_eta.reserve(nCollections_);
    histograms.den_ele_fake_eta.reserve(nCollections_);
    histograms.num_ele_fake_ene.reserve(nCollections_);
    histograms.den_ele_fake_ene.reserve(nCollections_);
    histograms.num_track_pt.reserve(nCollections_);
    histograms.num_track_eta.reserve(nCollections_);
    histograms.den_track_pt.reserve(nCollections_);
    histograms.den_track_eta.reserve(nCollections_);
    histograms.num_track_fake_pt.reserve(nCollections_);
    histograms.num_track_fake_eta.reserve(nCollections_);
    histograms.den_track_fake_pt.reserve(nCollections_);
    histograms.den_track_fake_eta.reserve(nCollections_);
    histograms.ele_mismatch_eta.reserve(nCollections_);
    histograms.ele_mismatch_ene.reserve(nCollections_);

    auto bookEta = [&](std::vector<dqm::reco::MonitorElement*>& histo, uint16_t i, const std::string& numDem) {
      histo.emplace_back(ibook.book1D(numDem + "_vs_eta_" + labels_[i], labels_[i] + " eta;#eta (GeV)", 50, -3, 3));
    };

    auto bookEnePt = [&](std::vector<dqm::reco::MonitorElement*>& histo, uint16_t i, const std::string& numDem, const std::string& var) {
      const std::string xtitle = (var=="pt") ? "#p_T [GeV]" : "Energy [GeV]";
      const int max =  (var=="pt") ? 100 : 500;
      histo.emplace_back(ibook.book1D(numDem + "_vs_" + var + "_" + labels_[i], labels_[i] + " " + var + ";" + xtitle, 50, 1, max));
    };

    for(auto i = 0; i < nCollections_; ++i) {
      bookEta(histograms.num_ele_eta, i, "num_ele");
      bookEta(histograms.num_eleDuplicate_eta, i, "numDuplicate_ele");
      bookEta(histograms.den_ele_eta, i, "den_ele");
      bookEta(histograms.num_ele_fake_eta, i, "numFake_ele");
      bookEta(histograms.den_ele_fake_eta, i, "denFake_ele");
      bookEta(histograms.num_track_eta, i, "num_track");
      bookEta(histograms.den_track_eta, i, "den_track");
      bookEta(histograms.num_track_fake_eta, i, "numFake_track");
      bookEta(histograms.den_track_fake_eta, i, "denFake_track");
      bookEta(histograms.ele_mismatch_eta, i, "mismatchSCvsGsf");

      bookEnePt(histograms.num_ele_ene, i, "num_ele", "energy");
      bookEnePt(histograms.num_eleDuplicate_ene, i, "numDuplicate_ele", "energy");
      bookEnePt(histograms.den_ele_ene, i, "den_ele", "energy");
      bookEnePt(histograms.num_ele_fake_ene, i, "numFake_ele", "energy");
      bookEnePt(histograms.den_ele_fake_ene, i, "denFake_ele", "energy");
      bookEnePt(histograms.num_track_pt, i, "num_track", "pt");
      bookEnePt(histograms.den_track_pt, i, "den_track", "pt");
      bookEnePt(histograms.num_track_fake_pt, i, "numFake_track", "pt");
      bookEnePt(histograms.den_track_fake_pt, i, "denFake_track", "pt");
      bookEnePt(histograms.ele_mismatch_ene, i, "mismatchSCvsGsf", "energy");
    }
}

  void ElectronValidator::fillHistos(const edm::Event& event,
                    const Histograms& histograms,
                    const truth::Graph& graph,
                    const truth::BranchHitAssociator& caloAssociator,
                    const std::vector<unsigned int>& selected,
                    const std::vector<bool>& assignable,
                    const TrackAssociationMap& branchToGsf,
                    const TrackAssociationMap& gsfToBranch,
                    const reco::GsfTrack* first
    ) const {

    constexpr uint32_t kInvalid = truth::BranchMatch::kInvalidRoot;
    std::vector<uint32_t> rootRow(graph.nParticles(), kInvalid);
    for (std::size_t i = 0; i < selected.size(); ++i)
      rootRow[selected[i]] = static_cast<uint32_t>(i);

    for (unsigned int i = 0; i < nCollections_; ++i) {
      const int offset = i - electronCoreTokens_.size();
      auto process = [&](auto const& electrons) {
        ClusterAssociationMap branchToCluster(static_cast<unsigned int>(selected.size()));
        for (unsigned int c = 0; c < electrons.size(); ++c) {
          const auto& scRef = electrons[c].superCluster();
          const auto hits = truth::recoHits(*scRef);
          if (hits.empty())
            continue;
          auto fake = true;
          uint32_t bestRootSc = kInvalid;
          float bestScore = 1;
          for (auto const& m : caloAssociator.bestBranches(std::span<const truth::RecoHit>(hits))) {
            if (m.rootParticleId >= assignable.size() || !assignable[m.rootParticleId])
              continue;
            const uint32_t row = rootRow[m.rootParticleId]; // row in selected
            branchToCluster.insert(row, c, m.sharedEnergy, m.reverseScore);
            if (m.score < clScoreCut_) fake = false;
            if (m.score < bestScore) {
              bestScore = m.score;
              bestRootSc = m.rootParticleId;
            }
          }
          histograms.den_ele_fake_eta[i]->Fill(scRef->eta());
          histograms.den_ele_fake_ene[i]->Fill(scRef->energy());
          if (fake) {
            histograms.num_ele_fake_eta[i]->Fill(scRef->eta());
            histograms.num_ele_fake_ene[i]->Fill(scRef->energy());
          }

          const auto& gsfRef = electrons[c].gsfTrack();
          const int gsfIdx = &(*gsfRef) - first;
          histograms.den_track_fake_eta[i]->Fill(gsfRef->eta());
          histograms.den_track_fake_pt[i]->Fill(gsfRef->pt());
          bool isFake = true;
          uint32_t bestRootGsf = kInvalid;
          bestScore = 1;
          for (auto const& m : gsfToBranch[gsfIdx]) {
            if ((1.f - m.fraction()) > minPurity_)
              isFake = false;
            if (m.fraction() < bestScore) {
              bestScore = m.fraction();
              bestRootGsf = selected[m.index()];
            }
          }
          if (isFake) {
            histograms.num_track_fake_eta[i]->Fill(gsfRef->eta());
            histograms.num_track_fake_pt[i]->Fill(gsfRef->pt());
          }
          // add mismatch plot where the best clusterToBranch and gsfToBranch point to different particles
          if (bestRootSc != kInvalid && bestRootGsf != kInvalid) {
            if (not (bestRootSc == bestRootGsf)) {
              histograms.ele_mismatch_eta[i]->Fill(scRef->eta());
              histograms.ele_mismatch_ene[i]->Fill(scRef->energy());
            }
          }
        }
        branchToCluster.sort(truth::byAscendingScore);

        for (uint32_t id = 0; id < graph.nParticles(); ++id) {
          const truth::Particle electron(&graph, id);
          // not assignable or not last copy
          if (!assignable[id] || electron.id() != electron.lastCopy().id() || std::abs(electron.pdgId()) != 11 || !electron.hasGen() || electron.status() != 1)
            continue;
          auto const& p4 = electron.momentum();
          if (p4.pt() < electronPtMin_ || std::abs(p4.eta()) > electronMaxAbsEta_)
            continue;

          // Let's use the branchToCluster to get the supercluster matched to the same truth electron
          auto nCl = 0;
          const auto& branchClusters = branchToCluster[rootRow[id]];
          for (const auto& cl : branchClusters) {
            if (cl.score() < clScoreCut_)
              nCl++;
          }
          const auto eta = p4.eta();
          const auto ene = p4.P();
          histograms.den_ele_eta[i]->Fill(eta);
          histograms.den_ele_ene[i]->Fill(ene);

          if (nCl >= 1) {
            histograms.num_ele_eta[i]->Fill(eta);
            histograms.num_ele_ene[i]->Fill(ene);
          }
          if (nCl > 1) {
            histograms.num_eleDuplicate_eta[i]->Fill(eta);
            histograms.num_eleDuplicate_ene[i]->Fill(ene);
          }

          // now check the track
          if (branchClusters.size() == 0) continue; // no reco electrons for this sim
          const auto pt = p4.pt();
          histograms.den_track_eta[i]->Fill(eta);
          histograms.den_track_pt[i]->Fill(pt);

          const auto& gsfRef = electrons[branchClusters[0].index()].gsfTrack(); // can loop over all above threshold
          const unsigned int gsfIdx = &(*gsfRef) - first;
          for (auto const& match : branchToGsf[rootRow[electron.id()]]) {
            if (match.index() == gsfIdx and match.fraction() < trackScoreCut_) {
              histograms.num_track_eta[i]->Fill(eta);
              histograms.num_track_pt[i]->Fill(pt);
              break;
            }
          }
        }
      };

      if (offset < 0)
        process(event.get(electronCoreTokens_[i]));
      else
        process(event.get(electronTokens_[offset]));
    }
  }
}

