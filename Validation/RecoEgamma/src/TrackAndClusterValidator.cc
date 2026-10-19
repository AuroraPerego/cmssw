#include "Validation/RecoEgamma/interface/TrackAndClusterValidator.h"

namespace egamma {

  TrackAndClusterValidator::TrackAndClusterValidator(
      edm::EDGetTokenT<std::vector<reco::GsfTrack>> tracks,
      std::vector<edm::EDGetTokenT<reco::SuperClusterCollection>> superClusters,
      const float mp, const float ptmin, const float maxeta, const float shEne, const float score, const float trackScore)
      : trackToken_(tracks), superclusterTokens_(superClusters), minPurity_(mp), electronPtMin_(ptmin), electronMaxAbsEta_(maxeta), shEnergyCut_(shEne), clScoreCut_(score), trackScoreCut_(trackScore) {};

  TrackAndClusterValidator::~TrackAndClusterValidator() {};

  void TrackAndClusterValidator::bookHistos(DQMStore::IBooker& ibook,
                                            Histograms& histograms,
                                            std::string baseDir) const {
    ibook.cd();
    ibook.setCurrentFolder(baseDir);

    auto bookEta = [&](const std::string& numDem) {
      return ibook.book1D(numDem + "_vs_eta", numDem + " eta;#eta (GeV)", 50, -3, 3);
    };
    auto bookEnePt = [&](const std::string& numDem, const std::string& var) {
      const std::string xtitle = (var == "pt") ? "#p_T [GeV]" : "Energy [GeV]";
      const int max = (var == "pt") ? 100 : 500;
      return ibook.book1D(numDem + "_vs_" + var, numDem + " " + var + ";" + xtitle, 50, 1, max);
    };

    histograms.den_SuperCluster_eta = bookEta("den_SuperCluster");
    histograms.num_SuperCluster_eta = bookEta("num_SuperCluster");
    histograms.num_SuperClusterDuplicate_eta = bookEta("num_SuperClusterDup");
    histograms.den_SuperCluster_ene = bookEnePt("den_SuperCluster", "ene");
    histograms.num_SuperCluster_ene = bookEnePt("num_SuperCluster", "ene");
    histograms.num_SuperClusterDuplicate_ene = bookEnePt("num_SuperClusterDup", "ene");

    histograms.den_SuperClusterFake_eta = bookEta("den_SuperClusterFake");
    histograms.num_SuperClusterFake_eta = bookEta("num_SuperClusterFake");
    histograms.den_SuperClusterFake_ene = bookEnePt("den_SuperClusterFake", "ene");
    histograms.num_SuperClusterFake_ene = bookEnePt("num_SuperClusterFake", "ene");

    histograms.den_GsfTrack_eta = bookEta("den_GsfTrack");
    histograms.num_GsfTrack_eta = bookEta("num_GsfTrack");
    histograms.den_GsfTrack_pt = bookEnePt("den_GsfTrack", "pt");
    histograms.num_GsfTrack_pt = bookEnePt("num_GsfTrack", "pt");

    histograms.den_GsfTrackFake_eta = bookEta("den_GsfTrackFake");
    histograms.num_GsfTrackFake_eta = bookEta("num_GsfTrackFake");
    histograms.den_GsfTrackFake_pt = bookEnePt("den_GsfTrackFake", "pt");
    histograms.num_GsfTrackFake_pt = bookEnePt("num_GsfTrackFake", "pt");

    histograms.track_purity = ibook.bookProfile("gsfTrackPurity", "track purity;#p_{T}^{reco};purity", 50, 1, 200, 0, 1, "s");
    histograms.track_response_pt = ibook.bookProfile("gsfTrackResponseVsPt", "gsf track response vs pt;#p_{T}^{sim} [GeV];#p_{T}^{reco}/p_{T}^{sim}", 50, 1., 200., 0., 2., "s");
  };

  void TrackAndClusterValidator::fillHistos(const edm::Event& event,
                                            const Histograms& histograms,
                                            const truth::Graph& graph,
                                            const truth::BranchHitAssociator& caloAssociator,
                                            const std::vector<unsigned int>& selected,
                                            const std::vector<bool>& assignable,
                                            const TrackAssociationMap& branchToGsf,
                                            const TrackAssociationMap& gsfToBranch) const {
    auto const& tracks = event.get(trackToken_);
    std::vector<edm::Handle<reco::SuperClusterCollection>> superClusters_h(superclusterTokens_.size());
    edm::MultiSpan<reco::SuperCluster> superClusters;
    for (unsigned int i = 0; i < superclusterTokens_.size(); ++i) {
      event.getByToken(superclusterTokens_[i], superClusters_h[i]);
      superClusters.add(*superClusters_h[i]);
    }

    std::vector<uint32_t> rootRow(graph.nParticles(), truth::BranchMatch::kInvalidRoot);
    for (std::size_t i = 0; i < selected.size(); ++i)
      rootRow[selected[i]] = static_cast<uint32_t>(i);

    ClusterAssociationMap branchToCluster(static_cast<unsigned int>(selected.size()));
    for (unsigned int c = 0; c < superClusters.size(); ++c) {
      const auto hits = truth::recoHits(superClusters[c]);
      if (hits.empty())
        continue;
      bool fake = true;
      for (auto const& m : caloAssociator.bestBranches(std::span<const truth::RecoHit>(hits))) {
        const uint32_t row = rootRow[m.rootParticleId];
        if (row == truth::BranchMatch::kInvalidRoot || m.rootParticleId >= assignable.size() || !assignable[m.rootParticleId])
          continue;
        branchToCluster.insert(row, c, m.sharedEnergy, m.reverseScore);
        if (m.sharedEnergy / superClusters[c].energy() > shEnergyCut_)
          fake = false;
      }
      const float scEta = superClusters[c].eta();
      const float scEne = superClusters[c].energy();
      histograms.den_SuperClusterFake_eta->Fill(scEta);
      histograms.den_SuperClusterFake_ene->Fill(scEne);
      if (fake) {
        histograms.num_SuperClusterFake_eta->Fill(scEta);
        histograms.num_SuperClusterFake_ene->Fill(scEne);
      }
    }
    branchToCluster.sort(truth::byAscendingScore);

    for (std::size_t t = 0; t < tracks.size(); ++t) {
      float bestPurity = 0.f;
      int bestRoot = -1;
      for (const auto& match : gsfToBranch[t]) {
        const uint32_t root = selected[match.index()];  // riga -> id della radice
        if (root >= assignable.size() || !assignable[root])
          continue;
        const float purity = 1.f - static_cast<float>(match.fraction());
        if (purity > bestPurity) {
          bestPurity = purity;
          bestRoot = static_cast<int>(root);
        }
      }

      const float trkEta = tracks[t].eta();
      const float trkPt = tracks[t].pt();
      histograms.den_GsfTrackFake_eta->Fill(trkEta);
      histograms.den_GsfTrackFake_pt->Fill(trkPt);

      if (bestRoot >= 0 && bestPurity >= minPurity_ && std::abs(truth::Particle(&graph, bestRoot).pdgId()) == 11) {
        histograms.track_purity->Fill(trkPt, bestPurity);
        const truth::Particle electron(&graph, bestRoot);
        histograms.track_response_pt->Fill(electron.momentum().pt(), tracks[t].pt() / electron.momentum().pt());
      } else {
        histograms.num_GsfTrackFake_eta->Fill(trkEta);
        histograms.num_GsfTrackFake_pt->Fill(trkPt);
      }
    }

    // ---------------- truth loop: efficiency and duplicates ----------------
    for (uint32_t id = 0; id < graph.nParticles(); ++id) {
      const truth::Particle electron(&graph, id);
      if (std::abs(electron.pdgId()) != 11 || !electron.hasGen() || electron.status() != 1)
        continue;
      auto const& p4 = electron.momentum();
      const float eta = p4.eta();
      if (p4.pt() < electronPtMin_ || std::abs(eta) > electronMaxAbsEta_)
        continue;
      const uint32_t row = rootRow[id];
      if (row == truth::BranchMatch::kInvalidRoot)
        continue;

      unsigned nCl = 0;
      for (const auto& cl : branchToCluster[row])
        if (cl.score() < clScoreCut_)
          ++nCl;

      histograms.den_SuperCluster_eta->Fill(eta);
      histograms.den_SuperCluster_ene->Fill(p4.P());
      if (nCl >= 1) {
        histograms.num_SuperCluster_eta->Fill(eta);
        histograms.num_SuperCluster_ene->Fill(p4.P());
      }
      if (nCl > 1) {
        histograms.num_SuperClusterDuplicate_eta->Fill(eta);
        histograms.num_SuperClusterDuplicate_ene->Fill(p4.P());
      }

      bool eff = false;
      for (const auto& tr : branchToGsf[row])
        if (tr.fraction() < trackScoreCut_)
          eff = true;

      histograms.den_GsfTrack_eta->Fill(eta);
      histograms.den_GsfTrack_pt->Fill(p4.pt());
      if (eff) {
        histograms.num_GsfTrack_eta->Fill(eta);
        histograms.num_GsfTrack_pt->Fill(p4.pt());
      }
    }
  };
}  // namespace egamma

