#include "MCFInterpretationAlgo.h"
#include "RecoHGCal/TICL/interface/MinCostFlow.h"
#include "RecoHGCal/TICL/plugins/TICLInterpretationPluginFactory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <unordered_map>

using namespace ticl;
using Vector = ticl::Trackster::Vector;

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
MCFInterpretationAlgo::MCFInterpretationAlgo(const edm::ParameterSet& conf,
                                             edm::ConsumesCollector cc)
    : TICLInterpretationAlgoBase<reco::Track>(conf, cc),
      drCut_(conf.getParameter<double>("drCut")),
      tsTsScoreShift_(conf.getParameter<double>("tsTsScoreShift")),
      trackTsScoreShift_(conf.getParameter<double>("trackTsScoreShift")),
      tsTsScoreWeight_(conf.getParameter<double>("tsTsScoreWeight")),
      trackTsScoreWeight_(conf.getParameter<double>("trackTsScoreWeight")),
      neutralPenalty_(conf.getParameter<double>("neutralPenalty")) {}

// ---------------------------------------------------------------------------
// initialize + buildLayers
// ---------------------------------------------------------------------------
void MCFInterpretationAlgo::initialize(const HGCalDDDConstants* hgcons,
                                       const hgcal::RecHitTools rhtools,
                                       const edm::ESHandle<MagneticField> bfieldH,
                                       const edm::ESHandle<Propagator> propH) {
  hgcons_ = hgcons;
  rhtools_ = rhtools;
  bfield_ = bfieldH;
  propagator_ = propH;
  buildLayers();
}

void MCFInterpretationAlgo::buildLayers() {
  float zVal = hgcons_->waferZ(1, true);
  std::pair<float, float> rMinMax = hgcons_->rangeR(zVal, true);
  for (int iSide = 0; iSide < 2; ++iSide) {
    float zSide = (iSide == 0) ? (-1.f * zVal) : zVal;
    firstDisk_[iSide] = std::make_unique<GeomDet>(
        Disk::build(Disk::PositionType(0, 0, zSide),
                    Disk::RotationType(),
                    SimpleDiskBounds(rMinMax.first, rMinMax.second, zSide - 0.5f, zSide + 0.5f))
            .get());
  }
}

// ---------------------------------------------------------------------------
// Scoring helpers
// ---------------------------------------------------------------------------
float MCFInterpretationAlgo::normTracks(float x, float y) const {
  y = std::abs(y);
  constexpr float a=2.29516386e-01f, b=-8.05371532e+02f, c=-6.45573586e-01f;
  constexpr float d=8.05370082e+02f, e=1.00000033e+00f,  f=-1.07458042e-04f;
  constexpr float g=4.79631755e-01f, h=1.21330763e+00f;
  if (x > 200.f) x = 200.f;
  if (y < 1.7f)  y = 1.7f;
  return std::max(a + b*x + c*y + d*std::pow(x,e) + f*x*y + g*std::pow(y,h), 0.008f);
}

float MCFInterpretationAlgo::normTracksters(float x, float y) const {
  y = std::abs(y);
  constexpr float a=1.49662496e+00f, b=4.04830636e+01f,  c=-1.19840947e+01f;
  constexpr float d=-4.04834456e+01f, e=9.99984896e-01f, f=-1.67329993e-03f;
  constexpr float g=1.06828890e+01f,  h=1.09862164e+00f;
  if (x > 200.f) x = 200.f;
  return a + b*x + c*y + d*std::pow(x,e) + g*std::pow(y,h) + f*x*y;
}

float MCFInterpretationAlgo::computeScore(float refPt, float refEta, float refPhi, float refP,
                                          float tsEta,  float tsPhi,  float tsEnergy) const {
  constexpr float wEp = 0.5f;
  float rEpNorm = (tsEnergy - refP) / refP;
  float dphi = reco::deltaPhi(refPhi, tsPhi);
  float deta = refEta - tsEta;
  float pullR2 = (deta*deta + dphi*dphi) / std::pow(normTracks(refPt, refEta), 2);
  return std::sqrt(wEp * rEpNorm*rEpNorm + (1.f - wEp) * pullR2);
}

// ---------------------------------------------------------------------------
// makeCandidates
// ---------------------------------------------------------------------------
void MCFInterpretationAlgo::makeCandidates(const Inputs& input,
                                           edm::Handle<MtdHostCollection> inputTimingh,
                                           std::vector<Trackster>& resultTracksters,
                                           std::vector<int>& resultCandidate) {
  const auto& tracksters = input.tracksters;
  const auto  tkH        = input.tracksHandle;
  const auto& tracks     = *tkH;
  const auto& maskTracks = input.maskedTracks;
  const float drCut2     = drCut_ * drCut_;

  auto bFieldProd = bfield_.product();
  const Propagator& prop = *propagator_;

  // -----------------------------------------------------------------------
  // 1. Propagate all tracksters to HGCal front face, fill tiles
  // -----------------------------------------------------------------------
  std::array<TICLLayerTile, 2> tracksterPropTiles = {};
  std::vector<Vector> tsAllProp;
  tsAllProp.reserve(tracksters.size());

  const float zVal_layer1 = hgcons_->waferZ(1, true);

  for (unsigned i = 0; i < tracksters.size(); ++i) {
    const Vector& baryc   = tracksters[i].barycenter();
    Vector        directnv = baryc.unit();
    float zVal = zVal_layer1 * (baryc.Z() > 0 ? 1.f : -1.f);
    float par  = (zVal - baryc.Z()) / directnv.Z();
    Vector tPoint(par * directnv.X() + baryc.X(),
                  par * directnv.Y() + baryc.Y(),
                  zVal);
    if      (tPoint.Eta() > 0) tracksterPropTiles[1].fill(tPoint.Eta(), tPoint.Phi(), i);
    else if (tPoint.Eta() < 0) tracksterPropTiles[0].fill(tPoint.Eta(), tPoint.Phi(), i);
    tsAllProp.emplace_back(tPoint);
  }

  // -----------------------------------------------------------------------
  // 2. Propagate valid tracks to HGCal front face
  // -----------------------------------------------------------------------
  struct TrackInfo { int origIdx; float eta, phi; double pt, p; };
  std::vector<TrackInfo> validTracks;
  validTracks.reserve(tracks.size());

  for (unsigned i = 0; i < tracks.size(); ++i) {
    if (!maskTracks[i]) continue;
    const auto& tk = tracks[i];
    if (tk.pt() < 1.f || tk.p() < 2.f) continue;
    if (std::abs(tk.eta()) < 1.5f || std::abs(tk.eta()) > 3.0f) continue;

    int iSide = int(tk.eta() > 0);
    FreeTrajectoryState fts = tk.outerOk()
        ? trajectoryStateTransform::outerFreeState(tk, bFieldProd)
        : trajectoryStateTransform::initialFreeState(tk, bFieldProd);

    const auto& tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
    if (!tsos.isValid()) continue;

    GlobalPoint pp = tsos.globalPosition();
    validTracks.push_back({static_cast<int>(i), pp.eta(), pp.phi(), tk.pt(), tk.p()});
  }
  const int nTracks = static_cast<int>(validTracks.size());

  // -----------------------------------------------------------------------
  // 3. findNeighbours lambda (tile-based, mirrors findTrackstersInWindow)
  // -----------------------------------------------------------------------
  auto findNeighbours = [&](float seed_eta, float seed_phi) -> std::vector<unsigned> {
    bool sideZ = seed_eta > 0;
    const TICLLayerTile& tile = tracksterPropTiles[sideZ];
    float eta_min = std::max(std::fabs(seed_eta) - drCut_, (float)TileConstants::minEta);
    float eta_max = std::min(std::fabs(seed_eta) + drCut_, (float)TileConstants::maxEta);
    auto search_box = tile.searchBoxEtaPhi(eta_min, eta_max, seed_phi - drCut_, seed_phi + drCut_);

    std::vector<unsigned> result;
    for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
      for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
        const auto& in_tile = tile[tile.globalBin(eta_i, phi_i % TileConstants::nPhiBins)];
        for (unsigned t_i : in_tile) {
          float deta = tsAllProp[t_i].Eta() - seed_eta;
          float dphi = reco::deltaPhi((float)tsAllProp[t_i].Phi(), seed_phi);
          if (deta*deta + dphi*dphi < drCut2)
            result.push_back(t_i);
        }
      }
    }
    return result;
  };

  // -----------------------------------------------------------------------
  // 4. Process each endcap side
  // -----------------------------------------------------------------------
  for (int side : {-1, +1}) {

    struct TsInfo { unsigned origIdx; float eta, phi, z, energy; };
    std::vector<TsInfo> ts;
    ts.reserve(tracksters.size());
    std::vector<int> globalToLocal(tracksters.size(), -1);

    for (unsigned i = 0; i < tracksters.size(); ++i) {
      const auto& t = tracksters[i];
      if (side > 0 && t.barycenter().eta() <= 0) continue;
      if (side < 0 && t.barycenter().eta() >= 0) continue;
      globalToLocal[i] = static_cast<int>(ts.size());
      ts.push_back({i,
                    static_cast<float>(t.barycenter().eta()),
                    static_cast<float>(t.barycenter().phi()),
                    static_cast<float>(t.barycenter().z()),
                    t.raw_energy()});
    }
    const int nTS = static_cast<int>(ts.size());
    if (nTS == 0) continue;

    std::vector<int> sideTrackIdx;
    for (int i = 0; i < nTracks; ++i) {
      if (side > 0 && validTracks[i].eta <= 0) continue;
      if (side < 0 && validTracks[i].eta >= 0) continue;
      sideTrackIdx.push_back(i);
    }
    const int nSideTracks = static_cast<int>(sideTrackIdx.size());

    // -----------------------------------------------------------------------
    // 5. Build edges
    // -----------------------------------------------------------------------
    struct Edge { int u, v; int64_t cost; };
    std::vector<Edge> trackTsEdges, tsTsEdges;

    // Track → TS
    for (int ti = 0; ti < nSideTracks; ++ti) {
      const auto& trk = validTracks[sideTrackIdx[ti]];
      for (unsigned globalJ : findNeighbours(trk.eta, trk.phi)) {
        int localJ = globalToLocal[globalJ];
        if (localJ < 0) continue;
        float score = computeScore(trk.pt, trk.eta, trk.phi, trk.p,
                                   ts[localJ].eta, ts[localJ].phi, ts[localJ].energy);
        int64_t cost = static_cast<int64_t>((score - trackTsScoreShift_) * trackTsScoreWeight_ * 1000);
        trackTsEdges.push_back({ti, localJ, cost});
      }
    }

    // TS → TS (outward direction only)
    for (int i = 0; i < nTS; ++i) {
      for (unsigned globalJ : findNeighbours(ts[i].eta, ts[i].phi)) {
        int localJ = globalToLocal[globalJ];
        if (localJ < 0 || localJ == i) continue;
        if (std::abs(ts[localJ].z) <= std::abs(ts[i].z)) continue;
        float deta = tsAllProp[ts[localJ].origIdx].Eta() - tsAllProp[ts[i].origIdx].Eta();
        float dphi = reco::deltaPhi((float)tsAllProp[ts[localJ].origIdx].Phi(),
                                    (float)tsAllProp[ts[i].origIdx].Phi());
        float dr   = std::sqrt(deta*deta + dphi*dphi);
        float en   = std::max(ts[i].energy, ts[localJ].energy);
        float eta  = (ts[i].energy >= ts[localJ].energy) ? ts[i].eta : ts[localJ].eta;
        int64_t cost = static_cast<int64_t>((dr / normTracksters(en, eta) - tsTsScoreShift_) * tsTsScoreWeight_ * 1000);
        tsTsEdges.push_back({i, localJ, cost});
      }
    }

    // -----------------------------------------------------------------------
    // 6. Build min-cost flow graph
    // Node layout:
    //   SRC=0
    //   TRACK nodes: [1, nSideTracks]
    //   TS_IN  nodes: [nSideTracks+1, nSideTracks+nTS]
    //   TS_OUT nodes: [nSideTracks+nTS+1, nSideTracks+2*nTS]
    //   SNK = nSideTracks+2*nTS+1
    // -----------------------------------------------------------------------
    const int SRC          = 0;
    const int TRACK_OFFSET = 1;
    const int TS_IN_OFFSET = TRACK_OFFSET + nSideTracks;
    const int TS_OUT_OFFSET= TS_IN_OFFSET + nTS;
    const int SNK          = TS_OUT_OFFSET + nTS;
    const int N_NODES      = SNK + 1;

    MinCostFlow mcf(N_NODES);

    // SRC → Track (large capacity, negative cost to incentivise track usage)
    for (int ti = 0; ti < nSideTracks; ++ti)
      mcf.addArc(SRC, TRACK_OFFSET + ti, nTS, -100000);  // -100 * 1000 scaling

    // SRC → TS_IN (neutral path, zero cost)
    for (int j = 0; j < nTS; ++j)
      mcf.addArc(SRC, TS_IN_OFFSET + j, 1, 0);

    // Track → TS_IN
    for (const auto& e : trackTsEdges)
      mcf.addArc(TRACK_OFFSET + e.u, TS_IN_OFFSET + e.v, 1, e.cost);

    // TS_IN → TS_OUT (capacity=1 enforces exclusivity)
    for (int j = 0; j < nTS; ++j)
      mcf.addArc(TS_IN_OFFSET + j, TS_OUT_OFFSET + j, 1, 0);

    // TS_OUT → TS_IN (TS→TS chaining)
    for (const auto& e : tsTsEdges)
      mcf.addArc(TS_OUT_OFFSET + e.u, TS_IN_OFFSET + e.v, 1, e.cost);

    // TS_OUT → SNK
    int64_t neutralCost = static_cast<int64_t>(neutralPenalty_ * 1000);
    for (int j = 0; j < nTS; ++j)
      mcf.addArc(TS_OUT_OFFSET + j, SNK, 1, neutralCost);

    // Track → SNK (track-only candidates)
    for (int ti = 0; ti < nSideTracks; ++ti)
      mcf.addArc(TRACK_OFFSET + ti, SNK, 1, neutralCost);

    // Supplies: push exactly nTS units through the network
    mcf.setNodeSupply(SRC,  nTS);
    mcf.setNodeSupply(SNK, -nTS);

    // -----------------------------------------------------------------------
    // 7. Solve
    // -----------------------------------------------------------------------
    if (mcf.solve() != MinCostFlow::OPTIMAL) {
      edm::LogWarning("MCFInterpretationAlgo") << "Min-cost flow did not find an optimal solution";
      continue;
    }

    // -----------------------------------------------------------------------
    // 8. Decode flow → adjacency map
    // -----------------------------------------------------------------------
//    std::unordered_map<int, std::vector<int>> usedOut;
//    for (int arc = 0; arc < mcf.numArcs(); ++arc) {
//      if (mcf.flow(arc) > 0)
//        usedOut[mcf.tail(arc)].push_back(mcf.head(arc));
//    }
//
//    // Follow a flow chain from startNode, collecting local TS indices
//    auto followChain = [&](int startNode) -> std::vector<int> {
//      std::vector<int> tsChain;
//      int cur = startNode;
//      while (cur != SNK) {
//        if (cur >= TS_IN_OFFSET && cur < TS_OUT_OFFSET)
//          tsChain.push_back(cur - TS_IN_OFFSET);
//        auto it = usedOut.find(cur);
//        if (it == usedOut.end() || it->second.empty()) break;
//        cur = it->second[0];
//      }
//      return tsChain;
//    };
//
//    // Helper: push merged or single trackster and set resultCandidate
//    auto pushCandidate = [&](int origTrackIdx, const std::vector<int>& localTsIndices) {
//      if (localTsIndices.size() == 1) {
//        if (origTrackIdx >= 0)
//          resultCandidate[origTrackIdx] = static_cast<int>(resultTracksters.size());
//        resultTracksters.push_back(tracksters[ts[localTsIndices[0]].origIdx]);
//      } else {
//        Trackster merged;
//        bool isHadron = false;
//        for (int idx : localTsIndices) {
//          merged.mergeTracksters(tracksters[ts[idx].origIdx]);
//          if (tracksters[ts[idx].origIdx].isHadronic()) isHadron = true;
//        }
//        if (origTrackIdx >= 0) {
//          resultCandidate[origTrackIdx] = static_cast<int>(resultTracksters.size());
//          merged.setIdProbability(
//              isHadron ? ticl::Trackster::ParticleType::charged_hadron
//                       : ticl::Trackster::ParticleType::electron, 1.f);
//        }
//        resultTracksters.push_back(merged);
//      }
//    };
//
//    // Charged candidates: SRC → Track → ...
//    for (int trkNode : usedOut[SRC]) {
//      if (trkNode < TRACK_OFFSET || trkNode >= TS_IN_OFFSET) continue;
//      int ti = trkNode - TRACK_OFFSET;
//      int origTrackIdx = validTracks[sideTrackIdx[ti]].origIdx;
//
//      std::vector<int> tsSet;
//      for (int v : usedOut[trkNode]) {
//        auto chain = followChain(v);
//        tsSet.insert(tsSet.end(), chain.begin(), chain.end());
//      }
//
//      if (tsSet.empty()) {
//        // Track-only: no trackster linked, record index but push nothing
//        resultCandidate[origTrackIdx] = static_cast<int>(resultTracksters.size());
//      } else {
//        pushCandidate(origTrackIdx, tsSet);
//      }
//    }
//
//    // Neutral candidates: SRC → TS_IN directly
//    for (int startNode : usedOut[SRC]) {
//      if (startNode < TS_IN_OFFSET || startNode >= TS_OUT_OFFSET) continue;
//      auto chain = followChain(startNode);
//      if (!chain.empty())
//        pushCandidate(-1, chain);  // -1 = no track
//    }


std::unordered_map<int, std::vector<int>> usedOut;
for (int arc = 0; arc < mcf.numArcs(); ++arc) {
    if (mcf.flow(arc) > 0)
        usedOut[mcf.tail(arc)].push_back(mcf.head(arc));
}

auto followChain = [&](int startNode) -> std::vector<int> {
    std::vector<int> tsChain;
    int cur = startNode;
    while (cur != SNK) {
        if (cur >= TS_IN_OFFSET && cur < TS_OUT_OFFSET)
            tsChain.push_back(cur - TS_IN_OFFSET);
        auto it = usedOut.find(cur);
        if (it == usedOut.end() || it->second.empty()) break;
        cur = it->second[0];
    }
    return tsChain;
};

auto pushCandidate = [&](int origTrackIdx, const std::vector<int>& localTsIndices) {
    if (localTsIndices.size() == 1) {
        const auto& inputTs = tracksters[ts[localTsIndices[0]].origIdx];

        if (inputTs.raw_energy() < 0 || inputTs.regressed_energy() < 0)
            std::cout << "[NEGATIVE] Single-TS candidate"
                      << " origTrackIdx=" << origTrackIdx
                      << " origTsIdx="    << ts[localTsIndices[0]].origIdx
                      << " raw_energy="   << inputTs.raw_energy()
                      << " regr_energy="  << inputTs.regressed_energy()
                      << " eta="          << inputTs.barycenter().eta()
                      << " phi="          << inputTs.barycenter().phi()
                      << std::endl;

        if (origTrackIdx >= 0)
            resultCandidate[origTrackIdx] = static_cast<int>(resultTracksters.size());
        resultTracksters.push_back(inputTs);

    } else {
        Trackster merged;
        bool isHadron = false;

        std::cout << "[MERGE] origTrackIdx=" << origTrackIdx
                  << " merging " << localTsIndices.size() << " tracksters:" << std::endl;

        for (int idx : localTsIndices) {
            const auto& inputTs = tracksters[ts[idx].origIdx];

            std::cout << "  + ts[" << idx << "] origIdx=" << ts[idx].origIdx
                      << " raw_energy="  << inputTs.raw_energy()
                      << " regr_energy=" << inputTs.regressed_energy()
                      << " eta="         << inputTs.barycenter().eta()
                      << " phi="         << inputTs.barycenter().phi()
                      << " time="        << inputTs.time()
                      << " isHadronic="  << inputTs.isHadronic()
                      << std::endl;

            float eBefore = merged.raw_energy();
            merged.mergeTracksters(inputTs);
            float eAfter = merged.raw_energy();

            if (eAfter < 0 || eAfter < eBefore)
                std::cout << "  [!] Energy changed after merge: "
                          << eBefore << " -> " << eAfter
                          << " (added ts raw_energy=" << inputTs.raw_energy() << ")"
                          << std::endl;
        }

        std::cout << "  => merged raw="  << merged.raw_energy()
                  << " regr=" << merged.regressed_energy()
                  << " isHadron=" << isHadron << std::endl;

        if (merged.raw_energy() < 0 || merged.regressed_energy() < 0)
            std::cout << "[NEGATIVE] Merged candidate"
                      << " origTrackIdx=" << origTrackIdx
                      << " nTs="          << localTsIndices.size()
                      << " raw_energy="   << merged.raw_energy()
                      << " regr_energy="  << merged.regressed_energy()
                      << std::endl;

        if (origTrackIdx >= 0) {
            resultCandidate[origTrackIdx] = static_cast<int>(resultTracksters.size());
            merged.setIdProbability(
                isHadron ? ticl::Trackster::ParticleType::charged_hadron
                         : ticl::Trackster::ParticleType::electron, 1.f);
        }
        resultTracksters.push_back(merged);
    }
};

// Charged candidates
for (int trkNode : usedOut[SRC]) {
    if (trkNode < TRACK_OFFSET || trkNode >= TS_IN_OFFSET) continue;
    int ti = trkNode - TRACK_OFFSET;
    int origTrackIdx = validTracks[sideTrackIdx[ti]].origIdx;

    std::vector<int> tsSet;
    for (int v : usedOut[trkNode]) {
        auto chain = followChain(v);
        tsSet.insert(tsSet.end(), chain.begin(), chain.end());
    }

    std::cout << "[CHARGED] origTrackIdx=" << origTrackIdx
              << " nLinkedTs=" << tsSet.size() << std::endl;

    if (tsSet.empty()) {
        std::cout << "  track-only candidate, no trackster" << std::endl;
        resultCandidate[origTrackIdx] = -1; // static_cast<int>(resultTracksters.size());
    } else {
        pushCandidate(origTrackIdx, tsSet);
    }
}

// Neutral candidates
for (int startNode : usedOut[SRC]) {
    if (startNode < TS_IN_OFFSET || startNode >= TS_OUT_OFFSET) continue;
    auto chain = followChain(startNode);

    std::cout << "[NEUTRAL] nLinkedTs=" << chain.size() << std::endl;

    if (!chain.empty())
        pushCandidate(-1, chain);
}


  }  // end side loop
}

// ---------------------------------------------------------------------------
// fillPSetDescription
// ---------------------------------------------------------------------------
void MCFInterpretationAlgo::fillPSetDescription(edm::ParameterSetDescription& desc) {
  desc.add<double>("drCut", 0.02);                // max dR for graph edges
  desc.add<double>("tsTsScoreShift", 1.0);        // shift applied to TS-TS edge cost
  desc.add<double>("trackTsScoreShift", 1.0);     // shift applied to track-TS edge cost
  desc.add<double>("tsTsScoreWeight", 1.0);       // scale for TS-TS edge cost
  desc.add<double>("trackTsScoreWeight", 1.0);    // scale for track-TS edge cost
  desc.add<double>("neutralPenalty", 1.0);        // cost for unlinked (neutral) flow
  TICLInterpretationAlgoBase<reco::Track>::fillPSetDescription(desc);
}


DEFINE_EDM_PLUGIN(TICLGeneralInterpretationPluginFactory,
                  ticl::MCFInterpretationAlgo,
                  "MCFInterpretationAlgo");
