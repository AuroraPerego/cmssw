#include "MCFwithNNInterpretationAlgo.h"
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
MCFwithNNInterpretationAlgo::MCFwithNNInterpretationAlgo(const edm::ParameterSet& conf,
                                                         TICLONNXGlobalCache const* cache)
    : TICLInterpretationAlgoBase<reco::Track>(conf, cache),
      drCut_(conf.getParameter<double>("drCut")),
      tsTsScoreShift_(conf.getParameter<double>("tsTsScoreShift")),
      trackTsScoreShift_(conf.getParameter<double>("trackTsScoreShift")),
      tsTsScoreWeight_(conf.getParameter<double>("tsTsScoreWeight")),
      trackTsScoreWeight_(conf.getParameter<double>("trackTsScoreWeight")),
      neutralPenalty_(conf.getParameter<int>("neutralPenalty")),
      tracksterInit_(conf.getParameter<int>("tracksterInit")),
      trackInit_(conf.getParameter<int>("trackInit")) {
  const std::string trackModel = conf.getParameter<std::string>("onnxTrackModel");
  const std::string tracksterModel = conf.getParameter<std::string>("onnxTracksterModel");

  if (cache_ != nullptr) {
    onnxSessionTracks_ = cache_->getByModelPathString(trackModel);
    onnxSessionTracksters_ = cache_->getByModelPathString(tracksterModel);
  }
}

// ---------------------------------------------------------------------------
// initialize + buildLayers
// ---------------------------------------------------------------------------
void MCFwithNNInterpretationAlgo::initialize(const HGCalDDDConstants* hgcons,
                                             const hgcal::RecHitTools rhtools,
                                             const edm::ESHandle<MagneticField> bfieldH,
                                             const edm::ESHandle<Propagator> propH) {
  hgcons_ = hgcons;
  rhtools_ = rhtools;
  bfield_ = bfieldH;
  propagator_ = propH;
  buildLayers();
}

void MCFwithNNInterpretationAlgo::buildLayers() {
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
float MCFwithNNInterpretationAlgo::normTracks(float x, float y) const {
  y = std::abs(y);
  constexpr float a = 2.29516386e-01f, b = -8.05371532e+02f, c = -6.45573586e-01f;
  constexpr float d = 8.05370082e+02f, e = 1.00000033e+00f, f = -1.07458042e-04f;
  constexpr float g = 4.79631755e-01f, h = 1.21330763e+00f;
  if (x > 200.f)
    x = 200.f;
  if (y < 1.7f)
    y = 1.7f;
  return std::max(a + b * x + c * y + d * std::pow(x, e) + f * x * y + g * std::pow(y, h), 0.008f);
}

float MCFwithNNInterpretationAlgo::normTracksters(float x, float y) const {
  y = std::abs(y);
  constexpr float a = 1.49662496e+00f, b = 4.04830636e+01f, c = -1.19840947e+01f;
  constexpr float d = -4.04834456e+01f, e = 9.99984896e-01f, f = -1.67329993e-03f;
  constexpr float g = 1.06828890e+01f, h = 1.09862164e+00f;
  if (x > 200.f)
    x = 200.f;
  return a + b * x + c * y + d * std::pow(x, e) + g * std::pow(y, h) + f * x * y;
}

float MCFwithNNInterpretationAlgo::computeScore(
    float refPt, float refEta, float refPhi, float refP, float tsEta, float tsPhi, float tsEnergy) const {
  constexpr float wEp = 0.5f;
  float rEpNorm = (tsEnergy - refP) / refP;
  float dphi = reco::deltaPhi(refPhi, tsPhi);
  float deta = refEta - tsEta;
  float pullR2 = (deta * deta + dphi * dphi) / std::pow(normTracks(refPt, refEta), 2);
  return std::sqrt(wEp * rEpNorm * rEpNorm + (1.f - wEp) * pullR2);
}

// ---------------------------------------------------------------------------
// makeCandidates
// ---------------------------------------------------------------------------
void MCFwithNNInterpretationAlgo::makeCandidates(const Inputs& input,
                                                 edm::Handle<MtdHostCollection> inputTimingh,
                                                 std::vector<Trackster>& resultTracksters,
                                                 std::vector<int>& resultCandidate) {
  const auto& tracksters = input.tracksters;
  const auto tkH = input.tracksHandle;
  const auto& tracks = *tkH;
  const auto& maskTracks = input.maskedTracks;
  const float drCut2 = drCut_ * drCut_;

  auto bFieldProd = bfield_.product();
  const Propagator& prop = *propagator_;

  // -----------------------------------------------------------------------
  // 1. Propagate all tracksters to HGCal front face, fill tiles
  // -----------------------------------------------------------------------
  struct TsInfo {
    unsigned origIdx;
    float eta, phi, energy;
    float x, y, z;
    float time, timeErr;
    int pid;
    // Vector tPoint;
  };
  std::array<TICLLayerTile, 2> tracksterPropTiles = {};
  std::vector<std::vector<TsInfo>> tsAllProp(2);
  tsAllProp[0].reserve(tracksters.size());
  tsAllProp[1].reserve(tracksters.size());

  const float zVal_layer1 = hgcons_->waferZ(1, true);

  for (unsigned i = 0; i < tracksters.size(); ++i) {
    const auto& t = tracksters[i];
    const Vector& baryc = t.barycenter();
    /*Vector directnv = baryc.unit();
    float zVal = zVal_layer1 * (baryc.Z() > 0 ? 1.f : -1.f);
    float par = (zVal - baryc.Z()) / directnv.Z();
    Vector tPoint(par * directnv.X() + baryc.X(), par * directnv.Y() + baryc.Y(), zVal);
    if (tPoint.Eta() > 0)
      tracksterPropTiles[1].fill(tPoint.Eta(), tPoint.Phi(), i);
    else if (tPoint.Eta() < 0)
      tracksterPropTiles[0].fill(tPoint.Eta(), tPoint.Phi(), i);
    tsAllProp.emplace_back(tPoint);*/
    const auto& probs = t.id_probabilities();
    const int pid = static_cast<int>(std::max_element(probs.begin(), probs.end()) - probs.begin());  // FIXME
    if (baryc.eta() >= 0) {
      tracksterPropTiles[1].fill(baryc.eta(), baryc.phi(), tsAllProp[1].size());
      tsAllProp[1].push_back({i,
                              t.barycenter().eta(),
                              t.barycenter().phi(),
                              t.raw_energy(),
                              t.barycenter().x(),
                              t.barycenter().y(),
                              t.barycenter().z(),
                              t.time(),
                              t.timeError(),
                              pid});
    } else {
      tracksterPropTiles[0].fill(baryc.eta(), baryc.phi(), tsAllProp[0].size());
      tsAllProp[0].push_back({i,
                              t.barycenter().eta(),
                              t.barycenter().phi(),
                              t.raw_energy(),
                              t.barycenter().x(),
                              t.barycenter().y(),
                              t.barycenter().z(),
                              t.time(),
                              t.timeError(),
                              pid});
    }
  }

  // -----------------------------------------------------------------------
  // 2. Propagate valid tracks to HGCal front face
  // -----------------------------------------------------------------------
  struct TrackInfo {
    int origIdx;
    float eta, phi;
    double pt, p;
  };
  std::vector<std::vector<TrackInfo>> validTracks(2);
  validTracks[0].reserve(tracks.size());
  validTracks[1].reserve(tracks.size());

  for (unsigned i = 0; i < tracks.size(); ++i) {
    if (!maskTracks[i])
      continue;
    const auto& tk = tracks[i];
    if (tk.pt() < 1.f || tk.p() < 2.f)
      continue;
    if (std::abs(tk.eta()) < 1.5f || std::abs(tk.eta()) > 3.0f)
      continue;

    int iSide = int(tk.eta() > 0);
    FreeTrajectoryState fts = tk.outerOk() ? trajectoryStateTransform::outerFreeState(tk, bFieldProd)
                                           : trajectoryStateTransform::initialFreeState(tk, bFieldProd);

    const auto& tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
    if (!tsos.isValid())
      continue;

    GlobalPoint pp = tsos.globalPosition();
    if (pp.eta() < 0)
      validTracks[0].push_back({static_cast<int>(i), pp.eta(), pp.phi(), tk.pt(), tk.p()});
    else
      validTracks[1].push_back({static_cast<int>(i), pp.eta(), pp.phi(), tk.pt(), tk.p()});
  }

  // -----------------------------------------------------------------------
  // 3. findNeighbours lambda (tile-based, mirrors findTrackstersInWindow)
  // -----------------------------------------------------------------------
  auto findNeighbours = [&](float seed_eta, float seed_phi, int side) -> std::vector<unsigned> {
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
          float deta = tsAllProp[side][t_i].eta - seed_eta;
          float dphi = reco::deltaPhi(tsAllProp[side][t_i].phi, seed_phi);
          if (deta * deta + dphi * dphi < drCut2)
            result.push_back(t_i);
        }
      }
    }
    return result;
  };
  // -----------------------------------------------------------------------
  // 4. Process each endcap side
  // -----------------------------------------------------------------------
  for (int side : {0, 1}) {
    struct TsInfo {
      unsigned origIdx;
      float eta, phi, energy;
      float x, y, z;
      float time, timeErr;
      int pid;
    };

    const auto& ts = tsAllProp[side];
    const int nTS = static_cast<int>(ts.size());
    if (nTS == 0)
      continue;

    const auto& sideTracks = validTracks[side];
    const int nSideTracks = static_cast<int>(sideTracks.size());

    constexpr float C_CM_PER_NS = 29.9792458f;
    bool useMTDTiming = inputTimingh.isValid();

    // -----------------------------------------------------------------------
    // 6. Build edges in two-pass: collect features, batch infer, fill costs
    // -----------------------------------------------------------------------
    struct Edge {
      unsigned int u, v;
      int64_t cost;
    };
    std::vector<Edge> trackTsEdges, tsTsEdges;

    constexpr int TRACK_TS_NFEAT = 18;
    constexpr int TS_TS_NFEAT = 21;

    std::vector<std::pair<unsigned int, unsigned int>> trkTsRaw;

    cms::Ort::FloatArrays trkTsFeatsInput(1);
    auto& trkTsFeats = trkTsFeatsInput[0];

    for (int ti = 0; ti < nSideTracks; ++ti) {
      const auto& trk = sideTracks[ti];
      //const auto& tkOrig  = tracks[trk.origIdx];

      float trkTime = 0.f, trkTimeErr = -1.f;
      float trkMtdX = 0.f, trkMtdY = 0.f, trkMtdZ = 0.f;
      if (useMTDTiming) {
        auto const& tv = (*inputTimingh).const_view();
        trkTime = tv.time()[trk.origIdx];
        trkTimeErr = tv.timeErr()[trk.origIdx];
        trkMtdX = tv.posInMTD_x()[trk.origIdx];
        trkMtdY = tv.posInMTD_y()[trk.origIdx];
        trkMtdZ = tv.posInMTD_z()[trk.origIdx];
      }

      for (unsigned localJ : findNeighbours(trk.eta, trk.phi, side)) {
        const auto& tsInf = ts[localJ];

        float deltaPhi = reco::deltaPhi(trk.phi, tsInf.phi);
        float deltaEta = trk.eta - tsInf.eta;
        float deltaE = (trk.p - tsInf.energy) / trk.p;
        float deltaR = std::sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);

        float deltaTime = 0.f;
        if (trkTimeErr > 0 && tsInf.timeErr > 0) {
          float dx = trkMtdX - tsInf.x, dy = trkMtdY - tsInf.y, dz = trkMtdZ - tsInf.z;
          float tof = std::sqrt(dx * dx + dy * dy + dz * dz) / C_CM_PER_NS;
          deltaTime = tsInf.time - trkTime - tof;
        }

        // Feature order must match training columns exactly
        trkTsFeats.insert(trkTsFeats.end(),
                          {
                              static_cast<float>(trk.pt),  // refPt
                              static_cast<float>(trk.p),   // refP
                              trk.eta,                     // refEta
                              std::sin(trk.phi),           // sin_refPhi
                              std::cos(trk.phi),           // cos_refPhi
                              trkTime,                     // trk_time
                              trkTimeErr,                  // trk_timeErr
                              tsInf.energy,                // tsEnergy
                              tsInf.eta,                   // tsEta
                              std::sin(tsInf.phi),         // sin_tsPhi
                              std::cos(tsInf.phi),         // cos_tsPhi
                              tsInf.time,                  // tsTime
                              tsInf.timeErr,               // tsTimeErr
                              deltaTime,                   // deltaTime
                              deltaE,                      // deltaE
                              deltaEta,                    // deltaEta
                              deltaPhi,                    // deltaPhi
                              deltaR                       // deltaR
                          });
        trkTsRaw.push_back({ti, localJ});
      }
    }

    std::vector<std::pair<unsigned int, unsigned int>> tsTsRaw;
    cms::Ort::FloatArrays tsTsFeatsInput(1);
    auto& tsTsFeats = tsTsFeatsInput[0];

    for (unsigned i = 0; (int)i < nTS; ++i) {
      for (unsigned localJ : findNeighbours(ts[i].eta, ts[i].phi, side)) {
        if (localJ == i)
          continue;
        if (std::abs(ts[localJ].z) <= std::abs(ts[i].z))
          continue;

        const auto& ts1 = ts[i];
        const auto& ts2 = ts[localJ];

        float deltaPhi = reco::deltaPhi(ts1.phi, ts2.phi);
        float deltaEta = ts1.eta - ts2.eta;
        float deltaR = std::sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
        float deltaE = ts1.energy - ts2.energy;
        float deltaZ = ts1.z - ts2.z;
        float samePid = (ts1.pid == ts2.pid) ? 1.f : 0.f;

        float deltaTime = 0.f;
        if (ts1.timeErr > 0 && ts2.timeErr > 0) {
          float dx = ts1.x - ts2.x, dy = ts1.y - ts2.y, dz = ts1.z - ts2.z;
          float tof = std::sqrt(dx * dx + dy * dy + dz * dz) / C_CM_PER_NS;
          deltaTime = std::abs(ts1.time - ts2.time) - tof;
        }

        // Feature order must match training columns exactly
        tsTsFeats.insert(tsTsFeats.end(),
                         {
                             ts1.energy,         // E1
                             ts1.eta,            // eta1
                             std::sin(ts1.phi),  // sin_phi1
                             std::cos(ts1.phi),  // cos_phi1
                             ts1.z,              // Z1
                             ts1.time,           // time1
                             ts1.timeErr,        // timeErr1
                             ts2.energy,         // E2
                             ts2.eta,            // eta2
                             std::sin(ts2.phi),  // sin_phi2
                             std::cos(ts2.phi),  // cos_phi2
                             ts2.z,              // Z2
                             ts2.time,           // time2
                             ts2.timeErr,        // timeErr2
                             deltaTime,          // deltaTime
                             samePid,            // samePid
                             deltaE,             // deltaE
                             deltaEta,           // deltaEta
                             deltaPhi,           // deltaPhi
                             deltaR,             // deltaR
                             deltaZ              // deltaZ
                         });
        tsTsRaw.push_back({i, localJ});
      }
    }

    const int nTrkTsEdges = static_cast<int>(trkTsRaw.size());
    const int nTsTsEdges = static_cast<int>(tsTsRaw.size());

    cms::Ort::FloatArrays trkTsScores;
    cms::Ort::FloatArrays tsTsScores;

    if (nTrkTsEdges > 0) {
      // Shape: [nEdges, 18]
      onnxSessionTracks_->runInto(
          {"input"}, trkTsFeatsInput, {{static_cast<int64_t>(nTrkTsEdges), TRACK_TS_NFEAT}}, {"score"}, trkTsScores);
    }

    if (nTsTsEdges > 0) {
      // Shape: [nEdges, 21]
      onnxSessionTracksters_->runInto(
          {"input"}, tsTsFeatsInput, {{static_cast<int64_t>(nTsTsEdges), TS_TS_NFEAT}}, {"score"}, tsTsScores);
    }

    trackTsEdges.reserve(nTrkTsEdges);
    for (int k = 0; k < nTrkTsEdges; ++k) {
      int64_t cost = static_cast<int64_t>(-trkTsScores[k][0] * trackTsScoreWeight_ + trackTsScoreShift_);
      trackTsEdges.push_back({trkTsRaw[k].first, trkTsRaw[k].second, cost});
    }

    tsTsEdges.reserve(nTsTsEdges);
    for (int k = 0; k < nTsTsEdges; ++k) {
      int64_t cost = static_cast<int64_t>(-tsTsScores[k][0] * tsTsScoreWeight_ + tsTsScoreShift_);
      tsTsEdges.push_back({tsTsRaw[k].second, tsTsRaw[k].second, cost});
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
    const int SRC = 0;
    const int TRACK_OFFSET = 1;
    const int TS_IN_OFFSET = TRACK_OFFSET + nSideTracks;
    const int TS_OUT_OFFSET = TS_IN_OFFSET + nTS;
    const int SNK = TS_OUT_OFFSET + nTS;
    const int N_NODES = SNK + 1;

    MinCostFlow mcf(N_NODES);

    // SRC -> Track (large capacity, negative cost to incentivise track usage)
    for (int ti = 0; ti < nSideTracks; ++ti)
      mcf.addArc(SRC, TRACK_OFFSET + ti, nTS, tracksterInit_);  // -100 * 1000 scaling

    // SRC -> TS_IN (neutral path, zero cost)
    for (int j = 0; j < nTS; ++j)
      mcf.addArc(SRC, TS_IN_OFFSET + j, 1, trackInit_);

    // Track -> TS_IN
    for (const auto& e : trackTsEdges)
      mcf.addArc(TRACK_OFFSET + e.u, TS_IN_OFFSET + e.v, 1, e.cost);

    // TS_IN -> TS_OUT (capacity=1 enforces exclusivity)
    for (int j = 0; j < nTS; ++j)
      mcf.addArc(TS_IN_OFFSET + j, TS_OUT_OFFSET + j, 1, 0);

    // TS_OUT -> TS_IN (TS->TS chaining)
    for (const auto& e : tsTsEdges)
      mcf.addArc(TS_OUT_OFFSET + e.u, TS_IN_OFFSET + e.v, 1, e.cost);

    // TS_OUT -> SNK
    for (int j = 0; j < nTS; ++j)
      mcf.addArc(TS_OUT_OFFSET + j, SNK, 1, neutralPenalty_);

    // Track -> SNK (track-only candidates)
    for (int ti = 0; ti < nSideTracks; ++ti)
      mcf.addArc(TRACK_OFFSET + ti, SNK, 1, neutralPenalty_);

    // Supplies: push exactly nTS units through the network
    mcf.setNodeSupply(SRC, nTS);
    mcf.setNodeSupply(SNK, -nTS);

    // -----------------------------------------------------------------------
    // 7. Solve
    // -----------------------------------------------------------------------
    if (mcf.solve() != MinCostFlow::OPTIMAL) {
      edm::LogWarning("MCFwithNNInterpretationAlgo") << "Min-cost flow did not find an optimal solution";
      continue;
    }

    // -----------------------------------------------------------------------
    // 8. Decode flow → adjacency map
    // -----------------------------------------------------------------------
    std::unordered_map<int, std::vector<int>> usedOut;
    for (int arc = 0; arc < mcf.numArcs(); ++arc) {
      if (mcf.flow(arc) > 0)
        usedOut[mcf.tail(arc)].push_back(mcf.head(arc));
    }

    // Follow a flow chain from startNode, collecting local TS indices
    auto followChain = [&](int startNode) -> std::vector<int> {
      std::vector<int> tsChain;
      int cur = startNode;
      while (cur != SNK) {
        if (cur >= TS_IN_OFFSET && cur < TS_OUT_OFFSET)
          tsChain.push_back(cur - TS_IN_OFFSET);
        auto it = usedOut.find(cur);
        if (it == usedOut.end() || it->second.empty())
          break;
        cur = it->second[0];
      }
      return tsChain;
    };

    // Helper: push merged or single trackster and set resultCandidate
    auto pushCandidate = [&](int origTrackIdx, const std::vector<int>& localTsIndices) {
      if (localTsIndices.size() == 1) {
        if (origTrackIdx >= 0)
          resultCandidate[origTrackIdx] = static_cast<int>(resultTracksters.size());
        resultTracksters.push_back(tracksters[ts[localTsIndices[0]].origIdx]);
      } else {
        Trackster merged;
        bool isHadron = false;
        for (int idx : localTsIndices) {
          merged.mergeTracksters(tracksters[ts[idx].origIdx]);
          if (tracksters[ts[idx].origIdx].isHadronic())
            isHadron = true;
        }
        if (origTrackIdx >= 0) {
          resultCandidate[origTrackIdx] = static_cast<int>(resultTracksters.size());
          merged.setIdProbability(
              isHadron ? ticl::Trackster::ParticleType::charged_hadron : ticl::Trackster::ParticleType::electron, 1.f);
        }
        resultTracksters.push_back(merged);
      }
    };

    // Charged candidates: SRC → Track → ...
    for (int trkNode : usedOut[SRC]) {
      if (trkNode < TRACK_OFFSET || trkNode >= TS_IN_OFFSET)
        continue;
      int ti = trkNode - TRACK_OFFSET;
      int origTrackIdx = sideTracks[ti].origIdx;

      std::vector<int> tsSet;
      for (int v : usedOut[trkNode]) {
        auto chain = followChain(v);
        tsSet.insert(tsSet.end(), chain.begin(), chain.end());
      }

      if (tsSet.empty()) {
        // Track-only: no trackster linked, record index but push nothing
        resultCandidate[origTrackIdx] = static_cast<int>(resultTracksters.size());
      } else {
        pushCandidate(origTrackIdx, tsSet);
      }
    }

    // Neutral candidates: SRC → TS_IN directly
    for (int startNode : usedOut[SRC]) {
      if (startNode < TS_IN_OFFSET || startNode >= TS_OUT_OFFSET)
        continue;
      auto chain = followChain(startNode);
      if (!chain.empty())
        pushCandidate(-1, chain);  // -1 = no track
    }
  }  // end side loop
}

// ---------------------------------------------------------------------------
// fillPSetDescription
// ---------------------------------------------------------------------------
void MCFwithNNInterpretationAlgo::fillPSetDescription(edm::ParameterSetDescription& desc) {
  desc.add<double>("drCut", 0.02);              // max dR for graph edges
  desc.add<double>("tsTsScoreShift", 1.0);      // shift applied to TS-TS edge cost
  desc.add<double>("trackTsScoreShift", 1.0);   // shift applied to track-TS edge cost
  desc.add<double>("tsTsScoreWeight", 1.0);     // scale for TS-TS edge cost
  desc.add<double>("trackTsScoreWeight", 1.0);  // scale for track-TS edge cost
  desc.add<int>("neutralPenalty", 1);           // cost for unlinked (neutral) flow
  desc.add<int>("tracksterInit", 0);            // cost to start a neutral
  desc.add<int>("trackInit", 0);                // cost to start a charged
  desc.add<std::string>("onnxTrackModel", "");
  desc.add<std::string>("onnxTracksterModel", "");
  TICLInterpretationAlgoBase<reco::Track>::fillPSetDescription(desc);
}

DEFINE_EDM_PLUGIN(TICLGeneralInterpretationPluginFactory,
                  ticl::MCFwithNNInterpretationAlgo,
                  "MCFwithNNInterpretationAlgo");
