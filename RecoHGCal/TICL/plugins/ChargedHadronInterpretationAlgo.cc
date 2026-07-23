#include "RecoHGCal/TICL/interface/TICLInterpretationAlgoBase.h"
#include "RecoHGCal/TICL/plugins/ChargedHadronInterpretationAlgo.h"
#include "RecoHGCal/TICL/interface/MinCostFlow.h"
#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <unordered_map>

using namespace ticl;

using Vector = ticl::Trackster::Vector;

ChargedHadronInterpretationAlgo::~ChargedHadronInterpretationAlgo() {}

ChargedHadronInterpretationAlgo::ChargedHadronInterpretationAlgo(const edm::ParameterSet &conf,
                                                                 TICLONNXGlobalCache const* cache)
    : TICLInterpretationAlgoBase<reco::Track>(conf, cache),
      del_tk_ts_layer1_(conf.getParameter<double>("delta_tk_ts_layer1")),
      del_tk_ts_int_(conf.getParameter<double>("delta_tk_ts_interface")),
      timing_quality_threshold_(conf.getParameter<double>("timing_quality_threshold")),
      energy_overshoot_fraction_(conf.getParameter<double>("energy_overshoot_fraction")),
      energy_overshoot_max_(conf.getParameter<double>("energy_overshoot_max")),
      drCut_(conf.getParameter<double>("drCut")),
      tsTsScoreShift_(conf.getParameter<double>("tsTsScoreShift")),
      trackTsScoreShift_(conf.getParameter<double>("trackTsScoreShift")),
      tsTsScoreWeight_(conf.getParameter<double>("tsTsScoreWeight")),
      trackTsScoreWeight_(conf.getParameter<double>("trackTsScoreWeight")),
      inputNames_({"input"}),
      outputNames_({"score"}) {
  const std::string trackModel = conf.getParameter<std::string>("onnxTrackModel");
  const std::string tracksterModel = conf.getParameter<std::string>("onnxTracksterModel");

  if (cache_ != nullptr) {
    onnxSessionTracks_ = cache_->getByModelPathString(trackModel);
    onnxSessionTracksters_ = cache_->getByModelPathString(tracksterModel);
  }
  if (onnxSessionTracks_==nullptr)
    std::cout << "ERROR onnxSessionTracks_ is nullptr!!!\n" ;
  if (onnxSessionTracksters_==nullptr)
    std::cout << "ERROR onnxSessionTracksters_ is nullptr!!!\n" ;
}

void ChargedHadronInterpretationAlgo::initialize(const HGCalDDDConstants *hgcons,
                                                 const hgcal::RecHitTools rhtools,
                                                 const edm::ESHandle<MagneticField> bfieldH,
                                                 const edm::ESHandle<Propagator> propH) {
  hgcons_ = hgcons;
  rhtools_ = rhtools;

  bfield_ = bfieldH;
  propagator_ = propH;

  buildLayers();
}

void ChargedHadronInterpretationAlgo::buildLayers() {
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
// makeCandidates - MCF+NN based linking
// ---------------------------------------------------------------------------
void ChargedHadronInterpretationAlgo::makeCandidates(const Inputs &input,
                                                     edm::Handle<MtdHostCollection> inputTiming_h,
                                                     std::vector<Trackster> &resultTracksters,
                                                     std::vector<int> &resultCandidate,
                                                     std::vector<bool> &maskedTracksters) {
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
  };
  std::array<TICLLayerTile, 2> tracksterPropTiles = {};
  std::vector<std::vector<TsInfo>> tsAllProp(2);
  tsAllProp[0].reserve(tracksters.size());
  tsAllProp[1].reserve(tracksters.size());

  const float zVal_layer1 = hgcons_->waferZ(1, true);

  for (unsigned i = 0; i < tracksters.size(); ++i) {
    const auto& t = tracksters[i];
    const Vector& baryc = t.barycenter();
    const auto& probs = t.id_probabilities();
    const int pid = static_cast<int>(std::max_element(probs.begin(), probs.end()) - probs.begin());
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
  // 3. findNeighbours lambda (tile-based)
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
    const auto& ts = tsAllProp[side];
    const int nTS = static_cast<int>(ts.size());
    if (nTS == 0)
      continue;

    const auto& sideTracks = validTracks[side];
    const int nSideTracks = static_cast<int>(sideTracks.size());

    constexpr float C_CM_PER_NS = 29.9792458f;
    bool useMTDTiming = inputTiming_h.isValid();

    // -----------------------------------------------------------------------
    // 5. Build edges: collect features, batch infer, fill costs
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

      float trkTime = 0.f, trkTimeErr = -1.f;
      float trkMtdX = 0.f, trkMtdY = 0.f, trkMtdZ = 0.f;
      if (useMTDTiming) {
        auto const& tv = (*inputTiming_h).const_view();
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
      onnxSessionTracks_->runInto(
          inputNames_, trkTsFeatsInput, {{static_cast<int64_t>(nTrkTsEdges), TRACK_TS_NFEAT}}, outputNames_, trkTsScores);
    }

    if (nTsTsEdges > 0) {
      onnxSessionTracksters_->runInto(
          inputNames_, tsTsFeatsInput, {{static_cast<int64_t>(nTsTsEdges), TS_TS_NFEAT}}, outputNames_, tsTsScores);
    }

    trackTsEdges.reserve(nTrkTsEdges);
    for (int k = 0; k < nTrkTsEdges; ++k) {
      int64_t cost = static_cast<int64_t>(-trkTsScores[k][0] * trackTsScoreWeight_ + trackTsScoreShift_);
      trackTsEdges.push_back({trkTsRaw[k].first, trkTsRaw[k].second, cost});
    }

    tsTsEdges.reserve(nTsTsEdges);
    for (int k = 0; k < nTsTsEdges; ++k) {
      int64_t cost = static_cast<int64_t>(-tsTsScores[k][0] * tsTsScoreWeight_ + tsTsScoreShift_);
      tsTsEdges.push_back({tsTsRaw[k].first, tsTsRaw[k].second, cost});
    }

    // -----------------------------------------------------------------------
    // 6. Apply masking: mark pre-consumed tracksters as unavailable
    // -----------------------------------------------------------------------
    std::vector<bool> chargedMask(tracksters.size(), true);
    for (size_t i = 0; i < tracksters.size() && i < maskedTracksters.size(); ++i)
      if (maskedTracksters[i])
        chargedMask[i] = false;

    // -----------------------------------------------------------------------
    // 7. Simple linking: for each track, collect nearby tracksters using NN scores
    // -----------------------------------------------------------------------
    for (int ti = 0; ti < nSideTracks; ++ti) {
      int origTrackIdx = sideTracks[ti].origIdx;
      std::vector<unsigned int> linkedTracksters;
      float totalEnergy = 0.f;

      // Find all track-trackster edges involving this track
      for (int k = 0; k < nTrkTsEdges; ++k) {
        if (trackTsEdges[k].u == (unsigned int)ti) {
          unsigned int localTsIdx = trackTsEdges[k].v;
          unsigned int origTsIdx = ts[localTsIdx].origIdx;

          // Check if trackster is available and passes basic criteria
          if (chargedMask[origTsIdx]) {
            linkedTracksters.push_back(origTsIdx);
            chargedMask[origTsIdx] = false;
            totalEnergy += tracksters[origTsIdx].raw_energy();
          }
        }
      }

      // Create candidate and update result
      if (!linkedTracksters.empty()) {
        if (linkedTracksters.size() == 1) {
          resultCandidate[origTrackIdx] = static_cast<int>(resultTracksters.size());
          resultTracksters.push_back(tracksters[linkedTracksters[0]]);
        } else {
          Trackster merged;
          bool isHadron = false;
          for (unsigned tsIdx : linkedTracksters) {
            merged.mergeTracksters(tracksters[tsIdx]);
            if (tracksters[tsIdx].isHadronic())
              isHadron = true;
          }
          resultCandidate[origTrackIdx] = static_cast<int>(resultTracksters.size());
          merged.setIdProbability(
              isHadron ? ticl::Trackster::ParticleType::charged_hadron : ticl::Trackster::ParticleType::electron, 1.f);
          resultTracksters.push_back(merged);
        }
      }
    }

    // -----------------------------------------------------------------------
    // 8. Add neutral tracksters (not linked to any track)
    // -----------------------------------------------------------------------
    for (size_t i = 0; i < tracksters.size(); ++i) {
      if (chargedMask[i]) {
        resultTracksters.push_back(tracksters[i]);
      }
    }

  }  // end side loop

  // Update maskedTracksters to reflect which tracksters were consumed
  for (size_t i = 0; i < tracksters.size(); ++i) {
    if (!chargedMask[i] && i < maskedTracksters.size()) {
      maskedTracksters[i] = true;
    }
  }
}

void ChargedHadronInterpretationAlgo::makeOpinions(const Inputs &input,
                                                   edm::Handle<MtdHostCollection> inputTiming_h,
                                                   std::vector<Trackster> &hypothesisTracksters,
                                                   std::vector<Hypothesis> &hypotheses) {
  // Reuse the linking from makeCandidates: run makeCandidates on local buffers and emit one
  // charged-hadron hypothesis per track-linked merged trackster, dropping the neutral
  // leftovers makeCandidates appends (the arbiter derives neutrals itself).
  std::vector<Trackster> localTracksters;
  std::vector<int> localCandidate(input.tracksHandle->size(), -1);
  std::vector<bool> localMask;  // nothing pre-consumed in opinion mode
  makeCandidates(input, inputTiming_h, localTracksters, localCandidate, localMask);

  for (size_t iTrack = 0; iTrack < localCandidate.size(); ++iTrack) {
    if (localCandidate[iTrack] < 0)
      continue;
    Hypothesis h;
    h.type = Hypothesis::Type::ChargedHadron;
    h.score = 0.5f;  // rule-based placeholder; type priority drives the arbitration
    h.trackIdx = static_cast<int>(iTrack);
    h.tracksterIdx = static_cast<int>(hypothesisTracksters.size());
    hypothesisTracksters.push_back(localTracksters[localCandidate[iTrack]]);
    hypotheses.push_back(h);
  }

  // Neutral-hadron opinions: one per input trackster. These COMPETE with the photon
  // hypotheses on the same energy in the same arbitration tier, so a K0L / neutron
  // shower with an EM-rich front is defended instead of being claimed as a photon
  // unopposed. The input tracksters carry no PID here (ticlTracksterLinks runs no
  // inference and merging zeroes the probabilities), so the score is left at 0 and
  // the producer re-scores the whole neutral tier from the PID inference it runs on
  // the hypothesis tracksters before arbitrating.
  const auto &tracksters = input.tracksters;
  for (unsigned iTs = 0; iTs < tracksters.size(); ++iTs) {
    const auto &ts = tracksters[iTs];
    if (ts.raw_energy() < 1.f)
      continue;
    Hypothesis h;
    h.type = Hypothesis::Type::NeutralHadron;
    h.score = 0.f;  // re-scored by the producer after inference
    h.tracksterIdx = static_cast<int>(hypothesisTracksters.size());
    hypothesisTracksters.push_back(ts);
    hypotheses.push_back(h);
  }
}

void ChargedHadronInterpretationAlgo::fillPSetDescription(edm::ParameterSetDescription &desc) {
  desc.add<double>("delta_tk_ts_layer1", 0.02);
  desc.add<double>("delta_tk_ts_interface", 0.03);
  desc.add<double>("timing_quality_threshold", 0.5);
  desc.add<double>("energy_overshoot_fraction", 0.2);
  desc.add<double>("energy_overshoot_max", 10.0);
  desc.add<double>("drCut", 0.02);              // max dR for graph edges
  desc.add<double>("tsTsScoreShift", 1.0);      // shift applied to TS-TS edge cost
  desc.add<double>("trackTsScoreShift", 1.0);   // shift applied to track-TS edge cost
  desc.add<double>("tsTsScoreWeight", 1.0);     // scale for TS-TS edge cost
  desc.add<double>("trackTsScoreWeight", 1.0);  // scale for track-TS edge cost
  desc.add<std::string>("onnxTrackModel", "");
  desc.add<std::string>("onnxTracksterModel", "");
  TICLInterpretationAlgoBase<reco::Track>::fillPSetDescription(desc);
}
