#include "MCFwithNNInterpretationAlgo.h"
#include "RecoHGCal/TICL/interface/MinCostFlow.h"
#include "RecoHGCal/TICL/plugins/TICLInterpretationPluginFactory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace ticl;
using Vector = ticl::Trackster::Vector;

namespace {
  struct RecoCandidateTmp {
    std::vector<int> tracks;
    std::vector<int> tracksters;

    bool operator==(const RecoCandidateTmp& other) const {
      return tracks == other.tracks && tracksters == other.tracksters;
    }
  };

  template <typename T>
  void sortUnique(std::vector<T>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
  }

  bool isSubset(const std::vector<int>& a, const std::vector<int>& b) {
    std::set<int> sa(a.begin(), a.end()), sb(b.begin(), b.end());
    return std::includes(sb.begin(), sb.end(), sa.begin(), sa.end());
  }

  float energyOf(const std::set<int>& ids, const std::vector<float>& tsEnergy) {
    float out = 0.f;
    for (int ts : ids) {
      if (ts >= 0 && ts < static_cast<int>(tsEnergy.size()))
        out += tsEnergy[ts];
    }
    return out;
  }

  float candidateEnergy(const RecoCandidateTmp& c, const std::vector<float>& tsEnergy) {
    float out = 0.f;
    for (int ts : c.tracksters) {
      if (ts >= 0 && ts < static_cast<int>(tsEnergy.size()))
        out += tsEnergy[ts];
    }
    return out;
  }

  float chargedP(const RecoCandidateTmp& c, const reco::TrackCollection& tracks) {
    float out = 0.f;
    for (int trk : c.tracks) {
      if (trk >= 0 && trk < static_cast<int>(tracks.size()))
        out += tracks[trk].p();
    }
    return out;
  }

  std::unordered_map<int, int> computeUsageCount(const std::vector<RecoCandidateTmp>& candidates) {
    std::unordered_map<int, int> usage;
    for (const auto& c : candidates) {
      std::unordered_set<int> uniq(c.tracksters.begin(), c.tracksters.end());
      for (int ts : uniq)
        ++usage[ts];
    }
    return usage;
  }

  bool allUsedOnce(const std::unordered_map<int, int>& usage) {
    for (const auto& [ts, count] : usage) {
      if (count > 1)
        return false;
    }
    return true;
  }

  std::vector<int> duplicatedTracksters(const std::unordered_map<int, int>& usage) {
    std::vector<int> out;
    for (const auto& [ts, count] : usage) {
      if (count > 1)
        out.push_back(ts);
    }
    return out;
  }

  std::vector<RecoCandidateTmp> mergeOverlappingChargedCandidates(const std::vector<RecoCandidateTmp>& chargedCandidates,
                                                                  const std::vector<float>& tsEnergy,
                                                                  const reco::TrackCollection& tracks,
                                                                  const std::unordered_map<int, int>& usageCount,
                                                                  float overlapThreshold = 0.5f) {
    std::vector<RecoCandidateTmp> merged;
    std::unordered_set<int> used;
    const auto dups = duplicatedTracksters(usageCount);
    const std::unordered_set<int> dupSet(dups.begin(), dups.end());

    for (size_t i = 0; i < chargedCandidates.size(); ++i) {
      if (used.count(i))
        continue;
      const auto& ca = chargedCandidates[i];
      if (ca.tracksters.empty()) {
        merged.push_back(ca);
        used.insert(i);
        continue;
      }

      std::set<int> setA(ca.tracksters.begin(), ca.tracksters.end());
      bool touchesDupA = false;
      for (int ts : setA) {
        if (dupSet.count(ts)) {
          touchesDupA = true;
          break;
        }
      }
      if (!touchesDupA) {
        merged.push_back(ca);
        used.insert(i);
        continue;
      }

      const float energyA = candidateEnergy(ca, tsEnergy);
      bool matched = false;

      for (size_t j = i + 1; j < chargedCandidates.size(); ++j) {
        if (used.count(j))
          continue;
        const auto& cb = chargedCandidates[j];
        if (cb.tracksters.empty())
          continue;

        std::set<int> setB(cb.tracksters.begin(), cb.tracksters.end());
        bool touchesDupB = false;
        for (int ts : setB) {
          if (dupSet.count(ts)) {
            touchesDupB = true;
            break;
          }
        }
        if (!touchesDupB)
          continue;

        std::set<int> shared;
        std::set_intersection(setA.begin(), setA.end(), setB.begin(), setB.end(), std::inserter(shared, shared.begin()));
        const float sharedE = energyOf(shared, tsEnergy);
        if (sharedE <= 0.f)
          continue;

        const float energyB = candidateEnergy(cb, tsEnergy);
        const float totalE = energyA + energyB - sharedE;
        const float overlap = sharedE / std::max(totalE, 1e-12f);

        if (overlap > overlapThreshold) {
          RecoCandidateTmp out;
          out.tracks = ca.tracks;
          out.tracks.insert(out.tracks.end(), cb.tracks.begin(), cb.tracks.end());
          out.tracksters = ca.tracksters;
          out.tracksters.insert(out.tracksters.end(), cb.tracksters.begin(), cb.tracksters.end());
          sortUnique(out.tracks);
          sortUnique(out.tracksters);
          merged.push_back(out);
        } else {
          const float pA = chargedP(ca, tracks);
          const float pB = chargedP(cb, tracks);

          std::set<int> onlyA, onlyB;
          std::set_difference(setA.begin(), setA.end(), shared.begin(), shared.end(), std::inserter(onlyA, onlyA.begin()));
          std::set_difference(setB.begin(), setB.end(), shared.begin(), shared.end(), std::inserter(onlyB, onlyB.begin()));

          const float eOnlyA = energyOf(onlyA, tsEnergy);
          const float eOnlyB = energyOf(onlyB, tsEnergy);

          const float diffA = std::abs(energyA - pA) + std::abs(eOnlyB - pB);
          const float diffB = std::abs(eOnlyA - pA) + std::abs(energyB - pB);

          RecoCandidateTmp caNew = ca;
          RecoCandidateTmp cbNew = cb;
          if (diffA <= diffB) {
            caNew.tracksters.assign(setA.begin(), setA.end());
            cbNew.tracksters.assign(onlyB.begin(), onlyB.end());
          } else {
            caNew.tracksters.assign(onlyA.begin(), onlyA.end());
            cbNew.tracksters.assign(setB.begin(), setB.end());
          }
          if (!caNew.tracksters.empty())
            merged.push_back(caNew);
          if (!cbNew.tracksters.empty())
            merged.push_back(cbNew);
        }
        used.insert(i);
        used.insert(j);
        matched = true;
        break;
      }

      if (!matched && !used.count(i)) {
        merged.push_back(ca);
        used.insert(i);
      }
    }

    return merged;
  }

  std::vector<RecoCandidateTmp> mergeMixedCandidates(const std::vector<RecoCandidateTmp>& candidates,
                                                     const std::vector<float>& tsEnergy,
                                                     const reco::TrackCollection& tracks,
                                                     const std::unordered_map<int, int>& usageCount,
                                                     float overlapThreshold = 0.5f) {
    auto isNeutral = [](const RecoCandidateTmp& c) { return c.tracks.empty(); };

    std::vector<RecoCandidateTmp> merged;
    std::unordered_set<int> used;
    const auto dups = duplicatedTracksters(usageCount);
    const std::unordered_set<int> dupSet(dups.begin(), dups.end());

    for (size_t i = 0; i < candidates.size(); ++i) {
      if (used.count(i))
        continue;
      const auto& ca = candidates[i];
      if (ca.tracksters.empty()) {
        used.insert(i);
        continue;
      }

      std::set<int> setA(ca.tracksters.begin(), ca.tracksters.end());
      bool touchesDupA = false;
      for (int ts : setA) {
        if (dupSet.count(ts)) {
          touchesDupA = true;
          break;
        }
      }
      if (!touchesDupA) {
        merged.push_back(ca);
        used.insert(i);
        continue;
      }

      const float energyA = candidateEnergy(ca, tsEnergy);
      bool matched = false;

      for (size_t j = i + 1; j < candidates.size(); ++j) {
        if (used.count(j))
          continue;
        const auto& cb = candidates[j];
        if (cb.tracksters.empty())
          continue;

        std::set<int> setB(cb.tracksters.begin(), cb.tracksters.end());
        bool touchesDupB = false;
        for (int ts : setB) {
          if (dupSet.count(ts)) {
            touchesDupB = true;
            break;
          }
        }
        if (!touchesDupB)
          continue;

        std::set<int> shared;
        std::set_intersection(setA.begin(), setA.end(), setB.begin(), setB.end(), std::inserter(shared, shared.begin()));
        if (shared.empty())
          continue;

        const float sharedE = energyOf(shared, tsEnergy);
        if (sharedE <= 0.f)
          continue;

        const bool aNeu = isNeutral(ca);
        const bool bNeu = isNeutral(cb);
        const float energyB = candidateEnergy(cb, tsEnergy);
        const float totalE = energyA + energyB - sharedE;
        const float overlap = sharedE / std::max(totalE, 1e-12f);

        std::set<int> onlyA, onlyB;
        std::set_difference(setA.begin(), setA.end(), shared.begin(), shared.end(), std::inserter(onlyA, onlyA.begin()));
        std::set_difference(setB.begin(), setB.end(), shared.begin(), shared.end(), std::inserter(onlyB, onlyB.begin()));

        if (aNeu && bNeu) {
          if (overlap > overlapThreshold) {
            RecoCandidateTmp out;
            out.tracksters.assign(setA.begin(), setA.end());
            out.tracksters.insert(out.tracksters.end(), setB.begin(), setB.end());
            sortUnique(out.tracksters);
            merged.push_back(out);
          } else {
            RecoCandidateTmp caNew = ca, cbNew = cb;
            const float eOnlyA = energyOf(onlyA, tsEnergy);
            const float eOnlyB = energyOf(onlyB, tsEnergy);
            if (eOnlyA >= eOnlyB) {
              caNew.tracksters.assign(setA.begin(), setA.end());
              cbNew.tracksters.assign(onlyB.begin(), onlyB.end());
            } else {
              caNew.tracksters.assign(onlyA.begin(), onlyA.end());
              cbNew.tracksters.assign(setB.begin(), setB.end());
            }
            if (!caNew.tracksters.empty())
              merged.push_back(caNew);
            if (!cbNew.tracksters.empty())
              merged.push_back(cbNew);
          }
          used.insert(i);
          used.insert(j);
          matched = true;
          break;
        }

        if (aNeu != bNeu) {
          const RecoCandidateTmp& charged = aNeu ? cb : ca;
          const RecoCandidateTmp& neutral = aNeu ? ca : cb;
          const std::set<int>& setC = aNeu ? setB : setA;
          const std::set<int>& setN = aNeu ? setA : setB;

          std::set<int> onlyC, onlyN;
          std::set_difference(setC.begin(), setC.end(), shared.begin(), shared.end(), std::inserter(onlyC, onlyC.begin()));
          std::set_difference(setN.begin(), setN.end(), shared.begin(), shared.end(), std::inserter(onlyN, onlyN.begin()));

          const float eOnlyC = energyOf(onlyC, tsEnergy);
          const float pC = chargedP(charged, tracks);
          const float diffWithout = std::abs(eOnlyC - pC);
          const float diffWith = std::abs(eOnlyC + sharedE - pC);

          RecoCandidateTmp chargedNew = charged;
          RecoCandidateTmp neutralNew = neutral;
          if (overlap > overlapThreshold || diffWith <= diffWithout) {
            chargedNew.tracksters.assign(setC.begin(), setC.end());
            neutralNew.tracksters.assign(onlyN.begin(), onlyN.end());
          } else {
            chargedNew.tracksters.assign(onlyC.begin(), onlyC.end());
            neutralNew.tracksters.assign(setN.begin(), setN.end());
          }

          if (aNeu) {
            if (!neutralNew.tracksters.empty())
              merged.push_back(neutralNew);
            if (!chargedNew.tracksters.empty())
              merged.push_back(chargedNew);
          } else {
            if (!chargedNew.tracksters.empty())
              merged.push_back(chargedNew);
            if (!neutralNew.tracksters.empty())
              merged.push_back(neutralNew);
          }

          used.insert(i);
          used.insert(j);
          matched = true;
          break;
        }
      }

      if (!matched && !used.count(i)) {
        merged.push_back(ca);
        used.insert(i);
      }
    }

    return merged;
  }
}  // namespace

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
MCFwithNNInterpretationAlgo::MCFwithNNInterpretationAlgo(const edm::ParameterSet& conf,
                                                         TICLONNXGlobalCache const* cache)
    : TICLInterpretationAlgoBase(conf, cache),
      drCut_(conf.getParameter<double>("drCut")),
      tsTsScoreShift_(conf.getParameter<double>("tsTsScoreShift")),
      trackTsScoreShift_(conf.getParameter<double>("trackTsScoreShift")),
      tsTsScoreWeight_(conf.getParameter<double>("tsTsScoreWeight")),
      trackTsScoreWeight_(conf.getParameter<double>("trackTsScoreWeight")),
      neutralPenalty_(conf.getParameter<int>("neutralPenalty")),
      tracksterInit_(conf.getParameter<int>("tracksterInit")),
      trackInit_(conf.getParameter<int>("trackInit")),
      manyPenalty_(conf.getParameter<int>("manyPenalty")),
      inputNames_({"input"}),
      outputNames_({"score"}) {
  const std::string trackModel = conf.getParameter<std::string>("onnxTrackModel");
  const std::string tracksterModel = conf.getParameter<std::string>("onnxTracksterModel");

  if (cache_ != nullptr) {
    onnxSessionTracks_ = cache_->getByModelPathString(trackModel);
    onnxSessionTracksters_ = cache_->getByModelPathString(tracksterModel);
  }
  if (onnxSessionTracks_==nullptr)
	  edm::LogError("MCFwithNNInterpretationAlgo") << "onnxSessionTracks_ is nullptr !!" ;
  if (onnxSessionTracksters_==nullptr)
	  edm::LogError("MCFwithNNInterpretationAlgo") << "onnxSessionTracksters_ is nullptr !!" ;
}

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

void MCFwithNNInterpretationAlgo::makeCandidates(const Inputs& input,
                                                 edm::Handle<MtdHostCollection> inputTimingh,
                                                 std::vector<Trackster>& resultTracksters,
                                                 std::vector<int>& resultCandidate,
                                                 std::vector<bool>& maskedTracksters) {
  const auto& tracksters = input.tracksters;
  const auto& tracks = *(input.tracksHandle);
  const auto& maskTracks = input.maskedTracks;
  const float drCut2 = drCut_ * drCut_;
  auto bFieldProd = bfield_.product();
  const Propagator& prop = *propagator_;

  struct TsInfo {
    unsigned origIdx;
    float eta, phi, energy;
    float x, y, z;
    float time, timeErr;
    int pid;
  };
  struct TrackInfo {
    int origIdx;
    float eta, phi;
    double pt, p;
  };
  struct Edge {
    unsigned int u, v;
    int64_t cost;
  };

  std::array<TICLLayerTile, 2> tracksterPropTiles = {};
  std::vector<std::vector<TsInfo>> tsAllProp(2);
  tsAllProp[0].reserve(tracksters.size());
  tsAllProp[1].reserve(tracksters.size());

  for (unsigned i = 0; i < tracksters.size(); ++i) {
    const auto& t = tracksters[i];
    const auto& baryc = t.barycenter();
    const auto& probs = t.id_probabilities();
    const int pid = static_cast<int>(std::max_element(probs.begin(), probs.end()) - probs.begin());
    const int side = baryc.eta() >= 0.f ? 1 : 0;
    tracksterPropTiles[side].fill(baryc.eta(), baryc.phi(), tsAllProp[side].size());
    tsAllProp[side].push_back({i,
                               baryc.eta(),
                               baryc.phi(),
                               t.raw_energy(),
                               baryc.x(),
                               baryc.y(),
                               baryc.z(),
                               t.time(),
                               t.timeError(),
                               pid});
  }

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

    const int iSide = int(tk.eta() > 0);
    FreeTrajectoryState fts = tk.outerOk() ? trajectoryStateTransform::outerFreeState(tk, bFieldProd)
                                           : trajectoryStateTransform::initialFreeState(tk, bFieldProd);
    const auto& tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
    if (!tsos.isValid())
      continue;

    GlobalPoint pp = tsos.globalPosition();
    const int outSide = pp.eta() >= 0.f ? 1 : 0;
    validTracks[outSide].push_back({static_cast<int>(i), pp.eta(), pp.phi(), tk.pt(), tk.p()});
  }

  for (int side : {0, 1}) {
    const auto& ts = tsAllProp[side];
    const int nTS = static_cast<int>(ts.size());
    if (nTS == 0)
      continue;

    const auto& sideTracks = validTracks[side];
    const int nSideTracks = static_cast<int>(sideTracks.size());
    const bool useMTDTiming = inputTimingh.isValid();
    constexpr float C_CM_PER_NS = 29.9792458f;

    auto findNeighbours = [&](float seedEta, float seedPhi, int sideIdx) {
      bool sideZ = seedEta > 0;
      const TICLLayerTile& tile = tracksterPropTiles[sideZ];
      float etaMin = std::max(std::fabs(seedEta) - drCut_, static_cast<float>(TileConstants::minEta));
      float etaMax = std::min(std::fabs(seedEta) + drCut_, static_cast<float>(TileConstants::maxEta));
      auto searchBox = tile.searchBoxEtaPhi(etaMin, etaMax, seedPhi - drCut_, seedPhi + drCut_);
      std::vector<unsigned> result;
      for (int etaI = searchBox[0]; etaI <= searchBox[1]; ++etaI) {
        for (int phiI = searchBox[2]; phiI <= searchBox[3]; ++phiI) {
          const auto& inTile = tile[tile.globalBin(etaI, phiI % TileConstants::nPhiBins)];
          for (unsigned ti : inTile) {
            float deta = tsAllProp[sideIdx][ti].eta - seedEta;
            float dphi = reco::deltaPhi(tsAllProp[sideIdx][ti].phi, seedPhi);
            if (deta * deta + dphi * dphi < drCut2)
              result.push_back(ti);
          }
        }
      }
      return result;
    };

    std::vector<Edge> trackTsEdges, tsTsEdges;
    constexpr int TRACK_TS_NFEAT = 18;
    constexpr int TS_TS_NFEAT = 21;
    std::vector<std::pair<unsigned, unsigned>> trkTsRaw, tsTsRaw;
    cms::Ort::FloatArrays trkTsFeatsInput(1), tsTsFeatsInput(1);
    auto& trkTsFeats = trkTsFeatsInput[0];
    auto& tsTsFeats = tsTsFeatsInput[0];

    for (int ti = 0; ti < nSideTracks; ++ti) {
      const auto& trk = sideTracks[ti];
      float trkTime = 0.f, trkTimeErr = -1.f;
      float trkMtdX = 0.f, trkMtdY = 0.f, trkMtdZ = 0.f;
      if (useMTDTiming) {
        auto const& tv = inputTimingh->const_view();
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
        if (trkTimeErr > 0.f && tsInf.timeErr > 0.f) {
          float dx = trkMtdX - tsInf.x, dy = trkMtdY - tsInf.y, dz = trkMtdZ - tsInf.z;
          float tof = std::sqrt(dx * dx + dy * dy + dz * dz) / C_CM_PER_NS;
          deltaTime = tsInf.time - trkTime - tof;
        }
        trkTsFeats.insert(trkTsFeats.end(),
                          {static_cast<float>(trk.pt),
                           static_cast<float>(trk.p),
                           trk.eta,
                           std::sin(trk.phi),
                           std::cos(trk.phi),
                           trkTime,
                           trkTimeErr,
                           tsInf.energy,
                           tsInf.eta,
                           std::sin(tsInf.phi),
                           std::cos(tsInf.phi),
                           tsInf.time,
                           tsInf.timeErr,
                           deltaTime,
                           deltaE,
                           deltaEta,
                           deltaPhi,
                           deltaR});
        trkTsRaw.push_back({static_cast<unsigned>(ti), localJ});
      }
    }

    for (unsigned i = 0; i < static_cast<unsigned>(nTS); ++i) {
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
        if (ts1.timeErr > 0.f && ts2.timeErr > 0.f) {
          float dx = ts1.x - ts2.x, dy = ts1.y - ts2.y, dz = ts1.z - ts2.z;
          float tof = std::sqrt(dx * dx + dy * dy + dz * dz) / C_CM_PER_NS;
          deltaTime = std::abs(ts1.time - ts2.time) - tof;
        }
        tsTsFeats.insert(tsTsFeats.end(),
                         {ts1.energy,
                          ts1.eta,
                          std::sin(ts1.phi),
                          std::cos(ts1.phi),
                          ts1.z,
                          ts1.time,
                          ts1.timeErr,
                          ts2.energy,
                          ts2.eta,
                          std::sin(ts2.phi),
                          std::cos(ts2.phi),
                          ts2.z,
                          ts2.time,
                          ts2.timeErr,
                          deltaTime,
                          samePid,
                          deltaE,
                          deltaEta,
                          deltaPhi,
                          deltaR,
                          deltaZ});
        tsTsRaw.push_back({i, localJ});
      }
    }

    cms::Ort::FloatArrays trkTsScores, tsTsScores;
    if (!trkTsRaw.empty())
      onnxSessionTracks_->runInto(inputNames_,
                                  trkTsFeatsInput,
                                  {{static_cast<int64_t>(trkTsRaw.size()), TRACK_TS_NFEAT}},
                                  outputNames_,
                                  trkTsScores);
    if (!tsTsRaw.empty())
      onnxSessionTracksters_->runInto(inputNames_,
                                      tsTsFeatsInput,
                                      {{static_cast<int64_t>(tsTsRaw.size()), TS_TS_NFEAT}},
                                      outputNames_,
                                      tsTsScores);

    trackTsEdges.reserve(trkTsRaw.size());
    for (size_t k = 0; k < trkTsRaw.size(); ++k) {
      int64_t cost = static_cast<int64_t>(-trkTsScores[k][0] * trackTsScoreWeight_ + trackTsScoreShift_);
      trackTsEdges.push_back({trkTsRaw[k].first, trkTsRaw[k].second, cost});
    }
    tsTsEdges.reserve(tsTsRaw.size());
    for (size_t k = 0; k < tsTsRaw.size(); ++k) {
      int64_t cost = static_cast<int64_t>(-tsTsScores[k][0] * tsTsScoreWeight_ + tsTsScoreShift_);
      tsTsEdges.push_back({tsTsRaw[k].first, tsTsRaw[k].second, cost});
    }

    const int SRC = 0;
    const int TRACK_OFFSET = 1;
    const int TS_IN_OFFSET = TRACK_OFFSET + nSideTracks;
    const int TS_OUT_OFFSET = TS_IN_OFFSET + nTS;
    const int SNK = TS_OUT_OFFSET + nTS;
    MinCostFlow mcf(SNK + 1);

    for (int ti = 0; ti < nSideTracks; ++ti)
      mcf.addArc(SRC, TRACK_OFFSET + ti, nTS, tracksterInit_);
    for (int j = 0; j < nTS; ++j)
      mcf.addArc(SRC, TS_IN_OFFSET + j, 1, trackInit_);
    for (const auto& e : trackTsEdges)
      mcf.addArc(TRACK_OFFSET + e.u, TS_IN_OFFSET + e.v, 1, e.cost);
    for (int j = 0; j < nTS; ++j) {
      mcf.addArc(TS_IN_OFFSET + j, TS_OUT_OFFSET + j, 1, 0);
      mcf.addArc(TS_IN_OFFSET + j, TS_OUT_OFFSET + j, 1, manyPenalty_);
    }
    for (const auto& e : tsTsEdges)
      mcf.addArc(TS_OUT_OFFSET + e.u, TS_IN_OFFSET + e.v, 1, e.cost);
    for (int j = 0; j < nTS; ++j)
      mcf.addArc(TS_OUT_OFFSET + j, SNK, 1, neutralPenalty_);
    for (int ti = 0; ti < nSideTracks; ++ti)
      mcf.addArc(TRACK_OFFSET + ti, SNK, 1, neutralPenalty_);

    mcf.setNodeSupply(SRC, nTS);
    mcf.setNodeSupply(SNK, -nTS);
    if (mcf.solve() != MinCostFlow::OPTIMAL) {
      edm::LogWarning("MCFwithNNInterpretationAlgo") << "Min-cost flow did not find an optimal solution";
      continue;
    }

    std::unordered_map<int, std::vector<int>> usedOut;
    for (int arc = 0; arc < mcf.numArcs(); ++arc) {
      if (mcf.flow(arc) > 0)
        usedOut[mcf.tail(arc)].push_back(mcf.head(arc));
    }

    std::function<std::vector<int>(int, std::unordered_set<int>&)> followChain =
        [&](int startNode, std::unordered_set<int>& visited) -> std::vector<int> {
      std::vector<int> tsChain;
      int cur = startNode;
      while (cur != SNK) {
        if (visited.count(cur))
          break;
        visited.insert(cur);
        if (cur >= TS_IN_OFFSET && cur < TS_OUT_OFFSET)
          tsChain.push_back(cur - TS_IN_OFFSET);
        auto it = usedOut.find(cur);
        if (it == usedOut.end() || it->second.empty())
          break;
        auto nexts = it->second;
        sortUnique(nexts);
        if (nexts.size() == 1) {
          cur = nexts[0];
        } else {
          for (int nextNode : nexts) {
            auto branch = followChain(nextNode, visited);
            tsChain.insert(tsChain.end(), branch.begin(), branch.end());
            break;
          }
          return tsChain;
        }
      }
      return tsChain;
    };

    std::vector<float> sideTsEnergy(nTS, 0.f);
    for (int i = 0; i < nTS; ++i)
      sideTsEnergy[i] = ts[i].energy;

    std::vector<RecoCandidateTmp> chargedCandidates;
    auto usageCount = std::unordered_map<int, int>{};

    auto srcIt = usedOut.find(SRC);
    if (srcIt != usedOut.end()) {
      for (int trkNode : srcIt->second) {
        if (trkNode < TRACK_OFFSET || trkNode >= TS_IN_OFFSET)
          continue;
        int ti = trkNode - TRACK_OFFSET;
        int origTrackIdx = sideTracks[ti].origIdx;
        std::vector<int> tsList;
        std::unordered_set<int> visited;
        auto trkIt = usedOut.find(trkNode);
        if (trkIt != usedOut.end()) {
          for (int v : trkIt->second) {
            auto chain = followChain(v, visited);
            tsList.insert(tsList.end(), chain.begin(), chain.end());
          }
        }
        sortUnique(tsList);

        if (tsList.empty()) {
          chargedCandidates.push_back({{origTrackIdx}, {}});
          continue;
        }

        bool absorbed = false;
        for (auto& cand : chargedCandidates) {
          if (!cand.tracksters.empty() && isSubset(tsList, cand.tracksters)) {
            cand.tracks.push_back(origTrackIdx);
            sortUnique(cand.tracks);
            absorbed = true;
            break;
          }
        }
        if (!absorbed) {
          chargedCandidates.push_back({{origTrackIdx}, tsList});
          std::unordered_set<int> uniq(tsList.begin(), tsList.end());
          for (int tsid : uniq)
            ++usageCount[tsid];
        }
      }
    }

    while (!allUsedOnce(usageCount)) {
      auto newCandidates = mergeOverlappingChargedCandidates(chargedCandidates, sideTsEnergy, tracks, usageCount);
      if (newCandidates == chargedCandidates)
        break;
      chargedCandidates = std::move(newCandidates);
      usageCount = computeUsageCount(chargedCandidates);
    }

    std::vector<std::vector<int>> usedTsSets;
    for (const auto& c : chargedCandidates)
      usedTsSets.push_back(c.tracksters);

    std::vector<RecoCandidateTmp> neutralCandidates;
    if (srcIt != usedOut.end()) {
      for (int startNode : srcIt->second) {
        if (startNode < TS_IN_OFFSET || startNode >= TS_OUT_OFFSET)
          continue;
        std::unordered_set<int> visited;
        std::vector<int> tsChain = followChain(startNode, visited);
        sortUnique(tsChain);
        bool skip = false;
        for (const auto& chargedTs : usedTsSets) {
          if (isSubset(tsChain, chargedTs)) {
            skip = true;
            break;
          }
        }
        if (skip || tsChain.empty())
          continue;
        neutralCandidates.push_back({{}, tsChain});
        usedTsSets.push_back(tsChain);
      }
    }

    std::vector<RecoCandidateTmp> candidates = chargedCandidates;
    candidates.insert(candidates.end(), neutralCandidates.begin(), neutralCandidates.end());
    usageCount = computeUsageCount(candidates);
    while (!allUsedOnce(usageCount)) {
      auto newCandidates = mergeMixedCandidates(candidates, sideTsEnergy, tracks, usageCount, 0.5f);
      if (newCandidates == candidates)
        break;
      candidates = std::move(newCandidates);
      usageCount = computeUsageCount(candidates);
    }

    auto pushFinalCandidate = [&](const RecoCandidateTmp& cand) {
      if (cand.tracksters.empty()) {
        for (int trk : cand.tracks)
          resultCandidate[trk] = static_cast<int>(resultTracksters.size());
        return;
      }

      std::vector<int> localTs = cand.tracksters;
      sortUnique(localTs);
      if (localTs.size() == 1) {
        const auto& outTs = tracksters[ts[localTs[0]].origIdx];
        for (int trk : cand.tracks)
          resultCandidate[trk] = static_cast<int>(resultTracksters.size());
        resultTracksters.push_back(outTs);
        return;
      }

      Trackster merged;
      bool isHadron = false;
      for (int idx : localTs) {
        merged.mergeTracksters(tracksters[ts[idx].origIdx]);
        if (tracksters[ts[idx].origIdx].isHadronic())
          isHadron = true;
      }
      merged.setIdProbability(isHadron ? ticl::Trackster::ParticleType::charged_hadron
                                       : ticl::Trackster::ParticleType::electron,
                              1.f);
      for (int trk : cand.tracks)
        resultCandidate[trk] = static_cast<int>(resultTracksters.size());
      resultTracksters.push_back(merged);
    };

    for (const auto& cand : candidates)
      pushFinalCandidate(cand);
  }
}

void MCFwithNNInterpretationAlgo::fillPSetDescription(edm::ParameterSetDescription& desc) {
  desc.add<double>("drCut", 0.02);
  desc.add<double>("tsTsScoreShift", 1.0);
  desc.add<double>("trackTsScoreShift", 1.0);
  desc.add<double>("tsTsScoreWeight", 1.0);
  desc.add<double>("trackTsScoreWeight", 1.0);
  desc.add<int>("neutralPenalty", 1);
  desc.add<int>("tracksterInit", 0);
  desc.add<int>("trackInit", 0);
  desc.add<int>("manyPenalty", 1);
  desc.add<std::string>("onnxTrackModel", "");
  desc.add<std::string>("onnxTracksterModel", "");
  TICLInterpretationAlgoBase::fillPSetDescription(desc);
}

DEFINE_EDM_PLUGIN(TICLGeneralInterpretationPluginFactory,
                  ticl::MCFwithNNInterpretationAlgo,
                  "MCFwithNNInterpretationAlgo");
