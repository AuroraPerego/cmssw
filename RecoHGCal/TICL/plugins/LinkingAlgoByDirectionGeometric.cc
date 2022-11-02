#include <cmath>
#include <string>
#include "RecoHGCal/TICL/plugins/LinkingAlgoByDirectionGeometric.h"

#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/Common.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"

#include "TrackstersPCA.h"
using namespace ticl;

LinkingAlgoByDirectionGeometric::LinkingAlgoByDirectionGeometric(const edm::ParameterSet &conf)
    : LinkingAlgoBase(conf),
      del_tk_ts_layer1_(conf.getParameter<double>("delta_tk_ts_layer1")),
      del_tk_ts_int_(conf.getParameter<double>("delta_tk_ts_interface")),
      del_ts_em_had_(conf.getParameter<double>("delta_ts_em_had")),
      del_ts_had_had_(conf.getParameter<double>("delta_ts_had_had")),
      separation_threshold_(conf.getParameter<double>("separation")),
      timing_quality_threshold_(conf.getParameter<double>("track_time_quality_threshold")),
      cutTk_(conf.getParameter<std::string>("cutTk")) {}

LinkingAlgoByDirectionGeometric::~LinkingAlgoByDirectionGeometric() {}

void LinkingAlgoByDirectionGeometric::initialize(const HGCalDDDConstants *hgcons,
                                                 const hgcal::RecHitTools rhtools,
                                                 const edm::ESHandle<MagneticField> bfieldH,
                                                 const edm::ESHandle<Propagator> propH) {
  hgcons_ = hgcons;
  rhtools_ = rhtools;
  buildLayers();

  bfield_ = bfieldH;
  propagator_ = propH;
}

void LinkingAlgoByDirectionGeometric::fillTrackstersTile(const std::vector<Trackster> &tracksters,
                                                         std::array<TICLLayerTile, 2> &tracksterTiles) {
  // needs only the positive Z co-ordinate of the surface to propagate to
  // the correct sign is calculated inside according to the barycenter of trackster

  unsigned int idx = 0;
  for (const auto &t : tracksters) {
    Vector const &baryc = t.barycenter();

    if (baryc.Eta() > 0)
      tracksterTiles[1].fill(baryc.eta(), baryc.phi(), idx);

    else if (baryc.Eta() < 0)
      tracksterTiles[0].fill(baryc.eta(), baryc.phi(), idx);
    idx++;
  }
}

void LinkingAlgoByDirectionGeometric::findTrackstersInWindow(const std::vector<Trackster> &tracksters,
                                                             const std::array<TICLLayerTile, 2> &tracksterTiles,
                                                             const float delta,
                                                             const float separation,
                                                             std::vector<std::vector<unsigned>> &resultCollection,
                                                             bool useMask = false) {
  std::vector<int> mask(tracksters.size(), 0);
  float delta2 = delta * delta;
  float separation2 = separation * separation;
  auto distance = [](float r0, float r1, float phi0, float phi1) {
    auto delta_phi = reco::deltaPhi(phi0, phi1);
    return (r0 - r1) * (r0 - r1) + r1 * r1 * delta_phi * delta_phi;
  };
  unsigned int seedId = 0;
  for (auto const &t : tracksters) {
    auto const &barycenter = t.barycenter();
    float trackster_eta = barycenter.eta();
    float trackster_phi = barycenter.phi();
    float trackster_z = barycenter.z();
    float trackster_r_over_absz =
        sqrt(barycenter.x() * barycenter.x() + barycenter.y() * barycenter.y()) / std::abs(trackster_z);
    auto sideZ = trackster_eta > 0;  //forward or backward region
    const TICLLayerTile &tile = tracksterTiles[sideZ];
    float eta_min = std::max(abs(trackster_eta) - delta, (float)TileConstants::minEta);
    float eta_max = std::min(abs(trackster_eta) + delta, (float)TileConstants::maxEta);
    // get range of bins touched by delta
    std::array<int, 4> search_box =
        tile.searchBoxEtaPhi(eta_min, eta_max, trackster_phi - delta, trackster_phi + delta);
    std::vector<unsigned> in_delta;
    std::vector<float> distances2;
    int trackster_n = 1;
    for (int eta_i = search_box[0]; eta_i <= search_box[1]; ++eta_i) {
      for (int phi_i = search_box[2]; phi_i <= search_box[3]; ++phi_i) {
        const auto &in_tile = tile[tile.globalBin(eta_i, (phi_i % TileConstants::nPhiBins))];
        for (const unsigned &t_i : in_tile) {
          if (t_i != seedId) {
            // calculate actual distances of tracksters to the seed for a more accurate cut
            auto const &tracksterInTile = tracksters[t_i];
            auto const &tracksterInTile_barycenter = tracksterInTile.barycenter();
            float tracksterInTile_eta = tracksterInTile_barycenter.eta();
            float tracksterInTile_phi = tracksterInTile_barycenter.phi();
            float tracksterInTile_z = tracksterInTile_barycenter.z();
            float tracksterInTile_r_over_absz = sqrt(tracksterInTile_barycenter.x() * tracksterInTile_barycenter.x() +
                                                     tracksterInTile_barycenter.y() * tracksterInTile_barycenter.y()) /
                                                std::abs(tracksterInTile_z);
            std::cout << "Trackster Name\tTrackster X\tTrackster Y\tTrackster Z\tTrackster ROverAbsZ\t Trackster R "
                      << std::endl;
            std::cout << "Trackster Ref \t" << barycenter.x() << "\t" << barycenter.y() << "\t" << barycenter.z()
                      << "\t" << trackster_r_over_absz << "\t" << trackster_r_over_absz * trackster_z << std::endl;

            std::cout << "Trackster " << trackster_n << "\t" << tracksterInTile_barycenter.x() << "\t"
                      << tracksterInTile_barycenter.y() << "\t" << tracksterInTile_barycenter.z() << "\t"
                      << tracksterInTile_r_over_absz << "\t" << tracksterInTile_r_over_absz * tracksterInTile_z
                      << std::endl;
            auto sep2 = distance(trackster_r_over_absz * trackster_z,
                                 tracksterInTile_r_over_absz * trackster_z,
                                 trackster_phi,
                                 tracksterInTile_phi);

            trackster_n++;

            if (sep2 < separation2) {
              in_delta.push_back(t_i);
              distances2.push_back(sep2);
            }
          }
        }
      }
    }
    // sort tracksters found in ascending order of their distances from the seed
    std::vector<unsigned> indices(in_delta.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](unsigned i, unsigned j) { return distances2[i] < distances2[j]; });

    // push back sorted tracksters in the result collection
    for (const unsigned &index : indices) {
      const auto &t_i = in_delta[index];
      if (!mask[t_i]) {
        resultCollection[seedId].push_back(t_i);
        if (useMask)
          mask[t_i] = 1;
      }
    }
    seedId++;
  }
  std::cout << "Result Collection size " << resultCollection.size() << std::endl;
}

bool LinkingAlgoByDirectionGeometric::timeAndEnergyCompatible(float &total_raw_energy,
                                                              const reco::Track &track,
                                                              const Trackster &trackster,
                                                              const float &tkT,
                                                              const float &tkTErr,
                                                              const float &tkTimeQual) {
  return true;
}

void LinkingAlgoByDirectionGeometric::recordTrackster(const unsigned ts,  //trackster index
                                                      const std::vector<Trackster> &tracksters,
                                                      const edm::Handle<std::vector<Trackster>> tsH,
                                                      std::vector<unsigned> &ts_mask,
                                                      float &energy_in_candidate,
                                                      TICLCandidate &candidate) {}

void LinkingAlgoByDirectionGeometric::buildLayers() {
  // build disks at HGCal front & EM-Had interface for track propagation

  float zVal = hgcons_->waferZ(1, true);
  std::pair<float, float> rMinMax = hgcons_->rangeR(zVal, true);

  float zVal_interface = rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z();
  std::pair<float, float> rMinMax_interface = hgcons_->rangeR(zVal_interface, true);

  for (int iSide = 0; iSide < 2; ++iSide) {
    float zSide = (iSide == 0) ? (-1. * zVal) : zVal;
    firstDisk_[iSide] =
        std::make_unique<GeomDet>(Disk::build(Disk::PositionType(0, 0, zSide),
                                              Disk::RotationType(),
                                              SimpleDiskBounds(rMinMax.first, rMinMax.second, zSide - 0.5, zSide + 0.5))
                                      .get());

    zSide = (iSide == 0) ? (-1. * zVal_interface) : zVal_interface;
    interfaceDisk_[iSide] = std::make_unique<GeomDet>(
        Disk::build(Disk::PositionType(0, 0, zSide),
                    Disk::RotationType(),
                    SimpleDiskBounds(rMinMax_interface.first, rMinMax_interface.second, zSide - 0.5, zSide + 0.5))
            .get());
  }
}
void LinkingAlgoByDirectionGeometric::dumpLinksFound(std::vector<std::vector<unsigned>> &resultCollection,
                                                     const char *label) const {
  //#ifdef EDM_ML_DEBUG
  //  if (!(LinkingAlgoBase::algo_verbosity_ > VerbosityLevel::Advanced))
  //    return;

  std::cout << "All links found - " << label << "\n";
  std::cout << "(seed can either be a track or trackster depending on the step)\n";
  for (unsigned i = 0; i < resultCollection.size(); ++i) {
    std::cout << "seed " << i << " - tracksters : ";
    const auto &links = resultCollection[i];
    for (unsigned j = 0; j < links.size(); ++j) {
      std::cout << j;
    }
    std::cout << "\n";
  }
  //#endif  // EDM_ML_DEBUG
}

void DFSVisits(std::vector<std::vector<unsigned>> &graph,
               std::vector<bool> &visits,
               Trackster &outTrackster,
               const std::vector<Trackster> &tracksters,
               const std::vector<reco::CaloCluster> &layerClusters,
               int u) {
  visits[u] = true;
  std::cout << "Visiting " << u << std::endl;
  auto updatedSize = outTrackster.vertices().size();
  for (auto j = graph[u].begin(); j != graph[u].end(); j++) {
    if (!visits[*j]) {
      std::cout << "\t Visiting " << *j << std::endl;
      auto const &thisTrackster = tracksters[*j];
      updatedSize += thisTrackster.vertices().size();
      outTrackster.vertices().reserve(updatedSize);
      outTrackster.vertex_multiplicity().reserve(updatedSize);
      std::copy(std::begin(thisTrackster.vertices()),
                std::end(thisTrackster.vertices()),
                std::back_inserter(outTrackster.vertices()));
      std::copy(std::begin(thisTrackster.vertex_multiplicity()),
                std::end(thisTrackster.vertex_multiplicity()),
                std::back_inserter(outTrackster.vertex_multiplicity()));

      DFSVisits(graph, visits, outTrackster, tracksters, layerClusters, *j);
    }
  }
}
// DFS to collect all the links
std::vector<Trackster> DFS(std::vector<std::vector<unsigned>> &graph,
                           const std::vector<Trackster> &tracksters,
                           const std::vector<reco::CaloCluster> &layerClusters) {
  int graphSize = graph.size();
  std::vector<bool> visits(graphSize, false);
  std::vector<Trackster> resultCollection;
  for (int i = 0; i < graphSize; i++) {
    if (!visits[i]) {
      Trackster outTrackster = tracksters[i];
      DFSVisits(graph, visits, outTrackster, tracksters, layerClusters, i);
      resultCollection.push_back(outTrackster);
    }
  }
  resultCollection.shrink_to_fit();
  return resultCollection;
}

void LinkingAlgoByDirectionGeometric::linkTracksters(const edm::Handle<std::vector<reco::Track>> tkH,
                                                     const edm::ValueMap<float> &tkTime,
                                                     const edm::ValueMap<float> &tkTimeErr,
                                                     const edm::ValueMap<float> &tkTimeQual,
                                                     const std::vector<reco::Muon> &muons,
                                                     const edm::Handle<std::vector<Trackster>> tsH,
                                                     const std::vector<reco::CaloCluster> &layerClusters,
                                                     const edm::ValueMap<std::pair<float, float>> layerClustersTimes,
                                                     std::vector<TICLCandidate> &resultLinked,
                                                     std::vector<TICLCandidate> &chargedHadronsFromTk,
                                                     const EnergyRegressionAndIDModel &model) {
  const auto &tracks = *tkH;
  const auto &tracksters = *tsH;

  auto bFieldProd = bfield_.product();
  const Propagator &prop = (*propagator_);

  std::vector<Trackster> resultTrackstersMerged;
  std::array<TICLLayerTile, 2> tracksterPropTiles = {};  // all Tracksters

  // auto isHadron = [&](const Trackster &t) -> bool {
  //   auto boundary_z = rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z();
  //   return (std::abs(t.barycenter().Z()) > boundary_z);
  // };

  if (LinkingAlgoBase::algo_verbosity_ > VerbosityLevel::Advanced)
    LogDebug("LinkingAlgoByDirectionGeometric") << "------- Geometric Linking ------- \n";

  // Propagate tracks
  std::vector<unsigned> candidateTrackIds;
  candidateTrackIds.reserve(tracks.size());
  std::vector<std::pair<Vector, unsigned>> trackPColl;     // propagated track points and index of track in collection
  std::vector<std::pair<Vector, unsigned>> tkPropIntColl;  // tracks propagated to lastLayerEE
  trackPColl.reserve(tracks.size());
  tkPropIntColl.reserve(tracks.size());
  for (unsigned i = 0; i < tracks.size(); ++i) {
    const auto &tk = tracks[i];
    reco::TrackRef trackref = reco::TrackRef(tkH, i);

    // veto tracks associated to muons
    int muId = PFMuonAlgo::muAssocToTrack(trackref, muons);

    if (LinkingAlgoBase::algo_verbosity_ > VerbosityLevel::Advanced)
      LogDebug("LinkingAlgoByDirectionGeometric")
          << "track " << i << " - eta " << tk.eta() << " phi " << tk.phi() << " time " << tkTime[reco::TrackRef(tkH, i)]
          << " time qual " << tkTimeQual[reco::TrackRef(tkH, i)] << "  muid " << muId << "\n";

    if (!cutTk_((tk)) or muId != -1)
      continue;

    // record tracks that can be used to make a ticlcandidate
    candidateTrackIds.push_back(i);

    // don't consider tracks below 2 GeV for linking
    if (std::sqrt(tk.p() * tk.p() + ticl::mpion2) < tkEnergyCut_)
      continue;

    int iSide = int(tk.eta() > 0);
    const auto &fts = trajectoryStateTransform::outerFreeState((tk), bFieldProd);
    // to the HGCal front
    const auto &tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
    if (tsos.isValid()) {
      Vector trackP(tsos.globalPosition().x(), tsos.globalPosition().y(), tsos.globalPosition().z());
      trackPColl.emplace_back(trackP, i);
    }
    // to lastLayerEE
    const auto &tsos_int = prop.propagate(fts, interfaceDisk_[iSide]->surface());
    if (tsos_int.isValid()) {
      Vector trackP(tsos_int.globalPosition().x(), tsos_int.globalPosition().y(), tsos_int.globalPosition().z());
      tkPropIntColl.emplace_back(trackP, i);
    }
  }  // Tracks
  tkPropIntColl.shrink_to_fit();
  trackPColl.shrink_to_fit();
  candidateTrackIds.shrink_to_fit();

  fillTrackstersTile(tracksters, tracksterPropTiles);
  std::vector<std::vector<unsigned>> trackstersResults(tracksters.size());
  findTrackstersInWindow(tracksters, tracksterPropTiles, 0.03, separation_threshold_, trackstersResults);
  dumpLinksFound(trackstersResults, "Tracksters Links");

  std::vector<int> maskTracksters(tracksters.size(), 1);
  // Merge included tracksters
  //Apply DFS here!
  //
  std::cout << " STARTING DFS " << std::endl;
  auto resultCollectionMerged = DFS(trackstersResults, tracksters, layerClusters);
  std::cout << " END DFS " << std::endl;
  std::cout << "Input Tracksters " << trackstersResults.size() << " Output Tracksters " << resultCollectionMerged.size() << std::endl;
  // Find duplicate LCs
    int countTracksters = 0;
    for(auto& outTrackster : resultCollectionMerged){
    auto &orig_vtx = outTrackster.vertices();
    auto vtx_sorted{orig_vtx};
    std::sort(std::begin(vtx_sorted), std::end(vtx_sorted));
    for (unsigned int iLC = 1; iLC < vtx_sorted.size(); ++iLC) {
      if (vtx_sorted[iLC] == vtx_sorted[iLC - 1]) {
        // Clean up duplicate LCs
        const auto lcIdx = vtx_sorted[iLC];
        const auto firstEl = std::find(orig_vtx.begin(), orig_vtx.end(), lcIdx);
        const auto firstPos = std::distance(std::begin(orig_vtx), firstEl);
        auto iDup = std::find(std::next(firstEl), orig_vtx.end(), lcIdx);
        while (iDup != orig_vtx.end()) {
          orig_vtx.erase(iDup);
          outTrackster.vertex_multiplicity().erase(outTrackster.vertex_multiplicity().begin() +
                                                   std::distance(std::begin(orig_vtx), iDup));
          outTrackster.vertex_multiplicity()[firstPos] -= 1;
          iDup = std::find(std::next(firstEl), orig_vtx.end(), lcIdx);
        };
      }
    }

    outTrackster.zeroProbabilities();
    if (!outTrackster.vertices().empty()) {
      resultTrackstersMerged.push_back(outTrackster);
    }
    std::vector<Trackster> vecOut = {outTrackster};
    assignPCAtoTracksters(
        vecOut, layerClusters, layerClustersTimes, rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z());
    std::cout << "Out Trackster " << countTracksters << " Energy " << vecOut[0].raw_energy() << " Position (eta-phi)"
              << vecOut[0].barycenter().eta() << " , " << vecOut[0].barycenter().phi() << " Z "
              << vecOut[0].barycenter().z() << std::endl;
    countTracksters++;
    std::cout << "-----------------------------------------" << std::endl;
  }

  // Create TrackstersMerge with assignPCA
  // Create vector of Linked Tracksters
  // Pass this vector to assignPCA to have the final TrackstersMerged?
  resultTrackstersMerged.shrink_to_fit();
  assignPCAtoTracksters(
      resultTrackstersMerged, layerClusters, layerClustersTimes, rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z());
  model.energyRegressionAndID(layerClusters, resultTrackstersMerged);

}  // linkTracksters
void LinkingAlgoByDirectionGeometric::fillPSetDescription(edm::ParameterSetDescription &desc) {
  desc.add<std::string>("cutTk",
                        "1.48 < abs(eta) < 3.0 && pt > 1. && quality(\"highPurity\") && "
                        "hitPattern().numberOfLostHits(\"MISSING_OUTER_HITS\") < 5");
  desc.add<double>("delta_tk_ts_layer1", 0.02);
  desc.add<double>("delta_tk_ts_interface", 0.03);
  desc.add<double>("delta_ts_em_had", 0.03);
  desc.add<double>("delta_ts_had_had", 0.03);
  desc.add<double>("separation", 4);  //cm
  desc.add<double>("track_time_quality_threshold", 0.5);
  LinkingAlgoBase::fillPSetDescription(desc);
}
