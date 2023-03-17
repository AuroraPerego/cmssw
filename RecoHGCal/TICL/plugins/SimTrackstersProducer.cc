// Author: Felice Pantaleo, Leonardo Cristella - felice.pantaleo@cern.ch, leonardo.cristella@cern.ch
// Date: 09/2021

// user include files

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "SimDataFormats/Associations/interface/LayerClusterToSimClusterAssociator.h"
#include "SimDataFormats/Associations/interface/LayerClusterToCaloParticleAssociator.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "SimDataFormats/Track/interface/UniqueSimTrackId.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

#include "RecoHGCal/TICL/interface/commons.h"

#include "TrackstersPCA.h"
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>

using namespace ticl;

class SimTrackstersProducer : public edm::stream::EDProducer<> {
public:
  explicit SimTrackstersProducer(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::Event&, const edm::EventSetup&) override;
  void addTrackster(const int index,
                    const std::vector<std::pair<edm::Ref<reco::CaloClusterCollection>, std::pair<float, float>>>& lcVec,
                    const std::vector<float>& inputClusterMask,
                    const float fractionCut_,
                    const float energy,
                    const int pdgId,
                    const int charge,
                    float time,
                    const edm::ProductID seed,
                    const Trackster::IterationIndex iter,
                    std::vector<float>& output_mask,
                    std::vector<Trackster>& result);

private:
  std::string detector_;
  const bool doNose_ = false;
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;
  const edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>> clustersTime_token_;
  const edm::EDGetTokenT<std::vector<float>> filtered_layerclusters_mask_token_;

  const edm::EDGetTokenT<std::vector<SimCluster>> simclusters_token_;
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloparticles_token_;

  const edm::EDGetTokenT<hgcal::SimToRecoCollectionWithSimClusters> associatorMapSimClusterToReco_token_;
  const edm::EDGetTokenT<hgcal::SimToRecoCollection> associatorMapCaloParticleToReco_token_;
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geom_token_;
  hgcal::RecHitTools rhtools_;
  const double fractionCut_;
  const edm::EDGetTokenT<std::vector<TrackingParticle>> trackingParticleToken_;
  const edm::EDGetTokenT<std::vector<SimVertex>> simVerticesToken_;

  const edm::EDGetTokenT<std::vector<reco::Track>> recoTracksToken_;
  const StringCutObjectSelector<reco::Track> cutTk_;

  const edm::EDGetTokenT<reco::SimToRecoCollection> associatormapStRsToken_;
  const edm::EDGetTokenT<reco::RecoToSimCollection> associatormapRtSsToken_;
  const edm::EDGetTokenT<SimTrackToTPMap> associationSimTrackToTPToken_;
  
  //timing of RecoTracks from MTD
  const edm::EDGetTokenT<edm::ValueMap<float>> MTDTimeToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> MTDTimeErrorToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> tmtdToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> sigmatmtdToken_;
};
DEFINE_FWK_MODULE(SimTrackstersProducer);

SimTrackstersProducer::SimTrackstersProducer(const edm::ParameterSet& ps)
    : detector_(ps.getParameter<std::string>("detector")),
      doNose_(detector_ == "HFNose"),
      clusters_token_(consumes(ps.getParameter<edm::InputTag>("layer_clusters"))),
      clustersTime_token_(consumes(ps.getParameter<edm::InputTag>("time_layerclusters"))),
      filtered_layerclusters_mask_token_(consumes(ps.getParameter<edm::InputTag>("filtered_mask"))),
      simclusters_token_(consumes(ps.getParameter<edm::InputTag>("simclusters"))),
      caloparticles_token_(consumes(ps.getParameter<edm::InputTag>("caloparticles"))),
      associatorMapSimClusterToReco_token_(
          consumes(ps.getParameter<edm::InputTag>("layerClusterSimClusterAssociator"))),
      associatorMapCaloParticleToReco_token_(
          consumes(ps.getParameter<edm::InputTag>("layerClusterCaloParticleAssociator"))),
      geom_token_(esConsumes()),
      fractionCut_(ps.getParameter<double>("fractionCut")),
      trackingParticleToken_(
          consumes<std::vector<TrackingParticle>>(ps.getParameter<edm::InputTag>("trackingParticles"))),
      simVerticesToken_(
          consumes<std::vector<SimVertex>>(ps.getParameter<edm::InputTag>("simVertices"))),
      recoTracksToken_(consumes<std::vector<reco::Track>>(ps.getParameter<edm::InputTag>("recoTracks"))),
      cutTk_(ps.getParameter<std::string>("cutTk")),
      associatormapStRsToken_(consumes(ps.getParameter<edm::InputTag>("tpToTrack"))),
      associationSimTrackToTPToken_(consumes(ps.getParameter<edm::InputTag>("simTrackToTPMap"))),
      MTDTimeToken_(consumes<edm::ValueMap<float>>(ps.getParameter<edm::InputTag>("trackTimeValueMap"))),
      MTDTimeErrorToken_(consumes<edm::ValueMap<float>>(ps.getParameter<edm::InputTag>("trackTimeErrorMap"))),
      tmtdToken_(consumes<edm::ValueMap<float>>(ps.getParameter<edm::InputTag>("tmtd"))),
      sigmatmtdToken_(consumes<edm::ValueMap<float>>(ps.getParameter<edm::InputTag>("sigmatmtd"))) {
  produces<TracksterCollection>();
  produces<std::vector<float>>();
  produces<TracksterCollection>("fromCPs");
  produces<std::vector<float>>("fromCPs");
  produces<std::map<uint, std::vector<uint>>>();
  produces<std::vector<TICLCandidate>>();
}

void SimTrackstersProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("detector", "HGCAL");
  desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalLayerClusters"));
  desc.add<edm::InputTag>("time_layerclusters", edm::InputTag("hgcalLayerClusters", "timeLayerCluster"));
  desc.add<edm::InputTag>("filtered_mask", edm::InputTag("filteredLayerClustersSimTracksters", "ticlSimTracksters"));
  desc.add<edm::InputTag>("simclusters", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("caloparticles", edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("layerClusterSimClusterAssociator",
                          edm::InputTag("layerClusterSimClusterAssociationProducer"));
  desc.add<edm::InputTag>("layerClusterCaloParticleAssociator",
                          edm::InputTag("layerClusterCaloParticleAssociationProducer"));
  desc.add<edm::InputTag>("recoTracks", edm::InputTag("generalTracks"));
  desc.add<std::string>("cutTk",
                        "1.48 < abs(eta) < 3.0 && pt > 1. && quality(\"highPurity\") && "
                        "hitPattern().numberOfLostHits(\"MISSING_OUTER_HITS\") < 5");
  desc.add<edm::InputTag>("tpToTrack", edm::InputTag("trackingParticleRecoTrackAsssociation"));

  desc.add<edm::InputTag>("trackingParticles", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("simVertices", edm::InputTag("g4SimHits"));

  desc.add<edm::InputTag>("simTrackToTPMap", edm::InputTag("simHitTPAssocProducer", "simTrackToTP"));
  desc.add<edm::InputTag>("trackTimeValueMap", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("trackTimeErrorMap", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("tmtd", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("sigmatmtd", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<double>("fractionCut", 0.);

  descriptions.addWithDefaultLabel(desc);
}

void SimTrackstersProducer::addTrackster(
    const int index,
    const std::vector<std::pair<edm::Ref<reco::CaloClusterCollection>, std::pair<float, float>>>& lcVec,
    const std::vector<float>& inputClusterMask,
    const float fractionCut_,
    const float energy,
    const int pdgId,
    const int charge,
    const float time,
    const edm::ProductID seed,
    const Trackster::IterationIndex iter,
    std::vector<float>& output_mask,
    std::vector<Trackster>& result) {
  if (lcVec.empty())
    return;

  Trackster tmpTrackster;
  tmpTrackster.vertices().reserve(lcVec.size());
  tmpTrackster.vertex_multiplicity().reserve(lcVec.size());
  for (auto const& [lc, energyScorePair] : lcVec) {
    if (inputClusterMask[lc.index()] > 0) {
      double fraction = energyScorePair.first / lc->energy();
      if (fraction < fractionCut_)
        continue;
      tmpTrackster.vertices().push_back(lc.index());
      output_mask[lc.index()] -= fraction;
      tmpTrackster.vertex_multiplicity().push_back(1. / fraction);
    }
  }

  tmpTrackster.setIdProbability(tracksterParticleTypeFromPdgId(pdgId, charge), 1.f);
  tmpTrackster.setIteration(iter);
  tmpTrackster.setSeed(seed, index);
  tmpTrackster.setBoundaryTimeAndError(time, 0.f);
  result.emplace_back(tmpTrackster);
}

void SimTrackstersProducer::produce(edm::Event& evt, const edm::EventSetup& es) {
  auto result = std::make_unique<TracksterCollection>();
  auto output_mask = std::make_unique<std::vector<float>>();
  auto result_fromCP = std::make_unique<TracksterCollection>();
  auto output_mask_fromCP = std::make_unique<std::vector<float>>();
  auto cpToSc_SimTrackstersMap = std::make_unique<std::map<uint, std::vector<uint>>>();
  auto result_ticlCandidates = std::make_unique<std::vector<TICLCandidate>>();

  const auto& layerClusters = evt.get(clusters_token_);
  const auto& layerClustersTimes = evt.get(clustersTime_token_);
  const auto& inputClusterMask = evt.get(filtered_layerclusters_mask_token_);
  output_mask->resize(layerClusters.size(), 1.f);
  output_mask_fromCP->resize(layerClusters.size(), 1.f);

  const auto& simclusters = evt.get(simclusters_token_);
  edm::Handle<std::vector<CaloParticle>> caloParticles_h;
  evt.getByToken(caloparticles_token_, caloParticles_h);
  const auto& caloparticles = *caloParticles_h;

  const auto& simClustersToRecoColl = evt.get(associatorMapSimClusterToReco_token_);
  const auto& caloParticlesToRecoColl = evt.get(associatorMapCaloParticleToReco_token_);
  const auto& simVertices = evt.get(simVerticesToken_);

  edm::Handle<std::vector<TrackingParticle>> trackingParticles_h;
  evt.getByToken(trackingParticleToken_, trackingParticles_h);
  edm::Handle<std::vector<reco::Track>> recoTracks_h;
  evt.getByToken(recoTracksToken_, recoTracks_h);
  const auto& TPtoRecoTrackMap = evt.get(associatormapStRsToken_);
  const auto& simTrackToTPMap = evt.get(associationSimTrackToTPToken_);
  const auto& recoTracks = *recoTracks_h;
  
  edm::Handle<edm::ValueMap<float>> MTDTime_h, MTDTimeError_h;
  evt.getByToken(MTDTimeToken_, MTDTime_h);
  evt.getByToken(MTDTimeErrorToken_, MTDTimeError_h);

  edm::Handle<edm::ValueMap<float>> tMtd_h, tMtdError_h;
  evt.getByToken(tmtdToken_, tMtd_h);
  evt.getByToken(sigmatmtdToken_, tMtdError_h);

  const auto& geom = es.getData(geom_token_);
  rhtools_.setGeometry(geom);
  const auto num_simclusters = simclusters.size();
  result->reserve(num_simclusters);  // Conservative size, will call shrink_to_fit later
  const auto num_caloparticles = caloparticles.size();
  result_fromCP->reserve(num_caloparticles);
  std::map<uint, uint> SimClusterToCaloParticleMap;
  for (const auto& [key, lcVec] : caloParticlesToRecoColl) {
    auto const& cp = *(key);
    auto cpIndex = &cp - &caloparticles[0];
    for (const auto& scRef : cp.simClusters()) {
      auto const& sc = *(scRef);
      auto const scIndex = &sc - &simclusters[0];
      SimClusterToCaloParticleMap[scIndex] = cpIndex;
    }
    auto regr_energy = cp.energy();
    std::vector<uint> scSimTracksterIdx;
    scSimTracksterIdx.reserve(cp.simClusters().size());

    // Create a Trackster from the object entering HGCal
    if (cp.g4Tracks()[0].crossedBoundary()) {
      regr_energy = cp.g4Tracks()[0].getMomentumAtBoundary().energy();
      float time = cp.g4Tracks()[0].getPositionAtBoundary().t();
      addTrackster(cpIndex,
                   lcVec,
                   inputClusterMask,
                   fractionCut_,
                   regr_energy,
                   cp.pdgId(),
                   cp.charge(),
                   time,
                   key.id(),
                   ticl::Trackster::SIM,
                   *output_mask,
                   *result);
    } else {
      for (const auto& scRef : cp.simClusters()) {
        const auto& it = simClustersToRecoColl.find(scRef);
        if (it == simClustersToRecoColl.end())
          continue;
        const auto& lcVec = it->val;
        auto const& sc = *(scRef);
        auto const scIndex = &sc - &simclusters[0];

        addTrackster(scIndex,
                     lcVec,
                     inputClusterMask,
                     fractionCut_,
                     sc.g4Tracks()[0].getMomentumAtBoundary().energy(),
                     sc.pdgId(),
                     sc.charge(),
                     sc.g4Tracks()[0].getPositionAtBoundary().t(),
                     scRef.id(),
                     ticl::Trackster::SIM,
                     *output_mask,
                     *result);

        if (result->empty())
          continue;
        const auto index = result->size() - 1;
        if (std::find(scSimTracksterIdx.begin(), scSimTracksterIdx.end(), index) == scSimTracksterIdx.end()) {
          scSimTracksterIdx.emplace_back(index);
        }
      }
      scSimTracksterIdx.shrink_to_fit();
    } 
    float time = simVertices[cp.g4Tracks()[0].vertIndex()].position().t();
    // Create a Trackster from any CP
    addTrackster(cpIndex,
                 lcVec,
                 inputClusterMask,
                 fractionCut_,
                 regr_energy,
                 cp.pdgId(),
                 cp.charge(),
                 time,
                 key.id(),
                 ticl::Trackster::SIM_CP,
                 *output_mask_fromCP,
                 *result_fromCP);

    if (result_fromCP->empty())
      continue;
    const auto index = result_fromCP->size() - 1;
    if (cpToSc_SimTrackstersMap->find(index) == cpToSc_SimTrackstersMap->end()) {
      (*cpToSc_SimTrackstersMap)[index] = scSimTracksterIdx;
    }
  }
  // TODO: remove time computation from PCA calculation and
  //       store time from boundary position in simTracksters
  ticl::assignPCAtoTracksters(
      *result, layerClusters, layerClustersTimes, rhtools_.getPositionLayer(rhtools_.lastLayerEE(doNose_)).z());
  result->shrink_to_fit();
  ticl::assignPCAtoTracksters(
      *result_fromCP, layerClusters, layerClustersTimes, rhtools_.getPositionLayer(rhtools_.lastLayerEE(doNose_)).z());
  result_fromCP->shrink_to_fit();

  auto simTrackToRecoTrack = [&](UniqueSimTrackId simTkId) -> std::pair<int, float> {
    int trackIdx = -1;
    float quality = 0.f;
    auto ipos = simTrackToTPMap.mapping.find(simTkId);
    if (ipos != simTrackToTPMap.mapping.end()) {
      auto jpos = TPtoRecoTrackMap.find((ipos->second));
      if (jpos != TPtoRecoTrackMap.end()) {
        auto& associatedRecoTracks = jpos->val;
        if (!associatedRecoTracks.empty()) {
          // associated reco tracks are sorted by decreasing quality
          if (associatedRecoTracks[0].second > 0.75f) {
            trackIdx = &(*associatedRecoTracks[0].first) - &recoTracks[0];
            quality = associatedRecoTracks[0].second;
          }
        }
      }
    }
    return {trackIdx, quality};
  };

  // Creating the map from TrackingParticle to SimTrackstersFromCP
  auto& simTrackstersFromCP = *result_fromCP;
  for (unsigned int i = 0; i < simTrackstersFromCP.size(); ++i) {
    const auto& simTrack = caloparticles[simTrackstersFromCP[i].seedIndex()].g4Tracks()[0];
    UniqueSimTrackId simTkIds(simTrack.trackId(), simTrack.eventId());
    auto bestAssociatedRecoTrack = simTrackToRecoTrack(simTkIds);
    if (bestAssociatedRecoTrack.first != -1 and bestAssociatedRecoTrack.second > 0.75f) {
      auto trackIndex = bestAssociatedRecoTrack.first;
      simTrackstersFromCP[i].setTrackIdx(trackIndex);
      simTrackstersFromCP[i].setMTDTimeAndError((*MTDTime_h)[edm::Ref<std::vector<reco::Track>>(recoTracks_h, trackIndex)], (*MTDTimeError_h)[edm::Ref<std::vector<reco::Track>>(recoTracks_h, trackIndex)]);
      simTrackstersFromCP[i].settMtdAndError((*tMtd_h)[edm::Ref<std::vector<reco::Track>>(recoTracks_h, trackIndex)], (*tMtdError_h)[edm::Ref<std::vector<reco::Track>>(recoTracks_h, trackIndex)]);
    }
  }

  auto& simTracksters = *result;
  // Creating the map from TrackingParticle to SimTrackster
  std::unordered_map<unsigned int, std::vector<unsigned int>> TPtoSimTracksterMap;
  for (unsigned int i = 0; i < simTracksters.size(); ++i) {
    const auto& simTrack = (simTracksters[i].seedID() == caloParticles_h.id())
                               ? caloparticles[simTracksters[i].seedIndex()].g4Tracks()[0]
                               : simclusters[simTracksters[i].seedIndex()].g4Tracks()[0];
    UniqueSimTrackId simTkIds(simTrack.trackId(), simTrack.eventId());
    auto bestAssociatedRecoTrack = simTrackToRecoTrack(simTkIds);
    if (bestAssociatedRecoTrack.first != -1 and bestAssociatedRecoTrack.second > 0.75f) {
      auto trackIndex = bestAssociatedRecoTrack.first;
      simTracksters[i].setTrackIdx(trackIndex);
//reco::TrackRef trackref(ticl_cand.trackPtr().id(), int(ticl_cand.trackPtr().key()), &evt.productGetter());
      simTracksters[i].setMTDTimeAndError((*MTDTime_h)[edm::Ref<std::vector<reco::Track>>(recoTracks_h, trackIndex)], (*MTDTimeError_h)[edm::Ref<std::vector<reco::Track>>(recoTracks_h, trackIndex)]);
      simTracksters[i].settMtdAndError((*tMtd_h)[edm::Ref<std::vector<reco::Track>>(recoTracks_h, trackIndex)], (*tMtdError_h)[edm::Ref<std::vector<reco::Track>>(recoTracks_h, trackIndex)]);
    }
  }

  edm::OrphanHandle<std::vector<Trackster>> simTracksters_h = evt.put(std::move(result));


  result_ticlCandidates->resize(caloparticles.size());
  for (size_t i =0; i< simTracksters_h->size(); ++i )
  {
    const auto& simTrackster = (*simTracksters_h)[i];
    int cp_index = (simTrackster.seedID() == caloParticles_h.id())
                               ? simTrackster.seedIndex()
                               : SimClusterToCaloParticleMap[simTrackster.seedIndex()];
    auto& cand = (*result_ticlCandidates)[cp_index];
    cand.addTrackster(edm::Ptr<Trackster>(simTracksters_h, i));
    auto trackIndex = (*result_fromCP)[cp_index].trackIdx();
    cand.setTime((*result_fromCP)[cp_index].time());
    cand.setTimeError(0);

    if(trackIndex != -1 && caloparticles[cp_index].charge()!=0)
      cand.setTrackPtr(edm::Ptr<reco::Track>(recoTracks_h, trackIndex));
  }

  for(size_t i = 0; i< result_ticlCandidates->size(); ++i)
  {
    auto& cand = (*result_ticlCandidates)[i];
    const auto& cp = caloparticles[i];
    auto const particleType = tracksterParticleTypeFromPdgId(cp.pdgId(), 1);
    cand.setIdProbability(particleType, 1.f);
    cand.setPdgId(cp.pdgId());
    cand.setCharge(cp.charge());
    float rawEnergy = 0.f;
    float regressedEnergy = 0.f;
    for (const auto& trackster : cand.tracksters()) {
      rawEnergy += trackster->raw_energy();
      regressedEnergy += trackster->regressed_energy();
    }
    cand.setRawEnergy(rawEnergy);

    if(!cand.trackPtr().isNull())
    {
      auto track = cand.trackPtr();
      math::XYZTLorentzVector p4(regressedEnergy * track->momentum().unit().x(),
                                       regressedEnergy * track->momentum().unit().y(),
                                       regressedEnergy * track->momentum().unit().z(),
                                       regressedEnergy);
      cand.setP4(p4);
    }
    else
    {
      const auto& simTracksterFromCP = (*result_fromCP)[i];
      float regressedEnergy = simTracksterFromCP.regressed_energy();
      math::XYZTLorentzVector p4(regressedEnergy * simTracksterFromCP.barycenter().unit().x(),
                                 regressedEnergy * simTracksterFromCP.barycenter().unit().y(),
                                 regressedEnergy * simTracksterFromCP.barycenter().unit().z(),
                                 regressedEnergy);
      cand.setP4(p4);

    }

  }

  evt.put(std::move(result_ticlCandidates));
  evt.put(std::move(output_mask));
  evt.put(std::move(result_fromCP), "fromCPs");
  evt.put(std::move(output_mask_fromCP), "fromCPs");
  evt.put(std::move(cpToSc_SimTrackstersMap));
}
