// Original Authors:  Philipp Zehetner, Wahid Redjeb

#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include <memory>  // unique_ptr
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/PluginDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/FTLRecHit/interface/FTLCluster.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/TICLGraph.h"
#include "DataFormats/HGCalReco/interface/TICLCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/MTDTrackingRecHit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ForwardDetId/interface/MTDDetId.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/HGCalReco/interface/Common.h"

#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimTrackster.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimTracksterFwd.h"
#include "SimDataFormats/Track/interface/UniqueSimTrackId.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "SimCalorimetry/HGCalAssociatorProducers/interface/AssociatorTools.h"
#include "SimDataFormats/Associations/interface/TracksterToSimTracksterHitLCAssociator.h"
#include "RecoHGCal/TICL/interface/commons.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

class TICLDumper : public edm::one::EDAnalyzer<edm::one::WatchRuns,
                                               edm::one::SharedResources> {
public:
  explicit TICLDumper(const edm::ParameterSet&);
  ~TICLDumper() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  typedef math::XYZVector Vector;
  typedef std::vector<double> Vec;

private:
  void beginJob() override;
  void beginRun(const edm::Run&, const edm::EventSetup&) override;

  void initialize(const HGCalDDDConstants* hgcons,
                  const hgcal::RecHitTools rhtools,
                  const edm::ESHandle<MagneticField> bfieldH,
                  const edm::ESHandle<Propagator> propH);
  void buildLayers();

  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override{};
  void endJob() override;

  // some options
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_token_;
  const edm::EDGetTokenT<std::vector<reco::CaloCluster>> layer_clusters_token_;
  const edm::EDGetTokenT<TICLGraph> ticl_graph_token_;
  const edm::EDGetTokenT<std::vector<TICLCandidate>> ticl_candidates_token_;
  const edm::EDGetTokenT<std::vector<reco::Track>> tracks_token_;
  const edm::EDGetTokenT<std::vector<bool>> tracks_mask_token_;
  const edm::EDGetTokenT<edm::ValueMap<float>> tracks_time_token_;
  const edm::EDGetTokenT<edm::ValueMap<float>> tracks_time_quality_token_;
  const edm::EDGetTokenT<edm::ValueMap<float>> tracks_time_err_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_x_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_y_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_z_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_eta_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_phi_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_px_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_py_token_;
  const edm::EDGetTokenT<std::vector<double>> hgcaltracks_pz_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksters_merged_token_;
  const edm::EDGetTokenT<std::vector<float>> layerClustersLocalDensity_token_;
  const edm::EDGetTokenT<std::vector<float>> layerClustersRadius_token_;
  const edm::EDGetTokenT<edm::ValueMap<std::pair<float, float>>>
      clustersTime_token_;
  const edm::EDGetTokenT<std::vector<int>> tracksterSeeds_token_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometry_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> simTracksters_SC_token_;
  const edm::EDGetTokenT<std::vector<ticl::Trackster>> simTracksters_CP_token_;
  const edm::EDGetTokenT<std::vector<TICLCandidate>> simTICLCandidate_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters>
      tsRecoToSimSC_token_;
  const edm::EDGetTokenT<hgcal::SimToRecoCollectionSimTracksters>
      tsSimToRecoSC_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters>
      tsRecoToSimCP_token_;
  const edm::EDGetTokenT<hgcal::SimToRecoCollectionSimTracksters>
      tsSimToRecoCP_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters>
      MergeRecoToSimSC_token_;
  const edm::EDGetTokenT<hgcal::SimToRecoCollectionSimTracksters>
      MergeSimToRecoSC_token_;
  const edm::EDGetTokenT<hgcal::RecoToSimCollectionSimTracksters>
      MergeRecoToSimCP_token_;
  const edm::EDGetTokenT<hgcal::SimToRecoCollectionSimTracksters>
      MergeSimToRecoCP_token_;
  const edm::EDGetTokenT<std::vector<SimCluster>> simclusters_token_;
  const edm::EDGetTokenT<std::vector<CaloParticle>> caloparticles_token_;
  const edm::EDGetTokenT<MtdSimTracksterCollection> MTDSimTrackstersToken_;
  const edm::EDGetTokenT<MtdSimLayerClusterCollection> MTDSimLayerClustersToken_;
  const edm::EDGetTokenT<reco::SimToRecoCollection> assocTpToTrackToken_;
  const edm::EDGetTokenT<SimTrackToTPMap> associationSimTrackToTPToken_;
  const edm::EDGetTokenT<reco::TrackCollection> mtdTracksToken_;
  const edm::EDGetTokenT<edm::ValueMap<int>> mtdTracksAssocToken_;
  const edm::EDGetTokenT<FTLRecHitCollection> btlRecHitsToken_;
  const edm::EDGetTokenT<FTLRecHitCollection> etlRecHitsToken_;
  const edm::EDGetTokenT<FTLClusterCollection> btlClustersToken_;
  const edm::EDGetTokenT<FTLClusterCollection> etlClustersToken_;

  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometry_token_;
  const std::string detector_;
  const std::string propName_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfield_token_;
  const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagator_token_;
  hgcal::RecHitTools rhtools_;
  edm::ESGetToken<HGCalDDDConstants, IdealGeometryRecord> hdc_token_;
  const HGCalDDDConstants* hgcons_;
  std::unique_ptr<GeomDet> firstDisk_[2];
  std::unique_ptr<GeomDet> interfaceDisk_[2];
  edm::ESHandle<MagneticField> bfield_;
  edm::ESHandle<Propagator> propagator_;

  // Output tree
  TTree* tree_;

  void clearVariables();

  unsigned int event_index;

  // Variables for branches
  unsigned int ev_event_;
  unsigned int ntracksters_;
  unsigned int nclusters_;
  unsigned int stsSC_ntracksters_;
  unsigned int stsCP_ntracksters_;
  unsigned int tracksters_merged_ntracksters_;
  size_t nsimTrackstersSC;
  size_t nsimTrackstersCP;

  std::vector<float> trackster_time;
  std::vector<float> trackster_timeError;
  std::vector<float> trackster_regressed_energy;
  std::vector<float> trackster_raw_energy;
  std::vector<float> trackster_raw_em_energy;
  std::vector<float> trackster_raw_pt;
  std::vector<float> trackster_raw_em_pt;
  std::vector<float> trackster_barycenter_x;
  std::vector<float> trackster_barycenter_y;
  std::vector<float> trackster_barycenter_z;
  std::vector<float> trackster_EV1;
  std::vector<float> trackster_EV2;
  std::vector<float> trackster_EV3;
  std::vector<float> trackster_eVector0_x;
  std::vector<float> trackster_eVector0_y;
  std::vector<float> trackster_eVector0_z;
  std::vector<float> trackster_sigmaPCA1;
  std::vector<float> trackster_sigmaPCA2;
  std::vector<float> trackster_sigmaPCA3;
  std::vector<float> trackster_barycenter_eta;
  std::vector<float> trackster_barycenter_phi;
  std::vector<std::vector<float>> trackster_id_probabilities;
  std::vector<std::vector<uint32_t>> trackster_vertices_indexes;
  std::vector<std::vector<float>> trackster_vertices_x;
  std::vector<std::vector<float>> trackster_vertices_y;
  std::vector<std::vector<float>> trackster_vertices_z;
  std::vector<std::vector<float>> trackster_vertices_time;
  std::vector<std::vector<float>> trackster_vertices_timeErr;
  std::vector<std::vector<float>> trackster_vertices_energy;
  std::vector<std::vector<float>> trackster_vertices_correctedEnergy;
  std::vector<std::vector<float>> trackster_vertices_correctedEnergyUncertainty;
  std::vector<std::vector<float>> trackster_vertices_multiplicity;

  std::vector<float> stsSC_trackster_time;
  std::vector<float> stsSC_trackster_timeError;
  std::vector<float> stsSC_trackster_t0Mtd;
  std::vector<float> stsSC_trackster_t0MtdError;
  std::vector<float> stsSC_trackster_tMtd;
  std::vector<float> stsSC_trackster_tMtdError;
  std::vector<float> stsSC_trackster_speedMtd;
  std::vector<float> stsSC_trackster_tMtdSim;
  std::vector<GlobalPoint> stsSC_trackster_tMtdPos;
  std::vector<float> stsSC_trackster_regressed_energy;
  std::vector<float> stsSC_trackster_raw_energy;
  std::vector<float> stsSC_trackster_raw_em_energy;
  std::vector<float> stsSC_trackster_raw_pt;
  std::vector<float> stsSC_trackster_raw_em_pt;
  std::vector<float> stsSC_trackster_barycenter_x;
  std::vector<float> stsSC_trackster_barycenter_y;
  std::vector<float> stsSC_trackster_barycenter_z;
  std::vector<float> stsSC_trackster_barycenter_eta;
  std::vector<float> stsSC_trackster_barycenter_phi;
  std::vector<float> stsSC_trackster_EV1;
  std::vector<float> stsSC_trackster_EV2;
  std::vector<float> stsSC_trackster_EV3;
  std::vector<float> stsSC_trackster_eVector0_x;
  std::vector<float> stsSC_trackster_eVector0_y;
  std::vector<float> stsSC_trackster_eVector0_z;
  std::vector<float> stsSC_trackster_sigmaPCA1;
  std::vector<float> stsSC_trackster_sigmaPCA2;
  std::vector<float> stsSC_trackster_sigmaPCA3;
  std::vector<int> stsSC_pdgID;
  std::vector<int> stsSC_trackIdx;
  std::vector<float> stsSC_trackTime;
  std::vector<float> stsSC_boundaryX;
  std::vector<float> stsSC_boundaryY;
  std::vector<float> stsSC_boundaryZ;
  std::vector<float> stsSC_boundaryEta;
  std::vector<float> stsSC_boundaryPhi;
  std::vector<float> stsSC_boundaryPx;
  std::vector<float> stsSC_boundaryPy;
  std::vector<float> stsSC_boundaryPz;
  std::vector<float> stsSC_trackIdx;
  std::vector<float> stsSC_track_boundaryX;
  std::vector<float> stsSC_track_boundaryY;
  std::vector<float> stsSC_track_boundaryZ;
  std::vector<float> stsSC_track_boundaryEta;
  std::vector<float> stsSC_track_boundaryPhi;
  std::vector<float> stsSC_track_boundaryPx;
  std::vector<float> stsSC_track_boundaryPy;
  std::vector<float> stsSC_track_boundaryPz;
  std::vector<std::vector<float>> stsSC_trackster_id_probabilities;
  std::vector<std::vector<uint32_t>> stsSC_trackster_vertices_indexes;
  std::vector<std::vector<float>> stsSC_trackster_vertices_x;
  std::vector<std::vector<float>> stsSC_trackster_vertices_y;
  std::vector<std::vector<float>> stsSC_trackster_vertices_z;
  std::vector<std::vector<float>> stsSC_trackster_vertices_time;
  std::vector<std::vector<float>> stsSC_trackster_vertices_timeErr;
  std::vector<std::vector<float>> stsSC_trackster_vertices_energy;
  std::vector<std::vector<float>> stsSC_trackster_vertices_correctedEnergy;
  std::vector<std::vector<float>>
      stsSC_trackster_vertices_correctedEnergyUncertainty;
  std::vector<std::vector<float>> stsSC_trackster_vertices_multiplicity;
  std::vector<float> stsCP_trackster_time;
  std::vector<float> stsCP_trackster_timeError;
  std::vector<float> stsCP_trackster_t0Mtd;
  std::vector<float> stsCP_trackster_t0MtdError;
  std::vector<float> stsCP_trackster_tMtd;
  std::vector<float> stsCP_trackster_tMtdError;
  std::vector<float> stsCP_trackster_speedMtd;
  std::vector<float> stsCP_trackster_tMtdSim;
  std::vector<GlobalPoint> stsCP_trackster_tMtdPos;
  std::vector<float> stsCP_trackster_regressed_energy;
  std::vector<float> stsCP_trackster_raw_energy;
  std::vector<float> stsCP_trackster_raw_em_energy;
  std::vector<float> stsCP_trackster_raw_pt;
  std::vector<float> stsCP_trackster_raw_em_pt;
  std::vector<float> stsCP_trackster_barycenter_x;
  std::vector<float> stsCP_trackster_barycenter_y;
  std::vector<float> stsCP_trackster_barycenter_z;
  std::vector<float> stsCP_trackster_barycenter_eta;
  std::vector<float> stsCP_trackster_barycenter_phi;
  std::vector<float> stsCP_trackster_EV1;
  std::vector<float> stsCP_trackster_EV2;
  std::vector<float> stsCP_trackster_EV3;
  std::vector<float> stsCP_trackster_eVector0_x;
  std::vector<float> stsCP_trackster_eVector0_y;
  std::vector<float> stsCP_trackster_eVector0_z;
  std::vector<float> stsCP_trackster_sigmaPCA1;
  std::vector<float> stsCP_trackster_sigmaPCA2;
  std::vector<float> stsCP_trackster_sigmaPCA3;
  std::vector<int> stsCP_pdgID;
  std::vector<int> stsCP_trackIdx;
  std::vector<float> stsCP_trackTime;
  std::vector<float> stsCP_boundaryX;
  std::vector<float> stsCP_boundaryY;
  std::vector<float> stsCP_boundaryZ;
  std::vector<float> stsCP_boundaryEta;
  std::vector<float> stsCP_boundaryPhi;
  std::vector<float> stsCP_boundaryPx;
  std::vector<float> stsCP_boundaryPy;
  std::vector<float> stsCP_boundaryPz;
  std::vector<float> stsCP_trackIdx;
  std::vector<float> stsCP_track_boundaryX;
  std::vector<float> stsCP_track_boundaryY;
  std::vector<float> stsCP_track_boundaryZ;
  std::vector<float> stsCP_track_boundaryEta;
  std::vector<float> stsCP_track_boundaryPhi;
  std::vector<float> stsCP_track_boundaryPx;
  std::vector<float> stsCP_track_boundaryPy;
  std::vector<float> stsCP_track_boundaryPz;
  std::vector<std::vector<float>> stsCP_trackster_id_probabilities;
  std::vector<std::vector<uint32_t>> stsCP_trackster_vertices_indexes;
  std::vector<std::vector<float>> stsCP_trackster_vertices_x;
  std::vector<std::vector<float>> stsCP_trackster_vertices_y;
  std::vector<std::vector<float>> stsCP_trackster_vertices_z;
  std::vector<std::vector<float>> stsCP_trackster_vertices_time;
  std::vector<std::vector<float>> stsCP_trackster_vertices_timeErr;
  std::vector<std::vector<float>> stsCP_trackster_vertices_energy;
  std::vector<std::vector<float>> stsCP_trackster_vertices_correctedEnergy;
  std::vector<std::vector<float>>
      stsCP_trackster_vertices_correctedEnergyUncertainty;
  std::vector<std::vector<float>> stsCP_trackster_vertices_multiplicity;

  // from TICLGraph
  std::vector<std::vector<uint32_t>> node_linked_inners;
  std::vector<std::vector<float>> node_linked_scores;
  std::vector<std::vector<uint32_t>> node_linked_outers;
  std::vector<bool> isRootTrackster;

  std::vector<float> simTICLCandidate_raw_energy;
  std::vector<float> simTICLCandidate_regressed_energy;
  std::vector<std::vector<int>> simTICLCandidate_simTracksterCPIndex;
  std::vector<float> simTICLCandidate_boundaryX;
  std::vector<float> simTICLCandidate_boundaryY;
  std::vector<float> simTICLCandidate_boundaryZ;
  std::vector<float> simTICLCandidate_boundaryPx;
  std::vector<float> simTICLCandidate_boundaryPy;
  std::vector<float> simTICLCandidate_boundaryPz;
  std::vector<float> simTICLCandidate_trackTime;
  std::vector<float> simTICLCandidate_trackBeta;
  std::vector<float> simTICLCandidate_caloParticleMass;
  std::vector<int> simTICLCandidate_pdgId;
  std::vector<int> simTICLCandidate_charge;
  std::vector<int> simTICLCandidate_track_in_candidate;

  // from TICLCandidate, product of linking
  size_t nCandidates;
  std::vector<int> candidate_charge;
  std::vector<int> candidate_pdgId;
  std::vector<float> candidate_energy;
  std::vector<double> candidate_px;
  std::vector<double> candidate_py;
  std::vector<double> candidate_pz;
  std::vector<float> candidate_time;
  std::vector<float> candidate_time_err;
  std::vector<std::vector<float>> candidate_id_probabilities;
  std::vector<std::vector<uint32_t>> tracksters_in_candidate;
  std::vector<int> track_in_candidate;

  // merged tracksters
  size_t nTrackstersMerged;
  std::vector<float> tracksters_merged_time;
  std::vector<float> tracksters_merged_timeError;
  std::vector<float> tracksters_merged_t0Mtd;
  std::vector<float> tracksters_merged_t0MtdError;
  std::vector<float> tracksters_merged_tMtd;
  std::vector<float> tracksters_merged_tMtdError;
  std::vector<float> tracksters_merged_speedMtd;
  std::vector<GlobalPoint> tracksters_merged_tMtdPos;
  std::vector<float> tracksters_merged_regressed_energy;
  std::vector<float> tracksters_merged_raw_energy;
  std::vector<float> tracksters_merged_raw_em_energy;
  std::vector<float> tracksters_merged_raw_pt;
  std::vector<float> tracksters_merged_raw_em_pt;
  std::vector<float> tracksters_merged_barycenter_x;
  std::vector<float> tracksters_merged_barycenter_y;
  std::vector<float> tracksters_merged_barycenter_z;
  std::vector<float> tracksters_merged_barycenter_eta;
  std::vector<float> tracksters_merged_barycenter_phi;
  std::vector<float> tracksters_merged_EV1;
  std::vector<float> tracksters_merged_EV2;
  std::vector<float> tracksters_merged_EV3;
  std::vector<float> tracksters_merged_eVector0_x;
  std::vector<float> tracksters_merged_eVector0_y;
  std::vector<float> tracksters_merged_eVector0_z;
  std::vector<float> tracksters_merged_sigmaPCA1;
  std::vector<float> tracksters_merged_sigmaPCA2;
  std::vector<float> tracksters_merged_sigmaPCA3;
  std::vector<std::vector<uint32_t>> tracksters_merged_vertices_indexes;
  std::vector<std::vector<float>> tracksters_merged_vertices_x;
  std::vector<std::vector<float>> tracksters_merged_vertices_y;
  std::vector<std::vector<float>> tracksters_merged_vertices_z;
  std::vector<std::vector<float>> tracksters_merged_vertices_time;
  std::vector<std::vector<float>> tracksters_merged_vertices_timeErr;
  std::vector<std::vector<float>> tracksters_merged_vertices_energy;
  std::vector<std::vector<float>> tracksters_merged_vertices_correctedEnergy;
  std::vector<std::vector<float>>
      tracksters_merged_vertices_correctedEnergyUncertainty;
  std::vector<std::vector<float>> tracksters_merged_vertices_multiplicity;
  std::vector<std::vector<float>> tracksters_merged_id_probabilities;

  // associations
  std::vector<std::vector<uint32_t>> trackstersCLUE3D_recoToSim_SC;
  std::vector<std::vector<float>> trackstersCLUE3D_recoToSim_SC_score;
  std::vector<std::vector<float>> trackstersCLUE3D_recoToSim_SC_sharedE;
  std::vector<std::vector<uint32_t>> trackstersCLUE3D_simToReco_SC;
  std::vector<std::vector<float>> trackstersCLUE3D_simToReco_SC_score;
  std::vector<std::vector<float>> trackstersCLUE3D_simToReco_SC_sharedE;

  std::vector<std::vector<uint32_t>> trackstersCLUE3D_recoToSim_CP;
  std::vector<std::vector<float>> trackstersCLUE3D_recoToSim_CP_score;
  std::vector<std::vector<float>> trackstersCLUE3D_recoToSim_CP_sharedE;
  std::vector<std::vector<uint32_t>> trackstersCLUE3D_simToReco_CP;
  std::vector<std::vector<float>> trackstersCLUE3D_simToReco_CP_score;
  std::vector<std::vector<float>> trackstersCLUE3D_simToReco_CP_sharedE;

  std::vector<std::vector<uint32_t>> MergeTracksters_recoToSim_SC;
  std::vector<std::vector<float>> MergeTracksters_recoToSim_SC_score;
  std::vector<std::vector<float>> MergeTracksters_recoToSim_SC_sharedE;
  std::vector<std::vector<uint32_t>> MergeTracksters_simToReco_SC;
  std::vector<std::vector<float>> MergeTracksters_simToReco_SC_score;
  std::vector<std::vector<float>> MergeTracksters_simToReco_SC_sharedE;

  std::vector<std::vector<uint32_t>> MergeTracksters_recoToSim_CP;
  std::vector<std::vector<float>> MergeTracksters_recoToSim_CP_score;
  std::vector<std::vector<float>> MergeTracksters_recoToSim_CP_sharedE;
  std::vector<std::vector<uint32_t>> MergeTracksters_simToReco_CP;
  std::vector<std::vector<float>> MergeTracksters_simToReco_CP_score;
  std::vector<std::vector<float>> MergeTracksters_simToReco_CP_sharedE;

  std::vector<uint32_t> cluster_seedID;
  std::vector<float> cluster_energy;
  std::vector<float> cluster_correctedEnergy;
  std::vector<float> cluster_correctedEnergyUncertainty;
  std::vector<float> cluster_position_x;
  std::vector<float> cluster_position_y;
  std::vector<float> cluster_position_z;
  std::vector<float> cluster_position_eta;
  std::vector<float> cluster_position_phi;
  std::vector<unsigned int> cluster_layer_id;
  std::vector<int> cluster_type;
  std::vector<float> cluster_time;
  std::vector<float> cluster_timeErr;
  std::vector<float> cluster_ld;
  std::vector<float> cluster_radius;
  std::vector<uint32_t> cluster_number_of_hits;

  std::vector<int> track_ev;
  std::vector<unsigned int> track_id;
  std::vector<float> track_hgcal_x;
  std::vector<float> track_hgcal_y;
  std::vector<float> track_hgcal_z;
  std::vector<float> track_hgcal_px;
  std::vector<float> track_hgcal_py;
  std::vector<float> track_hgcal_pz;
  std::vector<float> track_hgcal_eta;
  std::vector<float> track_hgcal_phi;
  std::vector<float> track_pt;
  std::vector<int> track_charge;
  std::vector<double> track_time;
  std::vector<float> track_time_quality;
  std::vector<float> track_time_err;
  std::vector<int> track_nhits;

  // validation of mtd clusters / tracksters
  int number_of_mtdSimTracksters_;
  std::vector<int> number_of_mtdSimLCinST_;

  std::vector<int> number_of_mtdHits_;
  std::vector<std::vector<int>> mtdHits_det_;
  std::vector<int> ST_simTrack_;
  std::vector<int> ST_recoTrack_;
  std::vector<float> simTrack_pt_;
  std::vector<float> simTrack_eta_;
  std::vector<float> simTrack_phi_;

  std::vector<std::vector<bool>> simLC_is_looper_;
  std::vector<std::vector<int>> simLC_idx_;
  std::vector<std::vector<int>> simLC_simTrack_;
  std::vector<std::vector<int>> simLC_recoTrack_;
  std::vector<std::vector<float>> simLC_time_;
  std::vector<std::vector<float>> simLC_posX_;
  std::vector<std::vector<float>> simLC_posY_;
  std::vector<std::vector<float>> recocluster_time_;
  std::vector<std::vector<float>> recocluster_timeErr_;
  std::vector<std::vector<float>> recocluster_posX_;
  std::vector<std::vector<float>> recocluster_posY_;
  std::vector<std::vector<bool>> simLC_matched_;
  std::vector<std::vector<bool>> simLC_CorrectMatch_;
  std::vector<std::vector<bool>> simLC_DirectMatch_;

  TTree* trackster_tree_;
  TTree* cluster_tree_;
  TTree* graph_tree_;
  TTree* candidate_tree_;
  TTree* tracksters_merged_tree_;
  TTree* associations_tree_;
  TTree* simtrackstersSC_tree_;
  TTree* simtrackstersCP_tree_;
  TTree* tracks_tree_;
  TTree* simTICLCandidate_tree;
  TTree* MTDclusters_tree;
};

void TICLDumper::clearVariables() {
  // event info
  ev_event_ = 0;
  ntracksters_ = 0;
  nclusters_ = 0;

  trackster_time.clear();
  trackster_timeError.clear();
  trackster_regressed_energy.clear();
  trackster_raw_energy.clear();
  trackster_raw_em_energy.clear();
  trackster_raw_pt.clear();
  trackster_raw_em_pt.clear();
  trackster_barycenter_x.clear();
  trackster_barycenter_y.clear();
  trackster_barycenter_z.clear();
  trackster_EV1.clear();
  trackster_EV2.clear();
  trackster_EV3.clear();
  trackster_eVector0_x.clear();
  trackster_eVector0_y.clear();
  trackster_eVector0_z.clear();
  trackster_sigmaPCA1.clear();
  trackster_sigmaPCA2.clear();
  trackster_sigmaPCA3.clear();
  trackster_barycenter_eta.clear();
  trackster_barycenter_phi.clear();
  trackster_id_probabilities.clear();
  trackster_vertices_indexes.clear();
  trackster_vertices_x.clear();
  trackster_vertices_y.clear();
  trackster_vertices_z.clear();
  trackster_vertices_time.clear();
  trackster_vertices_timeErr.clear();
  trackster_vertices_energy.clear();
  trackster_vertices_correctedEnergy.clear();
  trackster_vertices_correctedEnergyUncertainty.clear();
  trackster_vertices_multiplicity.clear();

  stsSC_trackster_time.clear();
  stsSC_trackster_timeError.clear();
  stsSC_trackster_t0Mtd.clear();
  stsSC_trackster_t0MtdError.clear();
  stsSC_trackster_tMtd.clear();
  stsSC_trackster_tMtdError.clear();
  stsSC_trackster_speedMtd.clear();
  stsSC_trackster_tMtdSim.clear();
  stsSC_trackster_tMtdPos.clear();
  stsSC_trackster_regressed_energy.clear();
  stsSC_trackster_raw_energy.clear();
  stsSC_trackster_raw_em_energy.clear();
  stsSC_trackster_raw_pt.clear();
  stsSC_trackster_raw_em_pt.clear();
  stsSC_trackster_barycenter_x.clear();
  stsSC_trackster_barycenter_y.clear();
  stsSC_trackster_barycenter_z.clear();
  stsSC_trackster_EV1.clear();
  stsSC_trackster_EV2.clear();
  stsSC_trackster_EV3.clear();
  stsSC_trackster_eVector0_x.clear();
  stsSC_trackster_eVector0_y.clear();
  stsSC_trackster_eVector0_z.clear();
  stsSC_trackster_sigmaPCA1.clear();
  stsSC_trackster_sigmaPCA2.clear();
  stsSC_trackster_sigmaPCA3.clear();
  stsSC_trackster_barycenter_eta.clear();
  stsSC_trackster_barycenter_phi.clear();
  stsSC_pdgID.clear();
  stsSC_trackIdx.clear();
  stsSC_trackTime.clear();
  stsSC_boundaryX.clear();
  stsSC_boundaryY.clear();
  stsSC_boundaryZ.clear();
  stsSC_boundaryEta.clear();
  stsSC_boundaryPhi.clear();
  stsSC_boundaryPx.clear();
  stsSC_boundaryPy.clear();
  stsSC_boundaryPz.clear();
  stsSC_trackIdx.clear();
  stsSC_track_boundaryX.clear();
  stsSC_track_boundaryY.clear();
  stsSC_track_boundaryZ.clear();
  stsSC_track_boundaryEta.clear();
  stsSC_track_boundaryPhi.clear();
  stsSC_track_boundaryPx.clear();
  stsSC_track_boundaryPy.clear();
  stsSC_track_boundaryPz.clear();
  stsSC_trackster_id_probabilities.clear();
  stsSC_trackster_vertices_indexes.clear();
  stsSC_trackster_vertices_x.clear();
  stsSC_trackster_vertices_y.clear();
  stsSC_trackster_vertices_z.clear();
  stsSC_trackster_vertices_time.clear();
  stsSC_trackster_vertices_timeErr.clear();
  stsSC_trackster_vertices_energy.clear();
  stsSC_trackster_vertices_correctedEnergy.clear();
  stsSC_trackster_vertices_correctedEnergyUncertainty.clear();
  stsSC_trackster_vertices_multiplicity.clear();

  stsCP_trackster_time.clear();
  stsCP_trackster_timeError.clear();
  stsCP_trackster_t0Mtd.clear();
  stsCP_trackster_t0MtdError.clear();
  stsCP_trackster_tMtd.clear();
  stsCP_trackster_tMtdError.clear();
  stsCP_trackster_speedMtd.clear();
  stsCP_trackster_tMtdSim.clear();
  stsCP_trackster_tMtdPos.clear();
  stsCP_trackster_regressed_energy.clear();
  stsCP_trackster_raw_energy.clear();
  stsCP_trackster_raw_em_energy.clear();
  stsCP_trackster_raw_pt.clear();
  stsCP_trackster_raw_em_pt.clear();
  stsCP_trackster_barycenter_x.clear();
  stsCP_trackster_barycenter_y.clear();
  stsCP_trackster_barycenter_z.clear();
  stsCP_trackster_sigmaPCA1.clear();
  stsCP_trackster_sigmaPCA2.clear();
  stsCP_trackster_sigmaPCA3.clear();
  stsCP_trackster_barycenter_eta.clear();
  stsCP_trackster_barycenter_phi.clear();
  stsCP_pdgID.clear();
  stsCP_trackIdx.clear();
  stsCP_trackTime.clear();
  stsCP_boundaryX.clear();
  stsCP_boundaryY.clear();
  stsCP_boundaryZ.clear();
  stsCP_boundaryEta.clear();
  stsCP_boundaryPhi.clear();
  stsCP_boundaryPx.clear();
  stsCP_boundaryPy.clear();
  stsCP_boundaryPz.clear();
  stsCP_trackIdx.clear();
  stsCP_track_boundaryX.clear();
  stsCP_track_boundaryY.clear();
  stsCP_track_boundaryZ.clear();
  stsCP_track_boundaryEta.clear();
  stsCP_track_boundaryPhi.clear();
  stsCP_track_boundaryPx.clear();
  stsCP_track_boundaryPy.clear();
  stsCP_track_boundaryPz.clear();
  stsCP_trackster_id_probabilities.clear();
  stsCP_trackster_vertices_indexes.clear();
  stsCP_trackster_vertices_x.clear();
  stsCP_trackster_vertices_y.clear();
  stsCP_trackster_vertices_z.clear();
  stsCP_trackster_vertices_time.clear();
  stsCP_trackster_vertices_timeErr.clear();
  stsCP_trackster_vertices_energy.clear();
  stsCP_trackster_vertices_correctedEnergy.clear();
  stsCP_trackster_vertices_correctedEnergyUncertainty.clear();
  stsCP_trackster_vertices_multiplicity.clear();

  node_linked_inners.clear();
  node_linked_scores.clear();
  node_linked_outers.clear();
  isRootTrackster.clear();

  simTICLCandidate_raw_energy.clear();
  simTICLCandidate_regressed_energy.clear();
  simTICLCandidate_simTracksterCPIndex.clear();
  simTICLCandidate_boundaryX.clear();
  simTICLCandidate_boundaryY.clear();
  simTICLCandidate_boundaryZ.clear();
  simTICLCandidate_boundaryPx.clear();
  simTICLCandidate_boundaryPy.clear();
  simTICLCandidate_boundaryPz.clear();
  simTICLCandidate_trackTime.clear();
  simTICLCandidate_trackBeta.clear();
  simTICLCandidate_caloParticleMass.clear();
  simTICLCandidate_pdgId.clear();
  simTICLCandidate_charge.clear();
  simTICLCandidate_track_in_candidate.clear();

  nCandidates = 0;
  candidate_charge.clear();
  candidate_pdgId.clear();
  candidate_energy.clear();
  candidate_px.clear();
  candidate_py.clear();
  candidate_pz.clear();
  candidate_time.clear();
  candidate_time_err.clear();
  candidate_id_probabilities.clear();
  tracksters_in_candidate.clear();
  track_in_candidate.clear();

  nTrackstersMerged = 0;
  tracksters_merged_time.clear();
  tracksters_merged_timeError.clear();
  tracksters_merged_t0Mtd.clear();
  tracksters_merged_t0MtdError.clear();
  tracksters_merged_tMtd.clear();
  tracksters_merged_tMtdError.clear();
  tracksters_merged_speedMtd.clear();
  tracksters_merged_tMtdPos.clear();
  tracksters_merged_regressed_energy.clear();
  tracksters_merged_raw_energy.clear();
  tracksters_merged_raw_em_energy.clear();
  tracksters_merged_raw_pt.clear();
  tracksters_merged_raw_em_pt.clear();
  tracksters_merged_barycenter_x.clear();
  tracksters_merged_barycenter_y.clear();
  tracksters_merged_barycenter_z.clear();
  tracksters_merged_barycenter_eta.clear();
  tracksters_merged_barycenter_phi.clear();
  tracksters_merged_EV1.clear();
  tracksters_merged_EV2.clear();
  tracksters_merged_EV3.clear();
  tracksters_merged_eVector0_x.clear();
  tracksters_merged_eVector0_y.clear();
  tracksters_merged_eVector0_z.clear();
  tracksters_merged_sigmaPCA1.clear();
  tracksters_merged_sigmaPCA2.clear();
  tracksters_merged_sigmaPCA3.clear();
  tracksters_merged_id_probabilities.clear();
  tracksters_merged_time.clear();
  tracksters_merged_timeError.clear();
  tracksters_merged_regressed_energy.clear();
  tracksters_merged_raw_energy.clear();
  tracksters_merged_raw_em_energy.clear();
  tracksters_merged_raw_pt.clear();
  tracksters_merged_raw_em_pt.clear();

  tracksters_merged_vertices_indexes.clear();
  tracksters_merged_vertices_x.clear();
  tracksters_merged_vertices_y.clear();
  tracksters_merged_vertices_z.clear();
  tracksters_merged_vertices_time.clear();
  tracksters_merged_vertices_timeErr.clear();
  tracksters_merged_vertices_energy.clear();
  tracksters_merged_vertices_correctedEnergy.clear();
  tracksters_merged_vertices_correctedEnergyUncertainty.clear();
  tracksters_merged_vertices_multiplicity.clear();

  trackstersCLUE3D_recoToSim_SC.clear();
  trackstersCLUE3D_recoToSim_SC_score.clear();
  trackstersCLUE3D_recoToSim_SC_sharedE.clear();
  trackstersCLUE3D_simToReco_SC.clear();
  trackstersCLUE3D_simToReco_SC_score.clear();
  trackstersCLUE3D_simToReco_SC_sharedE.clear();

  trackstersCLUE3D_recoToSim_CP.clear();
  trackstersCLUE3D_recoToSim_CP_score.clear();
  trackstersCLUE3D_recoToSim_CP_sharedE.clear();
  trackstersCLUE3D_simToReco_CP.clear();
  trackstersCLUE3D_simToReco_CP_score.clear();
  trackstersCLUE3D_simToReco_CP_sharedE.clear();

  MergeTracksters_recoToSim_SC.clear();
  MergeTracksters_recoToSim_SC_score.clear();
  MergeTracksters_recoToSim_SC_sharedE.clear();
  MergeTracksters_simToReco_SC.clear();
  MergeTracksters_simToReco_SC_score.clear();
  MergeTracksters_simToReco_SC_sharedE.clear();

  MergeTracksters_recoToSim_CP.clear();
  MergeTracksters_recoToSim_CP_score.clear();
  MergeTracksters_recoToSim_CP_sharedE.clear();
  MergeTracksters_simToReco_CP.clear();
  MergeTracksters_simToReco_CP_score.clear();
  MergeTracksters_simToReco_CP_sharedE.clear();

  nsimTrackstersSC = 0;

  cluster_seedID.clear();
  cluster_energy.clear();
  cluster_correctedEnergy.clear();
  cluster_correctedEnergyUncertainty.clear();
  cluster_position_x.clear();
  cluster_position_y.clear();
  cluster_position_z.clear();
  cluster_position_eta.clear();
  cluster_position_phi.clear();
  cluster_layer_id.clear();
  cluster_type.clear();
  cluster_time.clear();
  cluster_timeErr.clear();
  cluster_ld.clear();
  cluster_radius.clear();
  cluster_number_of_hits.clear();

  track_ev.clear();
  track_id.clear();
  track_hgcal_x.clear();
  track_hgcal_y.clear();
  track_hgcal_z.clear();
  track_hgcal_eta.clear();
  track_hgcal_phi.clear();
  track_hgcal_px.clear();
  track_hgcal_py.clear();
  track_hgcal_pz.clear();
  track_pt.clear();
  track_charge.clear();
  track_time.clear();
  track_time_quality.clear();
  track_time_err.clear();
  track_nhits.clear();

  number_of_mtdSimTracksters_ = 0;
  number_of_mtdSimLCinST_.clear();
  
  number_of_mtdHits_.clear();
  mtdHits_det_.clear();
  ST_simTrack_.clear();
  ST_recoTrack_.clear();
  simTrack_pt_.clear();
  simTrack_eta_.clear();
  simTrack_phi_.clear();
  
  simLC_is_looper_.clear();
  simLC_idx_.clear();
  simLC_simTrack_.clear();
  simLC_recoTrack_.clear();
  simLC_time_.clear();
  simLC_posX_.clear();
  simLC_posY_.clear();
  recocluster_time_.clear();
  recocluster_timeErr_.clear();
  recocluster_posX_.clear();
  recocluster_posY_.clear();
  simLC_matched_.clear();
  simLC_CorrectMatch_.clear();
  simLC_DirectMatch_.clear();
};

TICLDumper::TICLDumper(const edm::ParameterSet& ps)
    : tracksters_token_(consumes<std::vector<ticl::Trackster>>(
          ps.getParameter<edm::InputTag>("trackstersclue3d"))),
      layer_clusters_token_(consumes<std::vector<reco::CaloCluster>>(
          ps.getParameter<edm::InputTag>("layerClusters"))),
      ticl_graph_token_(
          consumes<TICLGraph>(ps.getParameter<edm::InputTag>("ticlgraph"))),
      ticl_candidates_token_(consumes<std::vector<TICLCandidate>>(
          ps.getParameter<edm::InputTag>("ticlcandidates"))),
      tracks_token_(consumes<std::vector<reco::Track>>(
          ps.getParameter<edm::InputTag>("tracks"))),
      tracks_time_token_(consumes<edm::ValueMap<float>>(
          ps.getParameter<edm::InputTag>("tracksTime"))),
      tracks_time_quality_token_(consumes<edm::ValueMap<float>>(
          ps.getParameter<edm::InputTag>("tracksTimeQual"))),
      tracks_time_err_token_(consumes<edm::ValueMap<float>>(
          ps.getParameter<edm::InputTag>("tracksTimeErr"))),
      tracksters_merged_token_(consumes<std::vector<ticl::Trackster>>(
          ps.getParameter<edm::InputTag>("trackstersmerged"))),
      clustersTime_token_(consumes<edm::ValueMap<std::pair<float, float>>>(
          ps.getParameter<edm::InputTag>("layer_clustersTime"))),
      caloGeometry_token_(esConsumes<CaloGeometry, CaloGeometryRecord,
                                     edm::Transition::BeginRun>()),
      simTracksters_SC_token_(consumes<std::vector<ticl::Trackster>>(
          ps.getParameter<edm::InputTag>("simtrackstersSC"))),
      simTracksters_CP_token_(consumes<std::vector<ticl::Trackster>>(
          ps.getParameter<edm::InputTag>("simtrackstersCP"))),
      simTICLCandidate_token_(consumes<std::vector<TICLCandidate>>(
          ps.getParameter<edm::InputTag>("simTICLCandidates"))),
      tsRecoToSimSC_token_(consumes<hgcal::RecoToSimCollectionSimTracksters>(
          ps.getParameter<edm::InputTag>("recoToSimAssociatorSC"))),
      tsSimToRecoSC_token_(consumes<hgcal::SimToRecoCollectionSimTracksters>(
          ps.getParameter<edm::InputTag>("simToRecoAssociatorSC"))),
      tsRecoToSimCP_token_(consumes<hgcal::RecoToSimCollectionSimTracksters>(
          ps.getParameter<edm::InputTag>("recoToSimAssociatorCP"))),
      tsSimToRecoCP_token_(consumes<hgcal::SimToRecoCollectionSimTracksters>(
          ps.getParameter<edm::InputTag>("simToRecoAssociatorCP"))),
      MergeRecoToSimSC_token_(consumes<hgcal::RecoToSimCollectionSimTracksters>(
          ps.getParameter<edm::InputTag>("MergerecoToSimAssociatorSC"))),
      MergeSimToRecoSC_token_(consumes<hgcal::SimToRecoCollectionSimTracksters>(
          ps.getParameter<edm::InputTag>("MergesimToRecoAssociatorSC"))),
      MergeRecoToSimCP_token_(consumes<hgcal::RecoToSimCollectionSimTracksters>(
          ps.getParameter<edm::InputTag>("MergerecoToSimAssociatorCP"))),
      MergeSimToRecoCP_token_(consumes<hgcal::SimToRecoCollectionSimTracksters>(
          ps.getParameter<edm::InputTag>("MergesimToRecoAssociatorCP"))),
      simclusters_token_(
          consumes(ps.getParameter<edm::InputTag>("simclusters"))),
      caloparticles_token_(
          consumes(ps.getParameter<edm::InputTag>("caloparticles"))),
      MTDSimTrackstersToken_(consumes<MtdSimTracksterCollection>(ps.getParameter<edm::InputTag>("MtdSimTracksters"))),
      MTDSimLayerClustersToken_(consumes<MtdSimLayerClusterCollection>(ps.getParameter<edm::InputTag>("MtdSimLayerCluster"))),
      assocTpToTrackToken_(consumes(ps.getParameter<edm::InputTag>("tpToTrack"))),
      associationSimTrackToTPToken_(consumes(ps.getParameter<edm::InputTag>("simTrackToTPMap"))),
      mtdTracksToken_(consumes<reco::TrackCollection>(ps.getParameter<edm::InputTag>("tracksWithMtd"))),
      mtdTracksAssocToken_(consumes<edm::ValueMap<int>>(ps.getParameter<edm::InputTag>("mtdtrackAssocSrc"))),
      btlRecHitsToken_(consumes<FTLRecHitCollection>(ps.getParameter<edm::InputTag>("BTLrecHits"))),
      etlRecHitsToken_(consumes<FTLRecHitCollection>(ps.getParameter<edm::InputTag>("ETLrecHits"))),
      btlClustersToken_(consumes<FTLClusterCollection>(ps.getParameter<edm::InputTag>("BTLclusters"))),
      etlClustersToken_(consumes<FTLClusterCollection>(ps.getParameter<edm::InputTag>("ETLclusters"))),
      geometry_token_(esConsumes<CaloGeometry, CaloGeometryRecord,
                                 edm::Transition::BeginRun>()),
      detector_(ps.getParameter<std::string>("detector")),
      propName_(ps.getParameter<std::string>("propagator")),
      bfield_token_(esConsumes<MagneticField, IdealMagneticFieldRecord,
                               edm::Transition::BeginRun>()),
      propagator_token_(esConsumes<Propagator, TrackingComponentsRecord,
                                   edm::Transition::BeginRun>(
          edm::ESInputTag("", propName_))) {
  std::string detectorName_ =
      (detector_ == "HFNose") ? "HGCalHFNoseSensitive" : "HGCalEESensitive";
  hdc_token_ =
      esConsumes<HGCalDDDConstants, IdealGeometryRecord,
                 edm::Transition::BeginRun>(edm::ESInputTag("", detectorName_));
};

TICLDumper::~TICLDumper() {
  clearVariables();
};

void TICLDumper::beginRun(edm::Run const&, edm::EventSetup const& es) {
  const CaloGeometry& geom = es.getData(caloGeometry_token_);
  rhtools_.setGeometry(geom);

  edm::ESHandle<HGCalDDDConstants> hdc = es.getHandle(hdc_token_);
  hgcons_ = hdc.product();
  edm::ESHandle<MagneticField> bfield_ = es.getHandle(bfield_token_);
  edm::ESHandle<Propagator> propagator = es.getHandle(propagator_token_);
  initialize(hgcons_, rhtools_, bfield_, propagator);
}

// Define tree and branches
void TICLDumper::beginJob() {
  edm::Service<TFileService> fs;
  trackster_tree_ = fs->make<TTree>("tracksters", "TICL tracksters");
  cluster_tree_ = fs->make<TTree>("clusters", "TICL tracksters");
  graph_tree_ = fs->make<TTree>("graph", "TICL graph");
  candidate_tree_ = fs->make<TTree>("candidates", "TICL candidates");
  tracksters_merged_tree_ =
      fs->make<TTree>("trackstersMerged", "TICL tracksters merged");
  associations_tree_ = fs->make<TTree>("associations", "Associations");
  simtrackstersSC_tree_ =
      fs->make<TTree>("simtrackstersSC", "TICL simTracksters SC");
  simtrackstersCP_tree_ =
      fs->make<TTree>("simtrackstersCP", "TICL simTracksters CP");
  tracks_tree_ = fs->make<TTree>("tracks", "Tracks");
  simTICLCandidate_tree =
      fs->make<TTree>("simTICLCandidate", "Sim TICL Candidate");
  MTDclusters_tree =
      fs->make<TTree>("MTDclusters", "MTD clusters");

  simTICLCandidate_tree->Branch("simTICLCandidate_raw_energy", &simTICLCandidate_raw_energy);
  simTICLCandidate_tree->Branch("simTICLCandidate_regressed_energy", &simTICLCandidate_regressed_energy);
  simTICLCandidate_tree->Branch("simTICLCandidate_simTracksterCPIndex", &simTICLCandidate_simTracksterCPIndex);
  simTICLCandidate_tree->Branch("simTICLCandidate_boundaryX", &simTICLCandidate_boundaryX);
  simTICLCandidate_tree->Branch("simTICLCandidate_boundaryY", &simTICLCandidate_boundaryY);
  simTICLCandidate_tree->Branch("simTICLCandidate_boundaryZ", &simTICLCandidate_boundaryZ);
  simTICLCandidate_tree->Branch("simTICLCandidate_boundaryPx", &simTICLCandidate_boundaryPx);
  simTICLCandidate_tree->Branch("simTICLCandidate_boundaryPy", &simTICLCandidate_boundaryPy);
  simTICLCandidate_tree->Branch("simTICLCandidate_boundaryPz", &simTICLCandidate_boundaryPz);
  simTICLCandidate_tree->Branch("simTICLCandidate_trackTime", &simTICLCandidate_trackTime);
  simTICLCandidate_tree->Branch("simTICLCandidate_trackBeta", &simTICLCandidate_trackBeta);
  simTICLCandidate_tree->Branch("simTICLCandidate_caloParticleMass", &simTICLCandidate_caloParticleMass);
  simTICLCandidate_tree->Branch("simTICLCandidate_pdgId", &simTICLCandidate_pdgId);
  simTICLCandidate_tree->Branch("simTICLCandidate_charge", &simTICLCandidate_charge);
  simTICLCandidate_tree->Branch("simTICLCandidate_track_in_candidate", &simTICLCandidate_track_in_candidate);

  trackster_tree_->Branch("event", &ev_event_);
  trackster_tree_->Branch("NClusters", &nclusters_);
  trackster_tree_->Branch("NTracksters", &ntracksters_);
  trackster_tree_->Branch("time", &trackster_time);
  trackster_tree_->Branch("timeError", &trackster_timeError);
  trackster_tree_->Branch("regressed_energy", &trackster_regressed_energy);
  trackster_tree_->Branch("raw_energy", &trackster_raw_energy);
  trackster_tree_->Branch("raw_em_energy", &trackster_raw_em_energy);
  trackster_tree_->Branch("raw_pt", &trackster_raw_pt);
  trackster_tree_->Branch("raw_em_pt", &trackster_raw_em_pt);
  trackster_tree_->Branch("barycenter_x", &trackster_barycenter_x);
  trackster_tree_->Branch("barycenter_y", &trackster_barycenter_y);
  trackster_tree_->Branch("barycenter_z", &trackster_barycenter_z);
  trackster_tree_->Branch("trackster_barycenter_eta",
                          &trackster_barycenter_eta);
  trackster_tree_->Branch("trackster_barycenter_phi",
                          &trackster_barycenter_phi);
  trackster_tree_->Branch("EV1", &trackster_EV1);
  trackster_tree_->Branch("EV2", &trackster_EV2);
  trackster_tree_->Branch("EV3", &trackster_EV3);
  trackster_tree_->Branch("eVector0_x", &trackster_eVector0_x);
  trackster_tree_->Branch("eVector0_y", &trackster_eVector0_y);
  trackster_tree_->Branch("eVector0_z", &trackster_eVector0_z);
  trackster_tree_->Branch("sigmaPCA1", &trackster_sigmaPCA1);
  trackster_tree_->Branch("sigmaPCA2", &trackster_sigmaPCA2);
  trackster_tree_->Branch("sigmaPCA3", &trackster_sigmaPCA3);
  trackster_tree_->Branch("id_probabilities", &trackster_id_probabilities);
  trackster_tree_->Branch("vertices_indexes", &trackster_vertices_indexes);
  trackster_tree_->Branch("vertices_x", &trackster_vertices_x);
  trackster_tree_->Branch("vertices_y", &trackster_vertices_y);
  trackster_tree_->Branch("vertices_z", &trackster_vertices_z);
  trackster_tree_->Branch("vertices_time", &trackster_vertices_time);
  trackster_tree_->Branch("vertices_timeErr", &trackster_vertices_timeErr);
  trackster_tree_->Branch("vertices_energy", &trackster_vertices_energy);
  trackster_tree_->Branch("vertices_correctedEnergy",
                          &trackster_vertices_correctedEnergy);
  trackster_tree_->Branch("vertices_correctedEnergyUncertainty",
                          &trackster_vertices_correctedEnergyUncertainty);
  trackster_tree_->Branch("vertices_multiplicity",
                          &trackster_vertices_multiplicity);  // NEW

  simtrackstersSC_tree_->Branch("event", &ev_event_);
  simtrackstersSC_tree_->Branch("NTracksters", &stsSC_ntracksters_);
  simtrackstersSC_tree_->Branch("time", &stsSC_trackster_time);
  simtrackstersSC_tree_->Branch("timeError", &stsSC_trackster_timeError);
  simtrackstersSC_tree_->Branch("t0Mtd", &stsSC_trackster_t0Mtd);
  simtrackstersSC_tree_->Branch("t0MtdError", &stsSC_trackster_t0MtdError);
  simtrackstersSC_tree_->Branch("tMtd", &stsSC_trackster_tMtd);
  simtrackstersSC_tree_->Branch("tMtdError", &stsSC_trackster_tMtdError);
  simtrackstersSC_tree_->Branch("speedMtd", &stsSC_trackster_speedMtd);
  simtrackstersSC_tree_->Branch("tMtdSim", &stsSC_trackster_tMtdSim);
  simtrackstersSC_tree_->Branch("tMtdPos", &stsSC_trackster_tMtdPos);
  simtrackstersSC_tree_->Branch("regressed_energy",
                                &stsSC_trackster_regressed_energy);
  simtrackstersSC_tree_->Branch("raw_energy", &stsSC_trackster_raw_energy);
  simtrackstersSC_tree_->Branch("raw_em_energy",
                                &stsSC_trackster_raw_em_energy);
  simtrackstersSC_tree_->Branch("raw_pt", &stsSC_trackster_raw_pt);
  simtrackstersSC_tree_->Branch("raw_em_pt", &stsSC_trackster_raw_em_pt);
  simtrackstersSC_tree_->Branch("barycenter_x", &stsSC_trackster_barycenter_x);
  simtrackstersSC_tree_->Branch("barycenter_y", &stsSC_trackster_barycenter_y);
  simtrackstersSC_tree_->Branch("barycenter_z", &stsSC_trackster_barycenter_z);
  simtrackstersSC_tree_->Branch("trackster_barycenter_eta",
                                &stsSC_trackster_barycenter_eta);
  simtrackstersSC_tree_->Branch("trackster_barycenter_phi",
                                &stsSC_trackster_barycenter_phi);
  simtrackstersSC_tree_->Branch("EV1", &stsSC_trackster_EV1);
  simtrackstersSC_tree_->Branch("EV2", &stsSC_trackster_EV2);
  simtrackstersSC_tree_->Branch("EV3", &stsSC_trackster_EV3);
  simtrackstersSC_tree_->Branch("eVector0_x", &stsSC_trackster_eVector0_x);
  simtrackstersSC_tree_->Branch("eVector0_y", &stsSC_trackster_eVector0_y);
  simtrackstersSC_tree_->Branch("eVector0_z", &stsSC_trackster_eVector0_z);
  simtrackstersSC_tree_->Branch("sigmaPCA1", &stsSC_trackster_sigmaPCA1);
  simtrackstersSC_tree_->Branch("sigmaPCA2", &stsSC_trackster_sigmaPCA2);
  simtrackstersSC_tree_->Branch("sigmaPCA3", &stsSC_trackster_sigmaPCA3);
  simtrackstersSC_tree_->Branch("pdgID", &stsSC_pdgID);
  simtrackstersSC_tree_->Branch("trackIdx", &stsSC_trackIdx);
  simtrackstersSC_tree_->Branch("trackTime", &stsSC_trackTime);
  simtrackstersSC_tree_->Branch("boundaryX", &stsSC_boundaryX);
  simtrackstersSC_tree_->Branch("boundaryY", &stsSC_boundaryY);
  simtrackstersSC_tree_->Branch("boundaryZ", &stsSC_boundaryZ);
  simtrackstersSC_tree_->Branch("boundaryEta", &stsSC_boundaryEta);
  simtrackstersSC_tree_->Branch("boundaryPhi", &stsSC_boundaryPhi);
  simtrackstersSC_tree_->Branch("boundaryPx", &stsSC_boundaryPx);
  simtrackstersSC_tree_->Branch("boundaryPy", &stsSC_boundaryPy);
  simtrackstersSC_tree_->Branch("boundaryPz", &stsSC_boundaryPz);
  simtrackstersSC_tree_->Branch("trackIdx", &stsSC_trackIdx);
  simtrackstersSC_tree_->Branch("track_boundaryX", &stsSC_track_boundaryX);
  simtrackstersSC_tree_->Branch("track_boundaryY", &stsSC_track_boundaryY);
  simtrackstersSC_tree_->Branch("track_boundaryZ", &stsSC_track_boundaryZ);
  simtrackstersSC_tree_->Branch("track_boundaryEta", &stsSC_track_boundaryEta);
  simtrackstersSC_tree_->Branch("track_boundaryPhi", &stsSC_track_boundaryPhi);
  simtrackstersSC_tree_->Branch("track_boundaryPx", &stsSC_track_boundaryPx);
  simtrackstersSC_tree_->Branch("track_boundaryPy", &stsSC_track_boundaryPy);
  simtrackstersSC_tree_->Branch("track_boundaryPz", &stsSC_track_boundaryPz);
  simtrackstersSC_tree_->Branch("id_probabilities",
                                &stsSC_trackster_id_probabilities);
  simtrackstersSC_tree_->Branch("vertices_indexes",
                                &stsSC_trackster_vertices_indexes);
  simtrackstersSC_tree_->Branch("vertices_x", &stsSC_trackster_vertices_x);
  simtrackstersSC_tree_->Branch("vertices_y", &stsSC_trackster_vertices_y);
  simtrackstersSC_tree_->Branch("vertices_z", &stsSC_trackster_vertices_z);
  simtrackstersSC_tree_->Branch("vertices_time",
                                &stsSC_trackster_vertices_time);
  simtrackstersSC_tree_->Branch("vertices_timeErr",
                                &stsSC_trackster_vertices_timeErr);
  simtrackstersSC_tree_->Branch("vertices_energy",
                                &stsSC_trackster_vertices_energy);
  simtrackstersSC_tree_->Branch("vertices_correctedEnergy",
                                &stsSC_trackster_vertices_correctedEnergy);
  simtrackstersSC_tree_->Branch(
      "vertices_correctedEnergyUncertainty",
      &stsSC_trackster_vertices_correctedEnergyUncertainty);
  simtrackstersSC_tree_->Branch("vertices_multiplicity",
                                &stsSC_trackster_vertices_multiplicity);
  simtrackstersSC_tree_->Branch("NsimTrackstersSC", &nsimTrackstersSC);

  simtrackstersCP_tree_->Branch("event", &ev_event_);
  simtrackstersCP_tree_->Branch("NTracksters", &stsCP_ntracksters_);
  simtrackstersCP_tree_->Branch("time", &stsCP_trackster_time);
  simtrackstersCP_tree_->Branch("timeError", &stsCP_trackster_timeError);
  simtrackstersCP_tree_->Branch("t0Mtd", &stsCP_trackster_t0Mtd);
  simtrackstersCP_tree_->Branch("t0MtdError", &stsCP_trackster_t0MtdError);
  simtrackstersCP_tree_->Branch("tMtd", &stsCP_trackster_tMtd);
  simtrackstersCP_tree_->Branch("tMtdError", &stsCP_trackster_tMtdError);
  simtrackstersCP_tree_->Branch("speedMtd", &stsCP_trackster_speedMtd);
  simtrackstersCP_tree_->Branch("tMtdSim", &stsCP_trackster_tMtdSim);
  simtrackstersCP_tree_->Branch("tMtdPos", &stsCP_trackster_tMtdPos);
  simtrackstersCP_tree_->Branch("regressed_energy",
                                &stsCP_trackster_regressed_energy);
  simtrackstersCP_tree_->Branch("raw_energy", &stsCP_trackster_raw_energy);
  simtrackstersCP_tree_->Branch("raw_em_energy",
                                &stsCP_trackster_raw_em_energy);
  simtrackstersCP_tree_->Branch("raw_pt", &stsCP_trackster_raw_pt);
  simtrackstersCP_tree_->Branch("raw_em_pt", &stsCP_trackster_raw_em_pt);
  simtrackstersCP_tree_->Branch("barycenter_x", &stsCP_trackster_barycenter_x);
  simtrackstersCP_tree_->Branch("barycenter_y", &stsCP_trackster_barycenter_y);
  simtrackstersCP_tree_->Branch("barycenter_z", &stsCP_trackster_barycenter_z);
  simtrackstersCP_tree_->Branch("trackster_barycenter_eta",
                                &stsCP_trackster_barycenter_eta);
  simtrackstersCP_tree_->Branch("trackster_barycenter_phi",
                                &stsCP_trackster_barycenter_phi);
  simtrackstersCP_tree_->Branch("pdgID", &stsCP_pdgID);
  simtrackstersCP_tree_->Branch("trackIdx", &stsCP_trackIdx);
  simtrackstersCP_tree_->Branch("trackTime", &stsCP_trackTime);
  simtrackstersCP_tree_->Branch("boundaryX", &stsCP_boundaryX);
  simtrackstersCP_tree_->Branch("boundaryY", &stsCP_boundaryY);
  simtrackstersCP_tree_->Branch("boundaryZ", &stsCP_boundaryZ);
  simtrackstersCP_tree_->Branch("boundaryEta", &stsCP_boundaryEta);
  simtrackstersCP_tree_->Branch("boundaryPhi", &stsCP_boundaryPhi);
  simtrackstersCP_tree_->Branch("boundaryPx", &stsCP_boundaryPx);
  simtrackstersCP_tree_->Branch("boundaryPy", &stsCP_boundaryPy);
  simtrackstersCP_tree_->Branch("boundaryPz", &stsCP_boundaryPz);
  simtrackstersCP_tree_->Branch("trackIdx", &stsCP_trackIdx);
  simtrackstersCP_tree_->Branch("track_boundaryX", &stsCP_track_boundaryX);
  simtrackstersCP_tree_->Branch("track_boundaryY", &stsCP_track_boundaryY);
  simtrackstersCP_tree_->Branch("track_boundaryZ", &stsCP_track_boundaryZ);
  simtrackstersCP_tree_->Branch("track_boundaryEta", &stsCP_track_boundaryEta);
  simtrackstersCP_tree_->Branch("track_boundaryPhi", &stsCP_track_boundaryPhi);
  simtrackstersCP_tree_->Branch("track_boundaryPx", &stsCP_track_boundaryPx);
  simtrackstersCP_tree_->Branch("track_boundaryPy", &stsCP_track_boundaryPy);
  simtrackstersCP_tree_->Branch("track_boundaryPz", &stsCP_track_boundaryPz);
  simtrackstersCP_tree_->Branch("EV1", &stsCP_trackster_EV1);
  simtrackstersCP_tree_->Branch("EV2", &stsCP_trackster_EV2);
  simtrackstersCP_tree_->Branch("EV3", &stsCP_trackster_EV3);
  simtrackstersCP_tree_->Branch("eVector0_x", &stsCP_trackster_eVector0_x);
  simtrackstersCP_tree_->Branch("eVector0_y", &stsCP_trackster_eVector0_y);
  simtrackstersCP_tree_->Branch("eVector0_z", &stsCP_trackster_eVector0_z);
  simtrackstersCP_tree_->Branch("sigmaPCA1", &stsCP_trackster_sigmaPCA1);
  simtrackstersCP_tree_->Branch("sigmaPCA2", &stsCP_trackster_sigmaPCA2);
  simtrackstersCP_tree_->Branch("sigmaPCA3", &stsCP_trackster_sigmaPCA3);
  simtrackstersCP_tree_->Branch("id_probabilities",
                                &stsCP_trackster_id_probabilities);
  simtrackstersCP_tree_->Branch("vertices_indexes",
                                &stsCP_trackster_vertices_indexes);
  simtrackstersCP_tree_->Branch("vertices_x", &stsCP_trackster_vertices_x);
  simtrackstersCP_tree_->Branch("vertices_y", &stsCP_trackster_vertices_y);
  simtrackstersCP_tree_->Branch("vertices_z", &stsCP_trackster_vertices_z);
  simtrackstersCP_tree_->Branch("vertices_time",
                                &stsCP_trackster_vertices_time);
  simtrackstersCP_tree_->Branch("vertices_timeErr",
                                &stsCP_trackster_vertices_timeErr);
  simtrackstersCP_tree_->Branch("vertices_energy",
                                &stsCP_trackster_vertices_energy);
  simtrackstersCP_tree_->Branch("vertices_correctedEnergy",
                                &stsCP_trackster_vertices_correctedEnergy);
  simtrackstersCP_tree_->Branch(
      "vertices_correctedEnergyUncertainty",
      &stsCP_trackster_vertices_correctedEnergyUncertainty);
  simtrackstersCP_tree_->Branch("vertices_multiplicity",
                                &stsCP_trackster_vertices_multiplicity);  // NEW

  graph_tree_->Branch("linked_inners", &node_linked_inners);
  graph_tree_->Branch("linked_outers", &node_linked_outers);
  graph_tree_->Branch("isRootTrackster", &isRootTrackster);

  candidate_tree_->Branch("NCandidates", &nCandidates);
  candidate_tree_->Branch("candidate_charge", &candidate_charge);
  candidate_tree_->Branch("candidate_pdgId", &candidate_pdgId);
  candidate_tree_->Branch("candidate_id_probabilities",
                          &candidate_id_probabilities);
  candidate_tree_->Branch("candidate_time", &candidate_time);
  candidate_tree_->Branch("candidate_timeErr", &candidate_time_err);
  candidate_tree_->Branch("candidate_energy", &candidate_energy);
  candidate_tree_->Branch("candidate_px", &candidate_px);
  candidate_tree_->Branch("candidate_py", &candidate_py);
  candidate_tree_->Branch("candidate_pz", &candidate_pz);
  candidate_tree_->Branch("track_in_candidate", &track_in_candidate);
  candidate_tree_->Branch("tracksters_in_candidate", &tracksters_in_candidate);

  tracksters_merged_tree_->Branch("event", &ev_event_);
  tracksters_merged_tree_->Branch("NTracksters",
                                  &tracksters_merged_ntracksters_);
  tracksters_merged_tree_->Branch("time", &tracksters_merged_time);
  tracksters_merged_tree_->Branch("timeError", &tracksters_merged_timeError);
  tracksters_merged_tree_->Branch("t0Mtd", &tracksters_merged_t0Mtd);
  tracksters_merged_tree_->Branch("t0MtdError", &tracksters_merged_t0MtdError);
  tracksters_merged_tree_->Branch("tMtd", &tracksters_merged_tMtd);
  tracksters_merged_tree_->Branch("tMtdError", &tracksters_merged_tMtdError);
  tracksters_merged_tree_->Branch("speedMtd", &tracksters_merged_speedMtd);
  tracksters_merged_tree_->Branch("tMtdPos", &tracksters_merged_tMtdPos);
  tracksters_merged_tree_->Branch("regressed_energy",
                                  &tracksters_merged_regressed_energy);
  tracksters_merged_tree_->Branch("raw_energy", &tracksters_merged_raw_energy);
  tracksters_merged_tree_->Branch("raw_em_energy",
                                  &tracksters_merged_raw_em_energy);
  tracksters_merged_tree_->Branch("raw_pt", &tracksters_merged_raw_pt);
  tracksters_merged_tree_->Branch("raw_em_pt", &tracksters_merged_raw_em_pt);
  tracksters_merged_tree_->Branch("NTrackstersMerged", &nTrackstersMerged);
  tracksters_merged_tree_->Branch("barycenter_x",
                                  &tracksters_merged_barycenter_x);
  tracksters_merged_tree_->Branch("barycenter_y",
                                  &tracksters_merged_barycenter_y);
  tracksters_merged_tree_->Branch("barycenter_z",
                                  &tracksters_merged_barycenter_z);
  tracksters_merged_tree_->Branch("barycenter_eta",
                                  &tracksters_merged_barycenter_eta);
  tracksters_merged_tree_->Branch("barycenter_phi",
                                  &tracksters_merged_barycenter_phi);
  tracksters_merged_tree_->Branch("EV1", &tracksters_merged_EV1);
  tracksters_merged_tree_->Branch("EV2", &tracksters_merged_EV2);
  tracksters_merged_tree_->Branch("EV3", &tracksters_merged_EV3);
  tracksters_merged_tree_->Branch("eVector0_x", &tracksters_merged_eVector0_x);
  tracksters_merged_tree_->Branch("eVector0_y", &tracksters_merged_eVector0_y);
  tracksters_merged_tree_->Branch("eVector0_z", &tracksters_merged_eVector0_z);
  tracksters_merged_tree_->Branch("sigmaPCA1", &tracksters_merged_sigmaPCA1);
  tracksters_merged_tree_->Branch("sigmaPCA2", &tracksters_merged_sigmaPCA2);
  tracksters_merged_tree_->Branch("sigmaPCA3", &tracksters_merged_sigmaPCA3);
  tracksters_merged_tree_->Branch("id_probabilities",
                                  &tracksters_merged_id_probabilities);
  tracksters_merged_tree_->Branch("vertices_indexes",
                                  &tracksters_merged_vertices_indexes);
  tracksters_merged_tree_->Branch("vertices_x", &tracksters_merged_vertices_x);
  tracksters_merged_tree_->Branch("vertices_y", &tracksters_merged_vertices_y);
  tracksters_merged_tree_->Branch("vertices_z", &tracksters_merged_vertices_z);
  tracksters_merged_tree_->Branch("vertices_time",
                                  &tracksters_merged_vertices_time);
  tracksters_merged_tree_->Branch("vertices_timeErr",
                                  &tracksters_merged_vertices_timeErr);
  tracksters_merged_tree_->Branch("vertices_energy",
                                  &tracksters_merged_vertices_energy);
  tracksters_merged_tree_->Branch("vertices_correctedEnergy",
                                  &tracksters_merged_vertices_correctedEnergy);
  tracksters_merged_tree_->Branch(
      "vertices_correctedEnergyUncertainty",
      &tracksters_merged_vertices_correctedEnergyUncertainty);
  tracksters_merged_tree_->Branch(
      "vertices_multiplicity",
      &tracksters_merged_vertices_multiplicity);  // NEW

  associations_tree_->Branch("tsCLUE3D_recoToSim_SC",
                             &trackstersCLUE3D_recoToSim_SC);
  associations_tree_->Branch("tsCLUE3D_recoToSim_SC_score",
                             &trackstersCLUE3D_recoToSim_SC_score);
  associations_tree_->Branch("tsCLUE3D_recoToSim_SC_sharedE",
                             &trackstersCLUE3D_recoToSim_SC_sharedE);
  associations_tree_->Branch("tsCLUE3D_simToReco_SC",
                             &trackstersCLUE3D_simToReco_SC);
  associations_tree_->Branch("tsCLUE3D_simToReco_SC_score",
                             &trackstersCLUE3D_simToReco_SC_score);
  associations_tree_->Branch("tsCLUE3D_simToReco_SC_sharedE",
                             &trackstersCLUE3D_simToReco_SC_sharedE);

  associations_tree_->Branch("tsCLUE3D_recoToSim_CP",
                             &trackstersCLUE3D_recoToSim_CP);
  associations_tree_->Branch("tsCLUE3D_recoToSim_CP_score",
                             &trackstersCLUE3D_recoToSim_CP_score);
  associations_tree_->Branch("tsCLUE3D_recoToSim_CP_sharedE",
                             &trackstersCLUE3D_recoToSim_CP_sharedE);
  associations_tree_->Branch("tsCLUE3D_simToReco_CP",
                             &trackstersCLUE3D_simToReco_CP);
  associations_tree_->Branch("tsCLUE3D_simToReco_CP_score",
                             &trackstersCLUE3D_simToReco_CP_score);
  associations_tree_->Branch("tsCLUE3D_simToReco_CP_sharedE",
                             &trackstersCLUE3D_simToReco_CP_sharedE);

  associations_tree_->Branch("Mergetstracksters_recoToSim_SC",
                             &MergeTracksters_recoToSim_SC);
  associations_tree_->Branch("Mergetstracksters_recoToSim_SC_score",
                             &MergeTracksters_recoToSim_SC_score);
  associations_tree_->Branch("Mergetstracksters_recoToSim_SC_sharedE",
                             &MergeTracksters_recoToSim_SC_sharedE);
  associations_tree_->Branch("Mergetstracksters_simToReco_SC",
                             &MergeTracksters_simToReco_SC);
  associations_tree_->Branch("Mergetstracksters_simToReco_SC_score",
                             &MergeTracksters_simToReco_SC_score);
  associations_tree_->Branch("Mergetstracksters_simToReco_SC_sharedE",
                             &MergeTracksters_simToReco_SC_sharedE);

  associations_tree_->Branch("Mergetracksters_recoToSim_CP",
                             &MergeTracksters_recoToSim_CP);
  associations_tree_->Branch("Mergetracksters_recoToSim_CP_score",
                             &MergeTracksters_recoToSim_CP_score);
  associations_tree_->Branch("Mergetracksters_recoToSim_CP_sharedE",
                             &MergeTracksters_recoToSim_CP_sharedE);
  associations_tree_->Branch("Mergetracksters_simToReco_CP",
                             &MergeTracksters_simToReco_CP);
  associations_tree_->Branch("Mergetracksters_simToReco_CP_score",
                             &MergeTracksters_simToReco_CP_score);
  associations_tree_->Branch("Mergetracksters_simToReco_CP_sharedE",
                             &MergeTracksters_simToReco_CP_sharedE);

  cluster_tree_->Branch("seedID", &cluster_seedID);
  cluster_tree_->Branch("energy", &cluster_energy);
  cluster_tree_->Branch("correctedEnergy", &cluster_correctedEnergy);
  cluster_tree_->Branch("correctedEnergyUncertainty",
                        &cluster_correctedEnergyUncertainty);
  cluster_tree_->Branch("position_x", &cluster_position_x);
  cluster_tree_->Branch("position_y", &cluster_position_y);
  cluster_tree_->Branch("position_z", &cluster_position_z);
  cluster_tree_->Branch("position_eta", &cluster_position_eta);
  cluster_tree_->Branch("position_phi", &cluster_position_phi);
  cluster_tree_->Branch("cluster_layer_id", &cluster_layer_id);
  cluster_tree_->Branch("cluster_type", &cluster_type);
  cluster_tree_->Branch("cluster_time", &cluster_time);
  cluster_tree_->Branch("cluster_timeErr", &cluster_timeErr);
  cluster_tree_->Branch("cluster_local_density", &cluster_ld);
  cluster_tree_->Branch("cluster_radius", &cluster_radius);
  cluster_tree_->Branch("cluster_number_of_hits", &cluster_number_of_hits);

  tracks_tree_->Branch("track_ev", &track_ev);
  tracks_tree_->Branch("track_id", &track_id);
  tracks_tree_->Branch("track_pt", &track_pt);
  tracks_tree_->Branch("track_charge", &track_charge);
  tracks_tree_->Branch("track_time", &track_time);
  tracks_tree_->Branch("track_time_quality", &track_time_quality);
  tracks_tree_->Branch("track_time_err", &track_time_err);
  tracks_tree_->Branch("track_nhits", &track_nhits);

  MTDclusters_tree->Branch("number_of_mtdSimTracksters", &number_of_mtdSimTracksters_);
  MTDclusters_tree->Branch("number_of_mtdSimLCinST", &number_of_mtdSimLCinST_);
  MTDclusters_tree->Branch("number_of_mtdHits", &number_of_mtdHits_);
  MTDclusters_tree->Branch("mtdHits_det", &mtdHits_det_);
  MTDclusters_tree->Branch("ST_simTrack", &ST_simTrack_);
  MTDclusters_tree->Branch("ST_recoTrack", &ST_recoTrack_);
  MTDclusters_tree->Branch("ST_simTrack_pt", &simTrack_pt_);
  MTDclusters_tree->Branch("ST_simTrack_eta_", &simTrack_eta_);
  MTDclusters_tree->Branch("ST_simTrack_phi_", &simTrack_phi_);
  MTDclusters_tree->Branch("simLC_is_looper", &simLC_is_looper_);
  MTDclusters_tree->Branch("simLC_idx", &simLC_idx_);
  MTDclusters_tree->Branch("simLC_simTrack", &simLC_simTrack_);
  MTDclusters_tree->Branch("simLC_recoTrack", &simLC_recoTrack_);
  MTDclusters_tree->Branch("simLC_time", &simLC_time_);
  MTDclusters_tree->Branch("simLC_posX", &simLC_posX_);
  MTDclusters_tree->Branch("simLC_posY", &simLC_posY_);
  MTDclusters_tree->Branch("recocluster_time", &recocluster_time_);
  MTDclusters_tree->Branch("recocluster_timeErr", &recocluster_timeErr_);
  MTDclusters_tree->Branch("recocluster_posX", &recocluster_posX_);
  MTDclusters_tree->Branch("recocluster_posY", &recocluster_posY_);
  MTDclusters_tree->Branch("sc_matched", &simLC_matched_);
  MTDclusters_tree->Branch("sc_CorrMatched", &simLC_CorrectMatch_);
  MTDclusters_tree->Branch("sc_DirectMatched", &simLC_DirectMatch_);

  event_index = 0;
}

void TICLDumper::buildLayers() {
  // build disks at HGCal front & EM-Had interface for track propagation

  float zVal = hgcons_->waferZ(1, true);
  std::pair<float, float> rMinMax = hgcons_->rangeR(zVal, true);

  float zVal_interface = rhtools_.getPositionLayer(rhtools_.lastLayerEE()).z();
  std::pair<float, float> rMinMax_interface =
      hgcons_->rangeR(zVal_interface, true);

  for (int iSide = 0; iSide < 2; ++iSide) {
    float zSide = (iSide == 0) ? (-1. * zVal) : zVal;
    firstDisk_[iSide] = std::make_unique<GeomDet>(
        Disk::build(Disk::PositionType(0, 0, zSide), Disk::RotationType(),
                    SimpleDiskBounds(rMinMax.first, rMinMax.second, zSide - 0.5,
                                     zSide + 0.5)).get());

    zSide = (iSide == 0) ? (-1. * zVal_interface) : zVal_interface;
    interfaceDisk_[iSide] = std::make_unique<GeomDet>(Disk::build(
        Disk::PositionType(0, 0, zSide), Disk::RotationType(),
        SimpleDiskBounds(rMinMax_interface.first, rMinMax_interface.second,
                         zSide - 0.5, zSide + 0.5)).get());
  }
}

void TICLDumper::initialize(const HGCalDDDConstants* hgcons,
                            const hgcal::RecHitTools rhtools,
                            const edm::ESHandle<MagneticField> bfieldH,
                            const edm::ESHandle<Propagator> propH) {
  hgcons_ = hgcons;
  rhtools_ = rhtools;
  buildLayers();

  bfield_ = bfieldH;
  propagator_ = propH;
}
void TICLDumper::analyze(const edm::Event& event,
                         const edm::EventSetup& setup) {
  event_index++;
  clearVariables();
  auto bFieldProd = bfield_.product();
  const Propagator& prop = (*propagator_);
  // get all the tracksters
  edm::Handle<std::vector<ticl::Trackster>> tracksters_handle;
  event.getByToken(tracksters_token_, tracksters_handle);
  const auto& tracksters = *tracksters_handle;

  // get all the layer clusters
  edm::Handle<std::vector<reco::CaloCluster>> layer_clusters_h;
  event.getByToken(layer_clusters_token_, layer_clusters_h);
  const auto& clusters = *layer_clusters_h;

  edm::Handle<edm::ValueMap<std::pair<float, float>>> clustersTime_h;
  event.getByToken(clustersTime_token_, clustersTime_h);
  const auto& layerClustersTimes = *clustersTime_h;

  // TICL Graph
  edm::Handle<TICLGraph> ticl_graph_h;
  event.getByToken(ticl_graph_token_, ticl_graph_h);
  const auto& graph = *ticl_graph_h;

  // TICL Candidate
  edm::Handle<std::vector<TICLCandidate>> candidates_h;
  event.getByToken(ticl_candidates_token_, candidates_h);
  const auto& ticlcandidates = *candidates_h;

  // Track
  edm::Handle<std::vector<reco::Track>> tracks_h;
  event.getByToken(tracks_token_, tracks_h);
  const auto& tracks = *tracks_h;

  edm::Handle<edm::ValueMap<float>> trackTime_h;
  event.getByToken(tracks_time_token_, trackTime_h);
  const auto& trackTime = *trackTime_h;

  edm::Handle<edm::ValueMap<float>> trackTimeErr_h;
  event.getByToken(tracks_time_err_token_, trackTimeErr_h);
  const auto& trackTimeErr = *trackTimeErr_h;

  edm::Handle<edm::ValueMap<float>> trackTimeQual_h;
  event.getByToken(tracks_time_quality_token_, trackTimeQual_h);
  const auto& trackTimeQual = *trackTimeQual_h;

  // Tracksters merged
  edm::Handle<std::vector<ticl::Trackster>> tracksters_merged_h;
  event.getByToken(tracksters_merged_token_, tracksters_merged_h);
  const auto& trackstersmerged = *tracksters_merged_h;

  // simTracksters from SC
  edm::Handle<std::vector<ticl::Trackster>> simTrackstersSC_h;
  event.getByToken(simTracksters_SC_token_, simTrackstersSC_h);
  const auto& simTrackstersSC = *simTrackstersSC_h;

  // simTracksters from CP
  edm::Handle<std::vector<ticl::Trackster>> simTrackstersCP_h;
  event.getByToken(simTracksters_CP_token_, simTrackstersCP_h);
  const auto& simTrackstersCP = *simTrackstersCP_h;

  edm::Handle<std::vector<TICLCandidate>> simTICLCandidates_h;
  event.getByToken(simTICLCandidate_token_, simTICLCandidates_h);
  const auto& simTICLCandidates = *simTICLCandidates_h;

  // trackster reco to sim SC
  edm::Handle<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimSC_h;
  event.getByToken(tsRecoToSimSC_token_, tsRecoToSimSC_h);
  auto const& tsRecoSimSCMap = *tsRecoToSimSC_h;

  // sim simTrackster SC to reco trackster
  edm::Handle<hgcal::SimToRecoCollectionSimTracksters> tsSimToRecoSC_h;
  event.getByToken(tsSimToRecoSC_token_, tsSimToRecoSC_h);
  auto const& tsSimToRecoSCMap = *tsSimToRecoSC_h;

  // trackster reco to sim CP
  edm::Handle<hgcal::RecoToSimCollectionSimTracksters> tsRecoToSimCP_h;
  event.getByToken(tsRecoToSimCP_token_, tsRecoToSimCP_h);
  auto const& tsRecoSimCPMap = *tsRecoToSimCP_h;

  // sim simTrackster CP to reco trackster
  edm::Handle<hgcal::SimToRecoCollectionSimTracksters> tsSimToRecoCP_h;
  event.getByToken(tsSimToRecoCP_token_, tsSimToRecoCP_h);
  auto const& tsSimToRecoCPMap = *tsSimToRecoCP_h;

  edm::Handle<hgcal::RecoToSimCollectionSimTracksters> mergetsRecoToSimSC_h;
  event.getByToken(MergeRecoToSimSC_token_, mergetsRecoToSimSC_h);
  auto const& MergetsRecoSimSCMap = *mergetsRecoToSimSC_h;

  // sim simTrackster SC to reco trackster
  edm::Handle<hgcal::SimToRecoCollectionSimTracksters> mergetsSimToRecoSC_h;
  event.getByToken(MergeSimToRecoSC_token_, mergetsSimToRecoSC_h);
  auto const& MergetsSimToRecoSCMap = *mergetsSimToRecoSC_h;

  // trackster reco to sim CP
  edm::Handle<hgcal::RecoToSimCollectionSimTracksters> mergetsRecoToSimCP_h;
  event.getByToken(MergeRecoToSimCP_token_, mergetsRecoToSimCP_h);
  auto const& MergetsRecoSimCPMap = *mergetsRecoToSimCP_h;

  // sim simTrackster CP to reco trackster
  edm::Handle<hgcal::SimToRecoCollectionSimTracksters> mergetsSimToRecoCP_h;
  event.getByToken(MergeSimToRecoCP_token_, mergetsSimToRecoCP_h);
  auto const& MergetsSimToRecoCPMap = *mergetsSimToRecoCP_h;

  edm::Handle<std::vector<CaloParticle>> caloparticles_h;
  event.getByToken(caloparticles_token_, caloparticles_h);
  const auto& caloparticles = *caloparticles_h;

  const auto& simclusters = event.get(simclusters_token_);

  edm::Handle<MtdSimTracksterCollection> mtdSimTracksters_h;
  event.getByToken(MTDSimTrackstersToken_, mtdSimTracksters_h);
  const auto& mtdSimTracksters = *mtdSimTracksters_h;

  edm::Handle<MtdSimLayerClusterCollection> mtdSimLCs_h;
  event.getByToken(MTDSimLayerClustersToken_, mtdSimLCs_h);
  const auto& mtdSimLCs = *mtdSimLCs_h;

  const auto& TPtoRecoTrackMap = event.get(assocTpToTrackToken_);
  const auto& simTrackToTPMap = event.get(associationSimTrackToTPToken_);

  const auto& trackAssoc = event.get(mtdTracksAssocToken_);

  edm::Handle<reco::TrackCollection> mtdTracks_h;
  event.getByToken(mtdTracksToken_, mtdTracks_h);

  edm::Handle<FTLRecHitCollection> btlRecHits, etlRecHits;
  event.getByToken(btlRecHitsToken_, btlRecHits);
  event.getByToken(etlRecHitsToken_, etlRecHits);
  std::array<edm::Handle<FTLRecHitCollection>, 2> RecHitsHandle{{btlRecHits, etlRecHits}};

  edm::Handle<FTLClusterCollection> btlClusters_h, etlClusters_h;
  event.getByToken(btlClustersToken_, btlClusters_h);
  event.getByToken(etlClustersToken_, etlClusters_h);

  std::vector<FTLCluster> MtdRecoClusters;
  MtdRecoClusters.reserve(btlClusters_h->size() + etlClusters_h->size());

  for (const auto& btlClusters : *btlClusters_h) 
    std::copy(btlClusters.begin(), btlClusters.end(), std::back_inserter(MtdRecoClusters));

  for (const auto& etlClusters : *etlClusters_h) 
    std::copy(etlClusters.begin(), etlClusters.end(), std::back_inserter(MtdRecoClusters));

 // FTLClusterCollection etlClusters = *etlClusters_h;
 // FTLClusterCollection btlClusters = *btlClusters_h;

 // std::copy(btlClusters.begin(), btlClusters.end(), std::back_inserter(MtdRecoClusters));
 // std::copy(etlClusters.begin(), etlClusters.end(), std::back_inserter(MtdRecoClusters));

  ev_event_ = event_index;
  ntracksters_ = tracksters.size();
  nclusters_ = clusters.size();

  int t_id = 0;
  for (auto trackster_iterator = tracksters.begin();
       trackster_iterator != tracksters.end(); ++trackster_iterator) {
    // per-trackster analysis
    trackster_time.push_back(trackster_iterator->time());
    trackster_timeError.push_back(trackster_iterator->timeError());
    trackster_regressed_energy.push_back(
        trackster_iterator->regressed_energy());
    trackster_raw_energy.push_back(trackster_iterator->raw_energy());
    trackster_raw_em_energy.push_back(trackster_iterator->raw_em_energy());
    trackster_raw_pt.push_back(trackster_iterator->raw_pt());
    trackster_raw_em_pt.push_back(trackster_iterator->raw_em_pt());
    trackster_barycenter_x.push_back(trackster_iterator->barycenter().x());
    trackster_barycenter_y.push_back(trackster_iterator->barycenter().y());
    trackster_barycenter_z.push_back(trackster_iterator->barycenter().z());
    trackster_barycenter_eta.push_back(trackster_iterator->barycenter().eta());
    trackster_barycenter_phi.push_back(trackster_iterator->barycenter().phi());
    trackster_EV1.push_back(trackster_iterator->eigenvalues()[0]);
    trackster_EV2.push_back(trackster_iterator->eigenvalues()[1]);
    trackster_EV3.push_back(trackster_iterator->eigenvalues()[2]);
    trackster_eVector0_x.push_back((trackster_iterator->eigenvectors()[0]).x());
    trackster_eVector0_y.push_back((trackster_iterator->eigenvectors()[0]).y());
    trackster_eVector0_z.push_back((trackster_iterator->eigenvectors()[0]).z());
    trackster_sigmaPCA1.push_back(trackster_iterator->sigmasPCA()[0]);
    trackster_sigmaPCA2.push_back(trackster_iterator->sigmasPCA()[1]);
    trackster_sigmaPCA3.push_back(trackster_iterator->sigmasPCA()[2]);
    std::vector<float> id_probs;
    for (size_t i = 0; i < 8; i++)
      id_probs.push_back(trackster_iterator->id_probabilities(i));
    trackster_id_probabilities.push_back(id_probs);

    // Clusters
    std::vector<uint32_t> vertices_indexes;
    std::vector<float> vertices_x;
    std::vector<float> vertices_y;
    std::vector<float> vertices_z;
    std::vector<float> vertices_time;
    std::vector<float> vertices_timeErr;
    std::vector<float> vertices_energy;
    std::vector<float> vertices_correctedEnergy;
    std::vector<float> vertices_correctedEnergyUncertainty;
    for (auto idx : trackster_iterator->vertices()) {
      vertices_indexes.push_back(idx);
      auto associated_cluster = (*layer_clusters_h)[idx];
      vertices_x.push_back(associated_cluster.x());
      vertices_y.push_back(associated_cluster.y());
      vertices_z.push_back(associated_cluster.z());
      vertices_energy.push_back(associated_cluster.energy());
      vertices_correctedEnergy.push_back(associated_cluster.correctedEnergy());
      vertices_correctedEnergyUncertainty.push_back(
          associated_cluster.correctedEnergyUncertainty());
      vertices_time.push_back(layerClustersTimes.get(idx).first);
      vertices_timeErr.push_back(layerClustersTimes.get(idx).second);
    }
    trackster_vertices_indexes.push_back(vertices_indexes);
    trackster_vertices_x.push_back(vertices_x);
    trackster_vertices_y.push_back(vertices_y);
    trackster_vertices_z.push_back(vertices_z);
    trackster_vertices_time.push_back(vertices_time);
    trackster_vertices_timeErr.push_back(vertices_timeErr);
    trackster_vertices_energy.push_back(vertices_energy);
    trackster_vertices_correctedEnergy.push_back(vertices_correctedEnergy);
    trackster_vertices_correctedEnergyUncertainty.push_back(
        vertices_correctedEnergyUncertainty);

    // Multiplicity
    std::vector<float> vertices_multiplicity;
    for (auto multiplicity : trackster_iterator->vertex_multiplicity()) {
      vertices_multiplicity.push_back(multiplicity);
    }
    trackster_vertices_multiplicity.push_back(vertices_multiplicity);
    t_id += 1;
  }

  stsSC_ntracksters_ = simTrackstersSC.size();
  nclusters_ = clusters.size();

  for (auto trackster_iterator = simTrackstersSC.begin();
       trackster_iterator != simTrackstersSC.end(); ++trackster_iterator) {
    // per-trackster analysis
    stsSC_trackster_time.push_back(trackster_iterator->time());
    stsSC_trackster_timeError.push_back(trackster_iterator->timeError());
    stsSC_trackster_t0Mtd.push_back(trackster_iterator->t0Mtd());
    stsSC_trackster_t0MtdError.push_back(trackster_iterator->t0MtdError());
    stsSC_trackster_tMtd.push_back(trackster_iterator->tMtd());
    stsSC_trackster_tMtdError.push_back(trackster_iterator->tMtdError());
    stsSC_trackster_speedMtd.push_back(trackster_iterator->speed());
    stsSC_trackster_tMtdSim.push_back(trackster_iterator->MTDSimTime());
    stsSC_trackster_tMtdPos.push_back(trackster_iterator->tMtdPos());
    stsSC_trackster_regressed_energy.push_back(
        trackster_iterator->regressed_energy());
    stsSC_trackster_raw_energy.push_back(trackster_iterator->raw_energy());
    stsSC_trackster_raw_em_energy.push_back(
        trackster_iterator->raw_em_energy());
    stsSC_trackster_raw_pt.push_back(trackster_iterator->raw_pt());
    stsSC_trackster_raw_em_pt.push_back(trackster_iterator->raw_em_pt());
    stsSC_trackster_barycenter_x.push_back(
        trackster_iterator->barycenter().x());
    stsSC_trackster_barycenter_y.push_back(
        trackster_iterator->barycenter().y());
    stsSC_trackster_barycenter_z.push_back(
        trackster_iterator->barycenter().z());
    stsSC_trackster_barycenter_eta.push_back(
        trackster_iterator->barycenter().eta());
    stsSC_trackster_barycenter_phi.push_back(
        trackster_iterator->barycenter().phi());
    stsSC_trackster_EV1.push_back(trackster_iterator->eigenvalues()[0]);
    stsSC_trackster_EV2.push_back(trackster_iterator->eigenvalues()[1]);
    stsSC_trackster_EV3.push_back(trackster_iterator->eigenvalues()[2]);
    stsSC_trackster_eVector0_x.push_back(
        (trackster_iterator->eigenvectors()[0]).x());
    stsSC_trackster_eVector0_y.push_back(
        (trackster_iterator->eigenvectors()[0]).y());
    stsSC_trackster_eVector0_z.push_back(
        (trackster_iterator->eigenvectors()[0]).z());
    stsSC_trackster_sigmaPCA1.push_back(trackster_iterator->sigmasPCA()[0]);
    stsSC_trackster_sigmaPCA2.push_back(trackster_iterator->sigmasPCA()[1]);
    stsSC_trackster_sigmaPCA3.push_back(trackster_iterator->sigmasPCA()[2]);
    stsSC_pdgID.push_back(simclusters[trackster_iterator->seedIndex()].pdgId());
    auto simTrack =
        trackster_iterator->seedID() == caloparticles_h.id()
            ? caloparticles[trackster_iterator->seedIndex()].g4Tracks()[0]
            : simclusters[trackster_iterator->seedIndex()].g4Tracks()[0];
    if (simTrack.crossedBoundary()) {
      stsSC_boundaryX.push_back(simTrack.getPositionAtBoundary().x());
      stsSC_boundaryY.push_back(simTrack.getPositionAtBoundary().y());
      stsSC_boundaryZ.push_back(simTrack.getPositionAtBoundary().z());
      stsSC_boundaryEta.push_back(simTrack.getPositionAtBoundary().eta());
      stsSC_boundaryPhi.push_back(simTrack.getPositionAtBoundary().phi());
      stsSC_boundaryPx.push_back(simTrack.getMomentumAtBoundary().x());
      stsSC_boundaryPy.push_back(simTrack.getMomentumAtBoundary().y());
      stsSC_boundaryPz.push_back(simTrack.getMomentumAtBoundary().z());
    } else {
      stsSC_boundaryX.push_back(-999);
      stsSC_boundaryY.push_back(-999);
      stsSC_boundaryZ.push_back(-999);
      stsSC_boundaryEta.push_back(-999);
      stsSC_boundaryPhi.push_back(-999);
      stsSC_boundaryPx.push_back(-999);
      stsSC_boundaryPy.push_back(-999);
      stsSC_boundaryPz.push_back(-999);
    }
    auto const trackIdx = trackster_iterator->trackIdx();
				stsSC_trackIdx.push_back(trackIdx);
    if (trackIdx != -1) {
      auto track = tracks[trackIdx];
      stsSC_trackIdx.push_back(trackIdx);

      int iSide = int(track.eta() > 0);

      const auto& fts =
          trajectoryStateTransform::outerFreeState((track), bFieldProd);
      // to the HGCal front
      const auto& tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
      if (tsos.isValid()) {
        const auto& globalPos = tsos.globalPosition();
        const auto& globalMom = tsos.globalMomentum();
        stsSC_track_boundaryX.push_back(globalPos.x());
        stsSC_track_boundaryY.push_back(globalPos.y());
        stsSC_track_boundaryZ.push_back(globalPos.z());
        stsSC_track_boundaryEta.push_back(globalPos.eta());
        stsSC_track_boundaryPhi.push_back(globalPos.phi());
        stsSC_track_boundaryPx.push_back(globalMom.x());
        stsSC_track_boundaryPy.push_back(globalMom.y());
        stsSC_track_boundaryPz.push_back(globalMom.z());
        stsSC_trackTime.push_back(track.t0());
      } else {
        stsSC_track_boundaryX.push_back(-999);
        stsSC_track_boundaryY.push_back(-999);
        stsSC_track_boundaryZ.push_back(-999);
        stsSC_track_boundaryEta.push_back(-999);
        stsSC_track_boundaryPhi.push_back(-999);
        stsSC_track_boundaryPx.push_back(-999);
        stsSC_track_boundaryPy.push_back(-999);
        stsSC_track_boundaryPz.push_back(-999);
      }
    } else {
      stsSC_track_boundaryX.push_back(-999);
      stsSC_track_boundaryY.push_back(-999);
      stsSC_track_boundaryZ.push_back(-999);
      stsSC_track_boundaryEta.push_back(-999);
      stsSC_track_boundaryPhi.push_back(-999);
      stsSC_track_boundaryPx.push_back(-999);
      stsSC_track_boundaryPy.push_back(-999);
      stsSC_track_boundaryPz.push_back(-999);
    }

    std::vector<float> id_probs;
    for (size_t i = 0; i < 8; i++)
      id_probs.push_back(trackster_iterator->id_probabilities(i));
    stsSC_trackster_id_probabilities.push_back(id_probs);

    // Clusters
    std::vector<uint32_t> vertices_indexes;
    std::vector<float> vertices_x;
    std::vector<float> vertices_y;
    std::vector<float> vertices_z;
    std::vector<float> vertices_time;
    std::vector<float> vertices_timeErr;
    std::vector<float> vertices_energy;
    std::vector<float> vertices_correctedEnergy;
    std::vector<float> vertices_correctedEnergyUncertainty;
    for (auto idx : trackster_iterator->vertices()) {
      vertices_indexes.push_back(idx);
      auto associated_cluster = (*layer_clusters_h)[idx];
      vertices_x.push_back(associated_cluster.x());
      vertices_y.push_back(associated_cluster.y());
      vertices_z.push_back(associated_cluster.z());
      vertices_energy.push_back(associated_cluster.energy());
      vertices_correctedEnergy.push_back(associated_cluster.correctedEnergy());
      vertices_correctedEnergyUncertainty.push_back(
          associated_cluster.correctedEnergyUncertainty());
      vertices_time.push_back(layerClustersTimes.get(idx).first);
      vertices_timeErr.push_back(layerClustersTimes.get(idx).second);
    }
    stsSC_trackster_vertices_indexes.push_back(vertices_indexes);
    stsSC_trackster_vertices_x.push_back(vertices_x);
    stsSC_trackster_vertices_y.push_back(vertices_y);
    stsSC_trackster_vertices_z.push_back(vertices_z);
    stsSC_trackster_vertices_time.push_back(vertices_time);
    stsSC_trackster_vertices_timeErr.push_back(vertices_timeErr);
    stsSC_trackster_vertices_energy.push_back(vertices_energy);
    stsSC_trackster_vertices_correctedEnergy.push_back(
        vertices_correctedEnergy);
    stsSC_trackster_vertices_correctedEnergyUncertainty.push_back(
        vertices_correctedEnergyUncertainty);

    // Multiplicity
    std::vector<float> vertices_multiplicity;
    for (auto multiplicity : trackster_iterator->vertex_multiplicity()) {
      vertices_multiplicity.push_back(multiplicity);
    }
    stsSC_trackster_vertices_multiplicity.push_back(vertices_multiplicity);
  }

  stsCP_ntracksters_ = simTrackstersCP.size();
  nclusters_ = clusters.size();
  for (auto trackster_iterator = simTrackstersCP.begin();
       trackster_iterator != simTrackstersCP.end(); ++trackster_iterator) {
    // per-trackster analysis
    stsCP_trackster_time.push_back(trackster_iterator->time());
    stsCP_trackster_timeError.push_back(trackster_iterator->timeError());
    stsCP_trackster_t0Mtd.push_back(trackster_iterator->t0Mtd());
    stsCP_trackster_t0MtdError.push_back(trackster_iterator->t0MtdError());
    stsCP_trackster_tMtd.push_back(trackster_iterator->tMtd());
    stsCP_trackster_tMtdError.push_back(trackster_iterator->tMtdError());
    stsCP_trackster_speedMtd.push_back(trackster_iterator->speed());
    stsCP_trackster_tMtdSim.push_back(trackster_iterator->MTDSimTime());
    stsCP_trackster_tMtdPos.push_back(trackster_iterator->tMtdPos());
    stsCP_trackster_regressed_energy.push_back(
        trackster_iterator->regressed_energy());
    stsCP_trackster_raw_energy.push_back(trackster_iterator->raw_energy());
    stsCP_trackster_raw_em_energy.push_back(
        trackster_iterator->raw_em_energy());
    stsCP_trackster_raw_pt.push_back(trackster_iterator->raw_pt());
    stsCP_trackster_raw_em_pt.push_back(trackster_iterator->raw_em_pt());
    stsCP_trackster_barycenter_x.push_back(
        trackster_iterator->barycenter().x());
    stsCP_trackster_barycenter_y.push_back(
        trackster_iterator->barycenter().y());
    stsCP_trackster_barycenter_z.push_back(
        trackster_iterator->barycenter().z());
    stsCP_trackster_barycenter_eta.push_back(
        trackster_iterator->barycenter().eta());
    stsCP_trackster_barycenter_phi.push_back(
        trackster_iterator->barycenter().phi());
    stsCP_trackster_EV1.push_back(trackster_iterator->eigenvalues()[0]);
    stsCP_trackster_EV2.push_back(trackster_iterator->eigenvalues()[1]);
    stsCP_trackster_EV3.push_back(trackster_iterator->eigenvalues()[2]);
    stsCP_trackster_eVector0_x.push_back(
        (trackster_iterator->eigenvectors()[0]).x());
    stsCP_trackster_eVector0_y.push_back(
        (trackster_iterator->eigenvectors()[0]).y());
    stsCP_trackster_eVector0_z.push_back(
        (trackster_iterator->eigenvectors()[0]).z());
    stsCP_trackster_sigmaPCA1.push_back(trackster_iterator->sigmasPCA()[0]);
    stsCP_trackster_sigmaPCA2.push_back(trackster_iterator->sigmasPCA()[1]);
    stsCP_trackster_sigmaPCA3.push_back(trackster_iterator->sigmasPCA()[2]);
    stsCP_pdgID.push_back(
        caloparticles[trackster_iterator->seedIndex()].pdgId());
    auto simTrack =
        trackster_iterator->seedID() == caloparticles_h.id()
            ? caloparticles[trackster_iterator->seedIndex()].g4Tracks()[0]
            : simclusters[trackster_iterator->seedIndex()].g4Tracks()[0];
    if (simTrack.crossedBoundary()) {
      stsCP_boundaryX.push_back(simTrack.getPositionAtBoundary().x());
      stsCP_boundaryY.push_back(simTrack.getPositionAtBoundary().y());
      stsCP_boundaryZ.push_back(simTrack.getPositionAtBoundary().z());
      stsCP_boundaryEta.push_back(simTrack.getPositionAtBoundary().eta());
      stsCP_boundaryPhi.push_back(simTrack.getPositionAtBoundary().phi());
      stsCP_boundaryPx.push_back(simTrack.getMomentumAtBoundary().x());
      stsCP_boundaryPy.push_back(simTrack.getMomentumAtBoundary().y());
      stsCP_boundaryPz.push_back(simTrack.getMomentumAtBoundary().z());
    } else {
      stsCP_boundaryX.push_back(-999);
      stsCP_boundaryY.push_back(-999);
      stsCP_boundaryZ.push_back(-999);
      stsCP_boundaryEta.push_back(-999);
      stsCP_boundaryPhi.push_back(-999);
      stsCP_boundaryPx.push_back(-999);
      stsCP_boundaryPy.push_back(-999);
      stsCP_boundaryPz.push_back(-999);
    }
    auto const trackIdx = trackster_iterator->trackIdx();
				stsCP_trackIdx.push_back(trackIdx);
    if (trackIdx != -1) {
      auto track = tracks[trackIdx];
      stsCP_trackIdx.push_back(trackIdx);

      int iSide = int(track.eta() > 0);

      const auto& fts =
          trajectoryStateTransform::outerFreeState((track), bFieldProd);
      // to the HGCal front
      const auto& tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
      if (tsos.isValid()) {
        const auto& globalPos = tsos.globalPosition();
        const auto& globalMom = tsos.globalMomentum();
        stsCP_track_boundaryX.push_back(globalPos.x());
        stsCP_track_boundaryY.push_back(globalPos.y());
        stsCP_track_boundaryZ.push_back(globalPos.z());
        stsCP_track_boundaryEta.push_back(globalPos.eta());
        stsCP_track_boundaryPhi.push_back(globalPos.phi());
        stsCP_track_boundaryPx.push_back(globalMom.x());
        stsCP_track_boundaryPy.push_back(globalMom.y());
        stsCP_track_boundaryPz.push_back(globalMom.z());
        stsCP_trackTime.push_back(track.t0());
      } else {
        stsCP_track_boundaryX.push_back(-999);
        stsCP_track_boundaryY.push_back(-999);
        stsCP_track_boundaryZ.push_back(-999);
        stsCP_track_boundaryEta.push_back(-999);
        stsCP_track_boundaryPhi.push_back(-999);
        stsCP_track_boundaryPx.push_back(-999);
        stsCP_track_boundaryPy.push_back(-999);
        stsCP_track_boundaryPz.push_back(-999);
      }
    } else {
      stsCP_track_boundaryX.push_back(-999);
      stsCP_track_boundaryY.push_back(-999);
      stsCP_track_boundaryZ.push_back(-999);
      stsCP_track_boundaryEta.push_back(-999);
      stsCP_track_boundaryPhi.push_back(-999);
      stsCP_track_boundaryPx.push_back(-999);
      stsCP_track_boundaryPy.push_back(-999);
      stsCP_track_boundaryPz.push_back(-999);
    }
    std::vector<float> id_probs;
    for (size_t i = 0; i < 8; i++)
      id_probs.push_back(trackster_iterator->id_probabilities(i));
    stsCP_trackster_id_probabilities.push_back(id_probs);

    // Clusters
    std::vector<uint32_t> vertices_indexes;
    std::vector<float> vertices_x;
    std::vector<float> vertices_y;
    std::vector<float> vertices_z;
    std::vector<float> vertices_time;
    std::vector<float> vertices_timeErr;
    std::vector<float> vertices_energy;
    std::vector<float> vertices_correctedEnergy;
    std::vector<float> vertices_correctedEnergyUncertainty;
    for (auto idx : trackster_iterator->vertices()) {
      vertices_indexes.push_back(idx);
      auto associated_cluster = (*layer_clusters_h)[idx];
      vertices_x.push_back(associated_cluster.x());
      vertices_y.push_back(associated_cluster.y());
      vertices_z.push_back(associated_cluster.z());
      vertices_energy.push_back(associated_cluster.energy());
      vertices_correctedEnergy.push_back(associated_cluster.correctedEnergy());
      vertices_correctedEnergyUncertainty.push_back(
          associated_cluster.correctedEnergyUncertainty());
      vertices_time.push_back(layerClustersTimes.get(idx).first);
      vertices_timeErr.push_back(layerClustersTimes.get(idx).second);
    }
    stsCP_trackster_vertices_indexes.push_back(vertices_indexes);
    stsCP_trackster_vertices_x.push_back(vertices_x);
    stsCP_trackster_vertices_y.push_back(vertices_y);
    stsCP_trackster_vertices_z.push_back(vertices_z);
    stsCP_trackster_vertices_time.push_back(vertices_time);
    stsCP_trackster_vertices_timeErr.push_back(vertices_timeErr);
    stsCP_trackster_vertices_energy.push_back(vertices_energy);
    stsCP_trackster_vertices_correctedEnergy.push_back(
        vertices_correctedEnergy);
    stsCP_trackster_vertices_correctedEnergyUncertainty.push_back(
        vertices_correctedEnergyUncertainty);

    // Multiplicity
    std::vector<float> vertices_multiplicity;
    for (auto multiplicity : trackster_iterator->vertex_multiplicity()) {
      vertices_multiplicity.push_back(multiplicity);
    }
    stsCP_trackster_vertices_multiplicity.push_back(vertices_multiplicity);
  }

 	simTICLCandidate_track_in_candidate.resize(simTICLCandidates.size(), -1);
  for (size_t i = 0; i < simTICLCandidates.size(); ++i) {
    auto const& cand = simTICLCandidates[i];
    // auto const& cp = caloparticles[simTrackstersCP[i].seedIndex()];

    simTICLCandidate_raw_energy.push_back(cand.rawEnergy());
    simTICLCandidate_regressed_energy.push_back(cand.p4().energy());
    simTICLCandidate_pdgId.push_back(cand.pdgId());
    simTICLCandidate_charge.push_back(cand.charge());
    std::vector<int> tmpIdxVec;
    for (auto const& simTS : cand.tracksters()) {
      auto trackster_idx =
          simTS.get() - (edm::Ptr<ticl::Trackster>(simTrackstersSC_h, 0)).get();
      tmpIdxVec.push_back(trackster_idx);
    }
    simTICLCandidate_simTracksterCPIndex.push_back(tmpIdxVec);
    tmpIdxVec.clear();
    auto const& trackPtr = cand.trackPtr();
    if (!trackPtr.isNull()) {
      auto const& track = *trackPtr;
      int iSide = int(track.eta() > 0);
    		int tk_idx = trackPtr.get() - (edm::Ptr<reco::Track>(tracks_h, 0)).get();
 		   simTICLCandidate_track_in_candidate[i] = tk_idx;
  

      const auto& fts =
          trajectoryStateTransform::outerFreeState((track), bFieldProd);
      // to the HGCal front
      const auto& tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
      if (tsos.isValid()) {
        const auto& globalPos = tsos.globalPosition();
        const auto& globalMom = tsos.globalMomentum();
        simTICLCandidate_boundaryX.push_back(globalPos.x());
        simTICLCandidate_boundaryY.push_back(globalPos.y());
        simTICLCandidate_boundaryZ.push_back(globalPos.z());
        simTICLCandidate_boundaryPx.push_back(globalMom.x());
        simTICLCandidate_boundaryPy.push_back(globalMom.y());
        simTICLCandidate_boundaryPz.push_back(globalMom.z());
        simTICLCandidate_trackTime.push_back(track.t0());
        simTICLCandidate_trackBeta.push_back(track.beta());
      } else {
        simTICLCandidate_boundaryX.push_back(-999);
        simTICLCandidate_boundaryY.push_back(-999);
        simTICLCandidate_boundaryZ.push_back(-999);
        simTICLCandidate_boundaryPx.push_back(-999);
        simTICLCandidate_boundaryPy.push_back(-999);
        simTICLCandidate_boundaryPz.push_back(-999);
        simTICLCandidate_trackTime.push_back(-999);
        simTICLCandidate_trackBeta.push_back(-999);
      }
    } else {
      simTICLCandidate_boundaryX.push_back(-999);
      simTICLCandidate_boundaryY.push_back(-999);
      simTICLCandidate_boundaryZ.push_back(-999);
      simTICLCandidate_boundaryPx.push_back(-999);
      simTICLCandidate_boundaryPy.push_back(-999);
      simTICLCandidate_boundaryPz.push_back(-999);
      simTICLCandidate_trackTime.push_back(-999);
      simTICLCandidate_trackBeta.push_back(-999);
    }
  }

  node_linked_inners.resize(tracksters.size());
  node_linked_scores.resize(tracksters.size());
  node_linked_outers.resize(tracksters.size());

  int c_id = 0;

  for (auto cluster_iterator = clusters.begin();
       cluster_iterator != clusters.end(); ++cluster_iterator) {
    auto lc_seed = cluster_iterator->seed();
    cluster_seedID.push_back(lc_seed);
    cluster_energy.push_back(cluster_iterator->energy());
    cluster_correctedEnergy.push_back(cluster_iterator->correctedEnergy());
    cluster_correctedEnergyUncertainty.push_back(
        cluster_iterator->correctedEnergyUncertainty());
    cluster_position_x.push_back(cluster_iterator->x());
    cluster_position_y.push_back(cluster_iterator->y());
    cluster_position_z.push_back(cluster_iterator->z());
    cluster_position_eta.push_back(cluster_iterator->eta());
    cluster_position_phi.push_back(cluster_iterator->phi());
    auto haf = cluster_iterator->hitsAndFractions();
    auto layerId = rhtools_.getLayerWithOffset(haf[0].first);
    cluster_layer_id.push_back(layerId);
    uint32_t number_of_hits = cluster_iterator->hitsAndFractions().size();
    cluster_number_of_hits.push_back(number_of_hits);
    cluster_type.push_back(ticl::returnIndex(lc_seed, rhtools_));
    cluster_timeErr.push_back(layerClustersTimes.get(c_id).second);
    cluster_time.push_back(layerClustersTimes.get(c_id).first);
    c_id += 1;
  }

  tracksters_in_candidate.resize(ticlcandidates.size());
  track_in_candidate.resize(ticlcandidates.size(), -1);
  nCandidates = ticlcandidates.size();
  for (int i = 0; i < static_cast<int>(ticlcandidates.size()); ++i) {
    const auto& candidate = ticlcandidates[i];
    candidate_charge.push_back(candidate.charge());
    candidate_pdgId.push_back(candidate.pdgId());
    candidate_energy.push_back(candidate.energy());
    candidate_px.push_back(candidate.px());
    candidate_py.push_back(candidate.py());
    candidate_pz.push_back(candidate.pz());
    candidate_time.push_back(candidate.time());
    candidate_time_err.push_back(candidate.timeError());
    std::vector<float> id_probs;
    for (int j = 0; j < 8; j++) {
      ticl::Trackster::ParticleType type = static_cast<ticl::Trackster::ParticleType>(j);
      id_probs.push_back(candidate.id_probability(type));
    }
    candidate_id_probabilities.push_back(id_probs);

    auto trackster_ptrs = candidate.tracksters();
    auto track_ptr = candidate.trackPtr();
    for (auto ts_ptr : trackster_ptrs) {
      auto ts_idx = ts_ptr.get() -
                    (edm::Ptr<ticl::Trackster>(tracksters_handle, 0)).get();
      tracksters_in_candidate[i].push_back(ts_idx);
    }
    if (track_ptr.isNull())
      continue;
    int tk_idx = track_ptr.get() - (edm::Ptr<reco::Track>(tracks_h, 0)).get();
    track_in_candidate[i] = tk_idx;
  }

  nTrackstersMerged = trackstersmerged.size();
  for (auto trackster_iterator = trackstersmerged.begin();
       trackster_iterator != trackstersmerged.end(); ++trackster_iterator) {
    tracksters_merged_time.push_back(trackster_iterator->time());
    tracksters_merged_timeError.push_back(trackster_iterator->timeError());
    tracksters_merged_t0Mtd.push_back(trackster_iterator->t0Mtd());
    tracksters_merged_t0MtdError.push_back(trackster_iterator->t0MtdError());
    tracksters_merged_tMtd.push_back(trackster_iterator->tMtd());
    tracksters_merged_tMtdError.push_back(trackster_iterator->tMtdError());
    tracksters_merged_speedMtd.push_back(trackster_iterator->speed());
    tracksters_merged_tMtdPos.push_back(trackster_iterator->tMtdPos());
    tracksters_merged_regressed_energy.push_back(
        trackster_iterator->regressed_energy());
    tracksters_merged_raw_energy.push_back(trackster_iterator->raw_energy());
    tracksters_merged_raw_em_energy.push_back(
        trackster_iterator->raw_em_energy());
    tracksters_merged_raw_pt.push_back(trackster_iterator->raw_pt());
    tracksters_merged_raw_em_pt.push_back(trackster_iterator->raw_em_pt());
    tracksters_merged_barycenter_x.push_back(
        trackster_iterator->barycenter().x());
    tracksters_merged_barycenter_y.push_back(
        trackster_iterator->barycenter().y());
    tracksters_merged_barycenter_z.push_back(
        trackster_iterator->barycenter().z());
    tracksters_merged_barycenter_eta.push_back(
        trackster_iterator->barycenter().eta());
    tracksters_merged_barycenter_phi.push_back(
        trackster_iterator->barycenter().phi());
    tracksters_merged_EV1.push_back(trackster_iterator->eigenvalues()[0]);
    tracksters_merged_EV2.push_back(trackster_iterator->eigenvalues()[1]);
    tracksters_merged_EV3.push_back(trackster_iterator->eigenvalues()[2]);
    tracksters_merged_eVector0_x.push_back(
        (trackster_iterator->eigenvectors()[0]).x());
    tracksters_merged_eVector0_y.push_back(
        (trackster_iterator->eigenvectors()[0]).y());
    tracksters_merged_eVector0_z.push_back(
        (trackster_iterator->eigenvectors()[0]).z());
    tracksters_merged_sigmaPCA1.push_back(trackster_iterator->sigmasPCA()[0]);
    tracksters_merged_sigmaPCA2.push_back(trackster_iterator->sigmasPCA()[1]);
    tracksters_merged_sigmaPCA3.push_back(trackster_iterator->sigmasPCA()[2]);

    std::vector<float> id_probs;
    for (size_t i = 0; i < 8; i++)
      id_probs.push_back(trackster_iterator->id_probabilities(i));
    tracksters_merged_id_probabilities.push_back(id_probs);

    std::vector<uint32_t> vertices_indexes;
    std::vector<float> vertices_x;
    std::vector<float> vertices_y;
    std::vector<float> vertices_z;
    std::vector<float> vertices_time;
    std::vector<float> vertices_timeErr;
    std::vector<float> vertices_energy;
    std::vector<float> vertices_correctedEnergy;
    std::vector<float> vertices_correctedEnergyUncertainty;
    for (auto idx : trackster_iterator->vertices()) {
      vertices_indexes.push_back(idx);
      auto associated_cluster = (*layer_clusters_h)[idx];
      vertices_x.push_back(associated_cluster.x());
      vertices_y.push_back(associated_cluster.y());
      vertices_z.push_back(associated_cluster.z());
      vertices_energy.push_back(associated_cluster.energy());
      vertices_correctedEnergy.push_back(associated_cluster.correctedEnergy());
      vertices_correctedEnergyUncertainty.push_back(
          associated_cluster.correctedEnergyUncertainty());
      vertices_time.push_back(layerClustersTimes.get(idx).first);
      vertices_timeErr.push_back(layerClustersTimes.get(idx).second);
    }
    tracksters_merged_vertices_indexes.push_back(vertices_indexes);
    tracksters_merged_vertices_x.push_back(vertices_x);
    tracksters_merged_vertices_y.push_back(vertices_y);
    tracksters_merged_vertices_z.push_back(vertices_z);
    tracksters_merged_vertices_time.push_back(vertices_time);
    tracksters_merged_vertices_timeErr.push_back(vertices_timeErr);
    tracksters_merged_vertices_energy.push_back(vertices_energy);
    tracksters_merged_vertices_correctedEnergy.push_back(
        vertices_correctedEnergy);
    tracksters_merged_vertices_correctedEnergyUncertainty.push_back(
        vertices_correctedEnergyUncertainty);
  }

  // Tackster reco->sim associations
  trackstersCLUE3D_recoToSim_SC.resize(tracksters.size());
  trackstersCLUE3D_recoToSim_SC_score.resize(tracksters.size());
  trackstersCLUE3D_recoToSim_SC_sharedE.resize(tracksters.size());
  for (size_t i = 0; i < tracksters.size(); ++i) {
    const edm::Ref<ticl::TracksterCollection> tsRef(tracksters_handle, i);

    // CLUE3D -> STS-SC
    const auto stsSC_iter = tsRecoSimSCMap.find(tsRef);
    if (stsSC_iter != tsRecoSimSCMap.end()) {
      const auto& stsSCassociated = stsSC_iter->val;
      for (auto& sts : stsSCassociated) {
        auto sts_id = (sts.first).get() - (edm::Ref<ticl::TracksterCollection>(
                                               simTrackstersSC_h, 0)).get();
        trackstersCLUE3D_recoToSim_SC[i].push_back(sts_id);
        trackstersCLUE3D_recoToSim_SC_score[i].push_back(sts.second.second);
        trackstersCLUE3D_recoToSim_SC_sharedE[i].push_back(sts.second.first);
      }
    }
  }

  // SimTracksters
  nsimTrackstersSC = simTrackstersSC.size();
  trackstersCLUE3D_simToReco_SC.resize(nsimTrackstersSC);
  trackstersCLUE3D_simToReco_SC_score.resize(nsimTrackstersSC);
  trackstersCLUE3D_simToReco_SC_sharedE.resize(nsimTrackstersSC);
  for (size_t i = 0; i < nsimTrackstersSC; ++i) {
    const edm::Ref<ticl::TracksterCollection> stsSCRef(simTrackstersSC_h, i);

    // STS-SC -> CLUE3D
    const auto ts_iter = tsSimToRecoSCMap.find(stsSCRef);
    if (ts_iter != tsSimToRecoSCMap.end()) {
      const auto& tsAssociated = ts_iter->val;
      for (auto& ts : tsAssociated) {
        auto ts_idx = (ts.first).get() - (edm::Ref<ticl::TracksterCollection>(
                                              tracksters_handle, 0)).get();
        trackstersCLUE3D_simToReco_SC[i].push_back(ts_idx);
        trackstersCLUE3D_simToReco_SC_score[i].push_back(ts.second.second);
        trackstersCLUE3D_simToReco_SC_sharedE[i].push_back(ts.second.first);
      }
    }
  }

  // Tackster reco->sim associations
  trackstersCLUE3D_recoToSim_CP.resize(tracksters.size());
  trackstersCLUE3D_recoToSim_CP_score.resize(tracksters.size());
  trackstersCLUE3D_recoToSim_CP_sharedE.resize(tracksters.size());
  for (size_t i = 0; i < tracksters.size(); ++i) {
    const edm::Ref<ticl::TracksterCollection> tsRef(tracksters_handle, i);

    // CLUE3D -> STS-CP
    const auto stsCP_iter = tsRecoSimCPMap.find(tsRef);
    if (stsCP_iter != tsRecoSimCPMap.end()) {
      const auto& stsCPassociated = stsCP_iter->val;
      for (auto& sts : stsCPassociated) {
        auto sts_id = (sts.first).get() - (edm::Ref<ticl::TracksterCollection>(
                                               simTrackstersCP_h, 0)).get();
        trackstersCLUE3D_recoToSim_CP[i].push_back(sts_id);
        trackstersCLUE3D_recoToSim_CP_score[i].push_back(sts.second.second);
        trackstersCLUE3D_recoToSim_CP_sharedE[i].push_back(sts.second.first);
      }
    }
  }

  // SimTracksters
  nsimTrackstersCP = simTrackstersCP.size();
  trackstersCLUE3D_simToReco_CP.resize(nsimTrackstersCP);
  trackstersCLUE3D_simToReco_CP_score.resize(nsimTrackstersCP);
  trackstersCLUE3D_simToReco_CP_sharedE.resize(nsimTrackstersCP);
  for (size_t i = 0; i < nsimTrackstersCP; ++i) {
    const edm::Ref<ticl::TracksterCollection> stsCPRef(simTrackstersCP_h, i);

    // STS-CP -> CLUE3D
    const auto ts_iter = tsSimToRecoCPMap.find(stsCPRef);
    if (ts_iter != tsSimToRecoCPMap.end()) {
      const auto& tsAssociated = ts_iter->val;
      for (auto& ts : tsAssociated) {
        auto ts_idx = (ts.first).get() - (edm::Ref<ticl::TracksterCollection>(
                                              tracksters_handle, 0)).get();
        trackstersCLUE3D_simToReco_CP[i].push_back(ts_idx);
        trackstersCLUE3D_simToReco_CP_score[i].push_back(ts.second.second);
        trackstersCLUE3D_simToReco_CP_sharedE[i].push_back(ts.second.first);
      }
    }
  }

  // Tackster reco->sim associations
  MergeTracksters_recoToSim_SC.resize(trackstersmerged.size());
  MergeTracksters_recoToSim_SC_score.resize(trackstersmerged.size());
  MergeTracksters_recoToSim_SC_sharedE.resize(trackstersmerged.size());
  for (size_t i = 0; i < trackstersmerged.size(); ++i) {
    const edm::Ref<ticl::TracksterCollection> tsRef(tracksters_merged_h, i);

    // CLUE3D -> STS-SC
    const auto stsSC_iter = MergetsRecoSimSCMap.find(tsRef);
    if (stsSC_iter != MergetsRecoSimSCMap.end()) {
      const auto& stsSCassociated = stsSC_iter->val;
      for (auto& sts : stsSCassociated) {
        auto sts_id = (sts.first).get() - (edm::Ref<ticl::TracksterCollection>(
                                               simTrackstersSC_h, 0)).get();
        MergeTracksters_recoToSim_SC[i].push_back(sts_id);
        MergeTracksters_recoToSim_SC_score[i].push_back(sts.second.second);
        MergeTracksters_recoToSim_SC_sharedE[i].push_back(sts.second.first);
      }
    }
  }

  // SimTracksters
  nsimTrackstersSC = simTrackstersSC.size();
  MergeTracksters_simToReco_SC.resize(nsimTrackstersSC);
  MergeTracksters_simToReco_SC_score.resize(nsimTrackstersSC);
  MergeTracksters_simToReco_SC_sharedE.resize(nsimTrackstersSC);
  for (size_t i = 0; i < nsimTrackstersSC; ++i) {
    const edm::Ref<ticl::TracksterCollection> stsSCRef(simTrackstersSC_h, i);

    // STS-SC -> CLUE3D
    const auto ts_iter = MergetsSimToRecoSCMap.find(stsSCRef);
    if (ts_iter != MergetsSimToRecoSCMap.end()) {
      const auto& tsAssociated = ts_iter->val;
      for (auto& ts : tsAssociated) {
        auto ts_idx = (ts.first).get() - (edm::Ref<ticl::TracksterCollection>(
                                              tracksters_merged_h, 0)).get();
        MergeTracksters_simToReco_SC[i].push_back(ts_idx);
        MergeTracksters_simToReco_SC_score[i].push_back(ts.second.second);
        MergeTracksters_simToReco_SC_sharedE[i].push_back(ts.second.first);
      }
    }
  }

  // Tackster reco->sim associations
  MergeTracksters_recoToSim_CP.resize(trackstersmerged.size());
  MergeTracksters_recoToSim_CP_score.resize(trackstersmerged.size());
  MergeTracksters_recoToSim_CP_sharedE.resize(trackstersmerged.size());
  for (size_t i = 0; i < trackstersmerged.size(); ++i) {
    const edm::Ref<ticl::TracksterCollection> tsRef(tracksters_merged_h, i);

    // CLUE3D -> STS-CP
    const auto stsCP_iter = MergetsRecoSimCPMap.find(tsRef);
    if (stsCP_iter != MergetsRecoSimCPMap.end()) {
      const auto& stsCPassociated = stsCP_iter->val;
      for (auto& sts : stsCPassociated) {
        auto sts_id = (sts.first).get() - (edm::Ref<ticl::TracksterCollection>(
                                               simTrackstersCP_h, 0)).get();
        MergeTracksters_recoToSim_CP[i].push_back(sts_id);
        MergeTracksters_recoToSim_CP_score[i].push_back(sts.second.second);
        MergeTracksters_recoToSim_CP_sharedE[i].push_back(sts.second.first);
      }
    }
  }

  // SimTracksters
  nsimTrackstersCP = simTrackstersCP.size();
  MergeTracksters_simToReco_CP.resize(nsimTrackstersCP);
  MergeTracksters_simToReco_CP_score.resize(nsimTrackstersCP);
  MergeTracksters_simToReco_CP_sharedE.resize(nsimTrackstersCP);
  for (size_t i = 0; i < nsimTrackstersCP; ++i) {
    const edm::Ref<ticl::TracksterCollection> stsCPRef(simTrackstersCP_h, i);

    // STS-CP -> CLUE3D
    const auto ts_iter = MergetsSimToRecoCPMap.find(stsCPRef);
    if (ts_iter != MergetsSimToRecoCPMap.end()) {
      const auto& tsAssociated = ts_iter->val;
      for (auto& ts : tsAssociated) {
        auto ts_idx = (ts.first).get() - (edm::Ref<ticl::TracksterCollection>(
                                              tracksters_merged_h, 0)).get();
        MergeTracksters_simToReco_CP[i].push_back(ts_idx);
        MergeTracksters_simToReco_CP_score[i].push_back(ts.second.second);
        MergeTracksters_simToReco_CP_sharedE[i].push_back(ts.second.first);
      }
    }
  }

  // Tracks
  for (size_t i = 0; i < tracks.size(); i++) {
    auto track = tracks[i];
    reco::TrackRef trackref = reco::TrackRef(tracks_h, i);
    int iSide = int(track.eta() > 0);
    const auto& fts =
        trajectoryStateTransform::outerFreeState((track), bFieldProd);
    // to the HGCal front
    const auto& tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
    if (tsos.isValid()) {
      const auto& globalPos = tsos.globalPosition();
      const auto& globalMom = tsos.globalMomentum();
      track_ev.push_back(event_index);
      track_id.push_back(i);
      track_hgcal_x.push_back(globalPos.x());
      track_hgcal_y.push_back(globalPos.y());
      track_hgcal_z.push_back(globalPos.z());
      track_hgcal_eta.push_back(globalPos.eta());
      track_hgcal_phi.push_back(globalPos.phi());
      track_hgcal_px.push_back(globalMom.x());
      track_hgcal_py.push_back(globalMom.y());
      track_hgcal_pz.push_back(globalMom.z());
      track_pt.push_back(globalMom.perp());
      track_charge.push_back(track.charge());
      track_time.push_back(trackTime[trackref]);
      track_time_quality.push_back(trackTimeQual[trackref]);
      track_time_err.push_back(trackTimeErr[trackref]);
      track_nhits.push_back(tracks[i].recHitsSize());
    }
  }

  /* MTD VALIDATION */

  auto simTrackToRecoTrack = [&](UniqueSimTrackId simTkId)
      ->std::pair<int, float> {
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
            trackIdx = &(*associatedRecoTracks[0].first) - &tracks[0];
            quality = associatedRecoTracks[0].second;
          }
        }
      }
    }
    return {trackIdx, quality};
  };

  auto MtdSCtoRC = [&] (MtdSimLayerCluster sc) -> int {
    auto hAndF = sc.hits_and_fractions();
    int index = -1;
    for (const auto& cluster : MtdRecoClusters) {
      index++;
      MTDDetId cluId = cluster.id();
      for (int ihit = 0; ihit < cluster.size(); ++ihit) {
        int hit_row = cluster.minHitRow() + cluster.hitOffset()[ihit * 2];
        int hit_col = cluster.minHitCol() + cluster.hitOffset()[ihit * 2 + 1];
        for (const auto& RecHitsColl : RecHitsHandle) {
          for (const auto& recHit : *RecHitsColl) {
            MTDDetId hitId(recHit.id().rawId());

            if (hitId.mtdSide() != cluId.mtdSide() || hitId.mtdRR() != cluId.mtdRR() || recHit.row() != hit_row || recHit.column() != hit_col)
              continue;

            if (recHit.energy() != cluster.hitENERGY()[ihit] || recHit.time() != cluster.hitTIME()[ihit])
              continue;

            uint64_t targetId = static_cast<uint64_t>(hitId.rawId()) << 32;
            targetId |= recHit.row() << 16;
            targetId |= recHit.column();
            auto found_clu = std::find_if( hAndF.begin(), hAndF.end(),
                [&targetId] (const std::pair<uint64_t, float>& elem) { return elem.first == targetId;} );
            if (found_clu != hAndF.end()) {
              return index;
            }
          }
        }
      }
    }
    return -1;
  };

  [[maybe_unused]] auto MtdRCtoSC = [&] (FTLCluster cluster) -> int {
      MTDDetId cluId = cluster.id();
      std::vector<uint64_t> rcDetIds;
      for (int ihit = 0; ihit < cluster.size(); ++ihit) {
        int hit_row = cluster.minHitRow() + cluster.hitOffset()[ihit * 2];
        int hit_col = cluster.minHitCol() + cluster.hitOffset()[ihit * 2 + 1];
        for (const auto& RecHitsColl : RecHitsHandle) {
          for (const auto& recHit : *RecHitsColl) {
            MTDDetId hitId(recHit.id().rawId());

            if (hitId.mtdSide() != cluId.mtdSide() || hitId.mtdRR() != cluId.mtdRR() || recHit.row() != hit_row || recHit.column() != hit_col)
              continue;

            if (recHit.energy() != cluster.hitENERGY()[ihit] || recHit.time() != cluster.hitTIME()[ihit])
              continue;

            uint64_t targetId = static_cast<uint64_t>(hitId.rawId()) << 32;
            targetId |= recHit.row() << 16;
            targetId |= recHit.column();

	    rcDetIds.push_back(targetId);
          }
        }
      }
      int index = -1;
      for (const auto& sc : mtdSimLCs) {
        index ++;
        const auto& hAndF = sc.hits_and_fractions();
        std::vector<int> scHits(hAndF.size());
        std::transform(hAndF.begin(), hAndF.end(), scHits.begin(), [](const std::pair<int, float>& pair) {
          return pair.first;
        });
        std::vector<uint64_t> intersection;  
        std::set_intersection(scHits.begin(), scHits.end(), rcDetIds.begin(), rcDetIds.end(), std::back_inserter(intersection));
    
   	if (!intersection.empty())
          return index;
      }
      return -1;
    };

  std::cout << "BEGIN\nThere are " << mtdSimTracksters.size() << " mtdSimTracksters"<< std::endl;
number_of_mtdSimTracksters_ = mtdSimTracksters.size();
  // loop on simTracksters
  for (const auto& st : mtdSimTracksters) {
    const auto& SCindices = st.clusters();
number_of_mtdSimLCinST_.push_back(SCindices.size());
    std::cout << "  number of LC in ST: " << SCindices.size() << std::endl;

    const auto& simTrack = st.g4Tracks()[0];
    int trackIndex = -1;
   std::cout << "  ST associated with simTrack " << simTrack.trackId();
    // association with recotrack
    UniqueSimTrackId simTkIds(simTrack.trackId(), simTrack.eventId());
    auto bestAssociatedRecoTrack = simTrackToRecoTrack(simTkIds);
std::vector<int> mtdHit_det;
      std::vector<const MTDTrackingRecHit*> mtdRecHits;
    if (bestAssociatedRecoTrack.first != -1 and bestAssociatedRecoTrack.second > 0.75f) {
      trackIndex = bestAssociatedRecoTrack.first;
   std::cout << " and associated with recoTrack " << trackIndex << std::endl;
      const reco::TrackRef trackref(tracks_h, trackIndex);
      if (trackAssoc[trackref] == -1) {
        std::cout << " ** Extended track not associated ** \n";
number_of_mtdHits_.push_back(0);
mtdHit_det.push_back({-1});
        continue;
      }

      const reco::TrackRef mtdTrackref = reco::TrackRef(mtdTracks_h, trackAssoc[trackref]);
      const reco::Track& track = *mtdTrackref;
std::cout << " looping on track hits, there are ";
      for (const auto hit : track.recHits()) {
        if (hit->isValid() == false)
          continue;
        MTDDetId Hit = hit->geographicalId();
        if ((Hit.det() == 6) && (Hit.subdetId() == 1)) {
          const MTDTrackingRecHit* mtdhit = static_cast<const MTDTrackingRecHit*>(hit);
          mtdRecHits.push_back(mtdhit);
    mtdHit_det.push_back(Hit.mtdSubDetector());
	}
      }
      number_of_mtdHits_.push_back(mtdRecHits.size());
std::cout << mtdRecHits.size() <<  " hits in MTD\n";
    } else {
number_of_mtdHits_.push_back(0);
mtdHit_det.push_back({-1});
	std::cout << std::endl;
    }
mtdHits_det_.push_back(mtdHit_det);
ST_simTrack_.push_back(simTrack.trackId());
ST_recoTrack_.push_back(trackIndex);
simTrack_phi_.push_back(st.phi());
simTrack_eta_.push_back(st.eta());
simTrack_pt_.push_back(st.pt());

std::vector<bool> is_looper;
std::vector<int> simLC_idx;
std::vector<int> LC_simTrack;
std::vector<int> LC_recoTrack;
std::vector<float> simLC_time;
std::vector<float> simLC_posX;
std::vector<float> simLC_posY;
std::vector<float> recocluster_time;
std::vector<float> recocluster_timeErr;
std::vector<float> recocluster_posX;
std::vector<float> recocluster_posY;
std::vector<bool> simLC_matched;
std::vector<bool> simLC_CorrectMatch;
std::vector<bool> simLC_DirectMatch;

    for (const auto& idx : SCindices) {
    const auto& sc = mtdSimLCs[idx];
      std::cout << "Considering LC " << idx << " with time " << sc.simLCTime();

simLC_idx.push_back(idx);
LC_simTrack.push_back(simTrack.trackId());
LC_recoTrack.push_back(trackIndex);
simLC_time.push_back(sc.simLCTime());
simLC_posX.push_back(sc.simLCPos().x());
simLC_posY.push_back(sc.simLCPos().y());

      const auto& hAndF = sc.hits_and_fractions();
      // flag per dire se e' looper o no in qualche modo
      [[maybe_unused]] bool isLooper = false;
      if (std::abs(sc.simLCTime() - st.time()) > 0.2)
         isLooper = true; 
std::cout << (isLooper ? " marked as looper\n" : "\n"); 
is_looper.push_back(isLooper);

std::cout << "There are " << mtdRecHits.size() << " mtd hits\n"; 
        bool found = false;
        for (const auto& mtdhit : mtdRecHits){
            const auto& cluster = mtdhit->mtdCluster();
            MTDDetId cluId = cluster.id();
            for (int ihit = 0; ihit < cluster.size(); ++ihit) {
              int hit_row = cluster.minHitRow() + cluster.hitOffset()[ihit * 2];
              int hit_col = cluster.minHitCol() + cluster.hitOffset()[ihit * 2 + 1];
              for (const auto& RecHitsColl : RecHitsHandle) {
                for (const auto& recHit : *RecHitsColl) {
                  MTDDetId hitId(recHit.id().rawId());

                  if (hitId.mtdSide() != cluId.mtdSide() || hitId.mtdRR() != cluId.mtdRR() || recHit.row() != hit_row || recHit.column() != hit_col)
                    continue;

                  if (recHit.energy() != cluster.hitENERGY()[ihit] || recHit.time() != cluster.hitTIME()[ihit])
                    continue;

                  uint64_t targetId = static_cast<uint64_t>(hitId.rawId()) << 32;
                  targetId |= recHit.row() << 16;
                  targetId |= recHit.column();
                  auto found_clu = std::find_if( hAndF.begin(), hAndF.end(),
                      [&targetId](const std::pair<uint64_t, float>& elem){ return elem.first == targetId;} );
                    if (found_clu != hAndF.end()) {
                      found = true;
std::cout << "simLC matched with reco cluster with time " << mtdhit->time();
recocluster_time.push_back(mtdhit->time());
recocluster_timeErr.push_back(mtdhit->timeError());
recocluster_posX.push_back(cluster.x());
recocluster_posY.push_back(cluster.y());
simLC_matched.push_back(true);
   int rc_directMatch = MtdSCtoRC(sc);
   if (rc_directMatch != -1) {
     FTLCluster rc = MtdRecoClusters[rc_directMatch];
     MTDDetId rcId = rc.id();
     if (cluId == rcId and rc.time() == cluster.time() and rc.energy() == cluster.energy() and rc.size() == cluster.size()) {
       simLC_CorrectMatch.push_back(true);
       simLC_DirectMatch.push_back(true);
std::cout << "direct match also with reco and correct\n";
     } else {
       simLC_DirectMatch.push_back(true);
       simLC_CorrectMatch.push_back(false);
std::cout << " direct match also with reco but not correct\n";
     }
   } else {
       simLC_CorrectMatch.push_back(false);
       simLC_DirectMatch.push_back(false); 
std::cout << "but not direct match with reco\n";
   }
                      break;
                    }
                  }
                  if (found)
		    break;
                }
                if (found)
		  break;
              } 
              if (found)
	        break;
            }
            if (not found){
recocluster_time.push_back(0);
recocluster_timeErr.push_back(-1);
recocluster_posX.push_back(-99.);
recocluster_posY.push_back(-99.);
simLC_matched.push_back(false);
std::cout << "simLC not matched with reco cluster\n";
   int rc_directMatch = MtdSCtoRC(sc);
   if (rc_directMatch != -1) {
     simLC_DirectMatch.push_back(true);
std::cout << " but there is a direct match\n";
   } else {
     simLC_DirectMatch.push_back(false);
std::cout << " and there is no direct match\n";
   }
   simLC_CorrectMatch.push_back(false);
            }
          }
// push_back nel vettorone
simLC_is_looper_.push_back(is_looper);
simLC_idx_.push_back(simLC_idx);
simLC_simTrack_.push_back(LC_simTrack);
simLC_recoTrack_.push_back(LC_recoTrack);
simLC_time_.push_back(simLC_time);
simLC_posX_.push_back(simLC_posX);
simLC_posY_.push_back(simLC_posY);
recocluster_time_.push_back(recocluster_time);
recocluster_timeErr_.push_back(recocluster_timeErr);
recocluster_posX_.push_back(recocluster_posX);
recocluster_posY_.push_back(recocluster_posY);
simLC_matched_.push_back(simLC_matched);
simLC_CorrectMatch_.push_back(simLC_CorrectMatch);
simLC_DirectMatch_.push_back(simLC_DirectMatch);
        }

std::cout << "------------------------------\n";

  node_linked_inners.resize(tracksters.size());
  node_linked_outers.resize(tracksters.size());
  isRootTrackster.resize(tracksters.size(), false);
  for (size_t i = 0; i < tracksters.size(); ++i) {
    const auto& node = graph.getNode((int)i);
    auto this_inners = node.getInner();
    auto this_outers = node.getOuter();
    node_linked_inners[i].insert(node_linked_inners[i].end(), this_inners.begin(), this_inners.end());
    node_linked_outers[i].insert(node_linked_outers[i].end(), this_outers.begin(), this_outers.end());
    if (node.getInner().empty())
      isRootTrackster[i] = true;
  }

  trackster_tree_->Fill();
  cluster_tree_->Fill();
  graph_tree_->Fill();
  candidate_tree_->Fill();
  tracksters_merged_tree_->Fill();
  associations_tree_->Fill();
  simtrackstersSC_tree_->Fill();
  simtrackstersCP_tree_->Fill();
  tracks_tree_->Fill();
  simTICLCandidate_tree->Fill();
  MTDclusters_tree->Fill();
}

void TICLDumper::endJob() {}

void TICLDumper::fillDescriptions(
    edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackstersclue3d",
                          edm::InputTag("ticlTrackstersCLUE3DHigh"));
  desc.add<edm::InputTag>("layerClusters",
                          edm::InputTag("hgcalMergeLayerClusters"));
  desc.add<edm::InputTag>(
      "layer_clustersTime",
      edm::InputTag("hgcalMergeLayerClusters", "timeLayerCluster"));
  desc.add<edm::InputTag>("ticlgraph", edm::InputTag("ticlGraph"));
  desc.add<edm::InputTag>("ticlcandidates",
                          edm::InputTag("ticlTrackstersMerge"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("tracksTime", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("tracksTimeQual",
                          edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<edm::InputTag>("tracksTimeErr", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("trackstersmerged",
                          edm::InputTag("ticlTrackstersMerge"));
  desc.add<edm::InputTag>("simtrackstersSC",
                          edm::InputTag("ticlSimTracksters"));
  desc.add<edm::InputTag>("simtrackstersCP",
                          edm::InputTag("ticlSimTracksters", "fromCPs"));
  desc.add<edm::InputTag>("simTICLCandidates",
                          edm::InputTag("ticlSimTracksters"));
  desc.add<edm::InputTag>(
      "recoToSimAssociatorSC",
      edm::InputTag("tracksterSimTracksterAssociationPRbyCLUE3D", "recoToSim"));
  desc.add<edm::InputTag>(
      "simToRecoAssociatorSC",
      edm::InputTag("tracksterSimTracksterAssociationPRbyCLUE3D", "simToReco"));
  desc.add<edm::InputTag>(
      "recoToSimAssociatorCP",
      edm::InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3D",
                    "recoToSim"));
  desc.add<edm::InputTag>(
      "simToRecoAssociatorCP",
      edm::InputTag("tracksterSimTracksterAssociationLinkingbyCLUE3D",
                    "simToReco"));
  desc.add<edm::InputTag>(
      "MergerecoToSimAssociatorSC",
      edm::InputTag("tracksterSimTracksterAssociationPR", "recoToSim"));
  desc.add<edm::InputTag>(
      "MergesimToRecoAssociatorSC",
      edm::InputTag("tracksterSimTracksterAssociationPR", "simToReco"));
  desc.add<edm::InputTag>(
      "MergerecoToSimAssociatorCP",
      edm::InputTag("tracksterSimTracksterAssociationLinking", "recoToSim"));
  desc.add<edm::InputTag>(
      "MergesimToRecoAssociatorCP",
      edm::InputTag("tracksterSimTracksterAssociationLinking", "simToReco"));
  desc.add<edm::InputTag>("simclusters",
                          edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("caloparticles",
                          edm::InputTag("mix", "MergedCaloTruth"));
  desc.add<edm::InputTag>("MtdSimTracksters",
                          edm::InputTag("mix", "MergedMtdTruthST"));
  desc.add<edm::InputTag>("MtdSimLayerCluster",
                          edm::InputTag("mix", "MergedMtdTruthLC"));
  desc.add<edm::InputTag>("tpToTrack",
                          edm::InputTag("trackingParticleRecoTrackAsssociation"));
  desc.add<edm::InputTag>("simTrackToTPMap",
                          edm::InputTag("simHitTPAssocProducer", "simTrackToTP"));
  desc.add<edm::InputTag>("tracksWithMtd",
                          edm::InputTag("trackExtenderWithMTD"));
  desc.add<edm::InputTag>("mtdtrackAssocSrc",
                          edm::InputTag("trackExtenderWithMTD:generalTrackassoc"));
  desc.add<edm::InputTag>("BTLrecHits",
                          edm::InputTag("mtdRecHits", "FTLBarrel"));
  desc.add<edm::InputTag>("ETLrecHits",
                          edm::InputTag("mtdRecHits", "FTLEndcap"));
  desc.add<edm::InputTag>("BTLclusters",
                          edm::InputTag("mtdClusters", "FTLBarrel"));
  desc.add<edm::InputTag>("ETLclusters",
                          edm::InputTag("mtdClusters", "FTLEndcap"));
  desc.add<std::string>("detector", "HGCAL");
  desc.add<std::string>("propagator", "PropagatorWithMaterial");
  descriptions.add("ticlDumper", desc);
}

DEFINE_FWK_MODULE(TICLDumper);
