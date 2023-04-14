import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.simTrackstersProducer_cfi import simTrackstersProducer as _simTrackstersProducer
from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingGlobal, ticlSeedingGlobalHFNose
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer
from RecoHGCal.TICL.fineSimTrackstersProducer_cfi import fineSimTrackstersProducer as _fineSimTrackstersProducer

from SimCalorimetry.HGCalSimProducers.hgcHitAssociation_cfi import lcAssocByEnergyScoreProducer, scAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToCPAssociation_cfi import layerClusterCaloParticleAssociation as layerClusterCaloParticleAssociationProducer
from SimCalorimetry.HGCalAssociatorProducers.LCToSCAssociation_cfi import layerClusterSimClusterAssociation as layerClusterSimClusterAssociationProducer
from SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi import trackingParticleRecoTrackAsssociation
from SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi import simHitTPAssocProducer
# CA - PATTERN RECOGNITION


filteredLayerClustersSimTracksters = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSize",
    algo_number = 8,
    min_cluster_size = 0, # inclusive
    iteration_label = "ticlSimTracksters"
)

ticlSimTracksters = _simTrackstersProducer.clone(
		filtered_mask = "filteredLayerClustersSimTracksters:ticlSimTracksters"
)

ticlFineSimTracksters = _fineSimTrackstersProducer.clone(
  detector = "HGCAL",
  layer_clusters = "hgcalLayerClusters",
  time_layerclusters = "hgcalLayerClusters:timeLayerCluster",
  filtered_mask = "ticlSimTracksters",
  patternRecognitionBy = "CLUE3D",
    pluginPatternRecognitionByCLUE3D = dict (
        criticalDensity = 0.6,
        criticalEtaPhiDistance = 0.025,
    ),
  patternRecognitionMIPBy = "CA",
        pluginPatternRecognitionMIPByCA = dict (
        skip_layers = 3,
        min_layers_per_trackster = 0,
        min_cos_theta = 0.50, # ~10 degrees
        min_cos_pointing = 0.50,
        out_in_dfs = False,
        algo_verbosity = 4,
        max_delta_time = -1
    )
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(ticlSimTracksters,
    simclusters = "mixData:MergedCaloTruth",
    caloparticles = "mixData:MergedCaloTruth",
)

premix_stage2.toModify(ticlFineSimTracksters,
    simclusters = "mixData:MergedCaloTruth",
    caloparticles = "mixData:MergedCaloTruth",
)

ticlSimTrackstersTask = cms.Task(simHitTPAssocProducer, trackingParticleRecoTrackAsssociation, lcAssocByEnergyScoreProducer, layerClusterCaloParticleAssociationProducer, scAssocByEnergyScoreProducer, layerClusterSimClusterAssociationProducer, ticlSeedingGlobal, filteredLayerClustersSimTracksters, ticlSimTracksters, ticlFineSimTracksters)
