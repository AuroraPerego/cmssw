import FWCore.ParameterSet.Config as cms

layerClusterToTracksterAssociation = cms.EDProducer("LCToTSAssociatorProducer",
    layer_clusters = cms.InputTag("hgcalMergeLayerClusters"),
    tracksters = cms.InputTag("ticlTracksters"),
)

from SimCalorimetry.HGCalAssociatorProducers.LCToTSAssociatorProducer_cfi import LCToTSAssociatorProducer

layerClusterToCLUE3DTracksterAssociation = LCToTSAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlTrackstersCLUE3DHigh")
)

layerClusterToTracksterMergeAssociation = LCToTSAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlTrackstersMerge")
)

layerClusterToTracksterSuperclusteringAssociation = LCToTSAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlTracksterLinksSuperclusteringDNN")
)

layerClusterToSimTracksterAssociation = LCToTSAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlSimTracksters")
)
layerClusterToSimTracksterAssociation2Hits = LCToTSAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlSimTracksters2Hits")
)

layerClusterToSimTracksterFromCPsAssociation = LCToTSAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlSimTracksters", "fromCPs")
)
layerClusterToSimTracksterFromCPsAssociation2Hits = LCToTSAssociatorProducer.clone(
    tracksters = cms.InputTag("ticlSimTracksters2Hits", "fromCPs")
)

from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5
from Configuration.ProcessModifiers.ticl_superclustering_mustache_ticl_cff import ticl_superclustering_mustache_ticl 
from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5
ticl_v5.toModify(layerClusterToTracksterMergeAssociation, tracksters = cms.InputTag("ticlCandidate"))
(ticl_v5 & ticl_superclustering_mustache_ticl).toModify(layerClusterToTracksterSuperclusteringAssociation, tracksters = cms.InputTag("ticlTracksterLinksSuperclusteringMustache"))
