import FWCore.ParameterSet.Config as cms

from RecoTICL.LayerClustersProducers.ticlLayerTileProducer_cfi import ticlLayerTileProducer

from RecoTICL.TrackstersProducers.CLUE3DEM_cff import *
from RecoTICL.TrackstersProducers.CLUE3DHAD_cff import *
from RecoTICL.CandidatesProducers.pfTICLProducerV5_cfi import pfTICLProducerV5 as _pfTICLProducerV5

from RecoTICL.LayerClustersProducers.ticlLayerTileProducer_cfi import ticlLayerTileProducer
from RecoTICL.TrackstersProducers.tracksterSelectionTf_cfi import *

from RecoTICL.LinkingProducers.tracksterLinksProducer_cfi import tracksterLinksProducer as _tracksterLinksProducer
from RecoTICL.CandidatesProducers.ticlCandidateProducer_cfi import ticlCandidateProducer as _ticlCandidateProducer
from RecoTICL.Configuration.RecoTICL_EventContent_cff import customiseForTICLv5EventContent
from RecoTICL.CandidatesProducers.iterativeTICL_cff import ticlIterLabels, ticlIterLabelsMerge
from RecoTICL.CandidatesProducers.ticlDumper_cfi import ticlDumper
from RecoTICL.LinkingProducers.mergedTrackstersProducer_cfi import mergedTrackstersProducer as _mergedTrackstersProducer
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import tracksterSimTracksterAssociationLinkingbyCLUE3D as _tracksterSimTracksterAssociationLinkingbyCLUE3D
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import tracksterSimTracksterAssociationPRbyCLUE3D  as _tracksterSimTracksterAssociationPRbyCLUE3D
from Validation.HGCalValidation.HGCalValidator_cff import hgcalValidator
from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import HGCalUncalibRecHit
from RecoTICL.TrackstersProducers.SimTracksters_cff import ticlSimTracksters

from RecoTICL.TrackstersProducers.FastJetStep_cff import ticlTrackstersFastJet
from RecoTICL.TrackstersProducers.EMStep_cff import ticlTrackstersEM, ticlTrackstersHFNoseEM
from RecoTICL.TrackstersProducers.TrkStep_cff import ticlTrackstersTrk, ticlTrackstersHFNoseTrk
from RecoTICL.TrackstersProducers.MIPStep_cff import ticlTrackstersMIP, ticlTrackstersHFNoseMIP
from RecoTICL.TrackstersProducers.HADStep_cff import ticlTrackstersHAD, ticlTrackstersHFNoseHAD
from RecoTICL.TrackstersProducers.CLUE3DEM_cff import ticlTrackstersCLUE3DEM
from RecoTICL.TrackstersProducers.CLUE3DHAD_cff import ticlTrackstersCLUE3DHAD
from RecoTICL.TrackstersProducers.CLUE3DHighStep_cff import ticlTrackstersCLUE3DHigh
from RecoTICL.TrackstersProducers.TrkEMStep_cff import ticlTrackstersTrkEM, filteredLayerClustersHFNoseTrkEM

from RecoTICL.CandidatesProducers.mtdSoAProducer_cfi import mtdSoAProducer as _mtdSoAProducer

def customiseTICLv5FromReco(process, enableDumper = False):
    # TensorFlow ESSource

    process.TFESSource = cms.Task(process.trackdnn_source)

    process.hgcalLayerClustersTask = cms.Task(process.hgcalLayerClustersEE,
                                              process.hgcalLayerClustersHSi,
                                              process.hgcalLayerClustersHSci,
                                              process.hgcalMergeLayerClusters)

    # Reconstruction

    process.ticlSimTracksters.computeLocalTime = cms.bool(True)

    process.ticlTrackstersCLUE3DHigh.pluginPatternRecognitionByCLUE3D.computeLocalTime = cms.bool(True)

    '''for future CLUE3D separate iterations
    process.ticlTrackstersCLUE3DHAD.pluginPatternRecognitionByCLUE3D.computeLocalTime = cms.bool(True)
    process.ticlTrackstersCLUE3DEM.pluginPatternRecognitionByCLUE3D.computeLocalTime = cms.bool(True)
    '''

    process.ticlLayerTileTask = cms.Task(ticlLayerTileProducer)

    process.ticlIterationsTask = cms.Task(
        process.ticlTrackstersCLUE3DHigh,
    )

    process.mtdSoA = _mtdSoAProducer.clone()
    process.mtdSoATask = cms.Task(process.mtdSoA)

    process.ticlTracksterLinks = _tracksterLinksProducer.clone()
    process.ticlTracksterLinks = _tracksterLinksProducer.clone(
            tracksters_collections = cms.VInputTag(
              'ticlTrackstersCLUE3DHigh'
            ),
    )

    process.ticlCandidate = _ticlCandidateProducer.clone()
    process.ticlCandidateTask = cms.Task(process.ticlCandidate)

    process.allTrackstersToSimTrackstersAssociationsByLCs = _allTrackstersToSimTrackstersAssociationsByLCs.clone()

    process.allTrackstersToSimTrackstersAssociationsByHits = _allTrackstersToSimTrackstersAssociationsByHits.clone()

    process.iterTICLTask = cms.Path(process.hgcalLayerClustersTask,
                            process.TFESSource,
                            process.ticlLayerTileTask,
                            process.mtdSoATask,
                            process.ticlIterationsTask,
                            process.ticlTracksterLinksTask,
                            process.ticlCandidateTask)

    process.particleFlowClusterHGCal.initialClusteringStep.tracksterSrc = "ticlCandidate"
    process.globalrecoTask.remove(process.ticlTrackstersMerge)


    process.mergeTICLTask = cms.Task()
    process.pfTICL = _pfTICLProducer.clone(
      ticlCandidateSrc = cms.InputTag('ticlCandidate'),
      isTICLv5 = cms.bool(True)
    )
    process.hgcalAssociators = cms.Task(process.recHitMapProducer, process.lcAssocByEnergyScoreProducer, process.layerClusterCaloParticleAssociationProducer,
                            process.scAssocByEnergyScoreProducer, process.layerClusterSimClusterAssociationProducer,
                            # FP 07/2024 new associators:
                            process.allLayerClusterToTracksterAssociations, process.allHitToTracksterAssociations,
                            process.allTrackstersToSimTrackstersAssociationsByLCs, process.allTrackstersToSimTrackstersAssociationsByHits,
                            process.hitToSimClusterCaloParticleAssociator, process.SimClusterToCaloParticleAssociation,
                            )

    if(enableDumper):
        process.ticlDumper = ticlDumper
        process.TFileService = cms.Service("TFileService",
                                           fileName=cms.string("histo.root")
                                           )

    process.FEVTDEBUGHLToutput_step = cms.EndPath(process.ticlDumper)

    process.TICL_Validation = cms.Path(process.ticlSimTrackstersTask, process.hgcalAssociators)

# Schedule definition
    process.schedule = cms.Schedule(process.iterTICLTask,
                                    process.TICL_Validation,
                                    process.FEVTDEBUGHLToutput_step)
    process = customiseForTICLv5EventContent(process)

    return process
