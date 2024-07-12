# Reconstruction
from RecoHGCal.TICL.iterativeTICL_cff import *
from RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cff import hgcalLayerClustersEE, hgcalLayerClustersHSi, hgcalLayerClustersHSci
from RecoLocalCalo.HGCalRecProducers.hgcalMergeLayerClusters_cfi import hgcalMergeLayerClusters
from RecoHGCal.TICL.ticlDumper_cfi import ticlDumper
# Validation
from Validation.HGCalValidation.HGCalValidator_cfi import *
from RecoLocalCalo.HGCalRecProducers.recHitMapProducer_cfi import recHitMapProducer

# Load DNN ESSource
from RecoTracker.IterativeTracking.iterativeTk_cff import trackdnn_source

# Automatic addition of the customisation function from RecoHGCal.Configuration.RecoHGCal_EventContent_cff
from RecoHGCal.Configuration.RecoHGCal_EventContent_cff import customiseHGCalOnlyEventContent
from SimCalorimetry.HGCalAssociatorProducers.simTracksterAssociatorByEnergyScore_cfi import simTracksterAssociatorByEnergyScore as simTsAssocByEnergyScoreProducer
from SimCalorimetry.HGCalAssociatorProducers.TSToSimTSAssociation_cfi import tracksterSimTracksterAssociationLinking, tracksterSimTracksterAssociationPR, tracksterSimTracksterAssociationLinkingbyCLUE3D, tracksterSimTracksterAssociationPRbyCLUE3D, tracksterSimTracksterAssociationLinkingPU, tracksterSimTracksterAssociationPRPU


def customiseTICLFromReco(process):
    # TensorFlow ESSource
    process.TFESSource = cms.Task(process.trackdnn_source)

    process.hgcalLayerClustersTask = cms.Task(process.hgcalLayerClustersEE,
                                              process.hgcalLayerClustersHSi,
                                              process.hgcalLayerClustersHSci,
                                              process.hgcalMergeLayerClusters)

# Reconstruction
    process.TICL = cms.Path(process.hgcalLayerClustersTask,
                            process.TFESSource,
                            process.ticlLayerTileTask,
                            process.ticlIterationsTask,
                            process.ticlTracksterMergeTask)
# Validation
    process.TICL_ValidationProducers = cms.Task(process.recHitMapProducer,
                                                process.lcAssocByEnergyScoreProducer,
                                                process.layerClusterCaloParticleAssociationProducer,
                                                process.scAssocByEnergyScoreProducer,
                                                process.layerClusterSimClusterAssociationProducer,
                                                process.simTsAssocByEnergyScoreProducer,
                                                process.simTracksterHitLCAssociatorByEnergyScoreProducer,
                                                process.tracksterSimTracksterAssociationLinking,
                                                process.tracksterSimTracksterAssociationPR,
                                                process.tracksterSimTracksterAssociationLinkingbyCLUE3D,
                                                process.tracksterSimTracksterAssociationPRbyCLUE3D,
                                                process.tracksterSimTracksterAssociationLinkingPU,
                                                process.tracksterSimTracksterAssociationPRPU
                                                )

    process.TICL_Validator = cms.Task(process.hgcalValidator)
    process.TICL_Validation = cms.Path(process.TICL_ValidationProducers,
                                       process.TICL_Validator
                                       )
# Path and EndPath definitions
    process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
    process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
    process.schedule = cms.Schedule(process.TICL,
                                    process.TICL_Validation,
                                    process.FEVTDEBUGHLToutput_step,
                                    process.DQMoutput_step)
# call to customisation function customiseHGCalOnlyEventContent imported from RecoHGCal.Configuration.RecoHGCal_EventContent_cff
    process = customiseHGCalOnlyEventContent(process)

    return process

def customiseTICLForDumper(process):

    process.ticlDumper = ticlDumper.clone(
        saveLCs=True,
        saveCLUE3DTracksters=True,
        saveTrackstersMerged=True,
        saveSimTrackstersSC=True,
        saveSimTrackstersCP=True,
        saveTICLCandidate=True,
        saveSimTICLCandidate=True,
        saveTracks=True,
        saveAssociations=True,
    )

    from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5
    ticl_v5.toModify(process.ticlDumper,
                     # trackstersclue3d = cms.InputTag('mergedTrackstersProducer'), # For future separate iterations
                     trackstersclue3d = cms.InputTag('ticlTrackstersCLUE3DHigh'),
                     ticlcandidates = cms.InputTag("ticlCandidate"),
                     trackstersmerged = cms.InputTag("ticlCandidate"),
                     trackstersInCand = cms.InputTag("ticlCandidate"))

    process.TFileService = cms.Service("TFileService",
                                       fileName=cms.string("histo.root")
                                       )
    process.FEVTDEBUGHLToutput_step = cms.EndPath(
        process.FEVTDEBUGHLToutput + process.ticlDumper)
    return process

def customiseTICLForDumperHLT(process):

    process.lcAssocByEnergyScoreProducerHLT = cms.EDProducer("HGCalLCToCPAssociatorByEnergyScoreProducer",
        hardScatterOnly = cms.bool(True),
        hitMapTag = cms.InputTag("recHitMapProducerHLT","hgcalRecHitMap"),
        hits = cms.VInputTag("HGCalRecHit:HGCEERecHits:HLT", "HGCalRecHit:HGCHEFRecHits:HLT", "HGCalRecHit:HGCHEBRecHits:HLT"),
        mightGet = cms.optional.untracked.vstring
    )

    process.layerClusterCaloParticleAssociationProducerHLT = cms.EDProducer("LCToCPAssociatorEDProducer",
        associator = cms.InputTag("lcAssocByEnergyScoreProducerHLT"),
        label_cp = cms.InputTag("mix","MergedCaloTruth"),
        label_lc = cms.InputTag("hgcalMergeLayerClusters", '', 'HLT')
    )

    process.scAssocByEnergyScoreProducerHLT = cms.EDProducer("HGCalLCToSCAssociatorByEnergyScoreProducer",
        hardScatterOnly = cms.bool(True),
        hitMapTag = cms.InputTag("recHitMapProducerHLT","hgcalRecHitMap"),
        hits = cms.VInputTag("HGCalRecHit:HGCEERecHits:HLT", "HGCalRecHit:HGCHEFRecHits:HLT", "HGCalRecHit:HGCHEBRecHits:HLT"),
        mightGet = cms.optional.untracked.vstring
    )

    process.layerClusterSimClusterAssociationProducerHLT = cms.EDProducer("LCToSCAssociatorEDProducer",
        associator = cms.InputTag("scAssocByEnergyScoreProducerHLT"),
        label_lcl = cms.InputTag("hgcalMergeLayerClusters", '', 'HLT'),
        label_scl = cms.InputTag("mix","MergedCaloTruth")
    )



    process.tpClusterProducerHLT = cms.EDProducer("ClusterTPAssociationProducer",
        mightGet = cms.optional.untracked.vstring,
        phase2OTClusterSrc = cms.InputTag("siPhase2Clusters", "", "HLT"),
        phase2OTSimLinkSrc = cms.InputTag("simSiPixelDigis","Tracker", "HLT"),
        pixelClusterSrc = cms.InputTag("siPixelClusters", "", "HLT"),
        pixelSimLinkSrc = cms.InputTag("simSiPixelDigis","Pixel", "HLT"),
        simTrackSrc = cms.InputTag("g4SimHits"),
        stripClusterSrc = cms.InputTag("siStripClusters"),
        stripSimLinkSrc = cms.InputTag("simSiStripDigis"),
        throwOnMissingCollections = cms.bool(True),
        trackingParticleSrc = cms.InputTag("mix","MergedTrackTruth")
    )

    process.quickTrackAssociatorByHitsHLT = cms.EDProducer("QuickTrackAssociatorByHitsProducer",
        AbsoluteNumberOfHits = cms.bool(False),
        Cut_RecoToSim = cms.double(0.75),
        PixelHitWeight = cms.double(1.0),
        Purity_SimToReco = cms.double(0.75),
        Quality_SimToReco = cms.double(0.5),
        SimToRecoDenominator = cms.string('reco'),
        ThreeHitTracksAreSpecial = cms.bool(True),
        cluster2TPSrc = cms.InputTag("tpClusterProducerHLT"),
        useClusterTPAssociation = cms.bool(True)
    )



    process.trackingParticleRecoTrackAsssociationHLT = cms.EDProducer("TrackAssociatorEDProducer",
        associator = cms.InputTag("quickTrackAssociatorByHitsHLT"),
        ignoremissingtrackcollection = cms.untracked.bool(False),
        label_tp = cms.InputTag("mix","MergedTrackTruth"),
        label_tr = cms.InputTag("generalTracks", '', 'HLT')
    )

    process.filteredLayerClustersSimTrackstersHLT = cms.EDProducer("FilteredLayerClustersProducer",
        LayerClusters = cms.InputTag("hgcalMergeLayerClusters", '', 'HLT'),
        LayerClustersInputMask = cms.InputTag("hgcalMergeLayerClusters","InitialLayerClustersMask", 'HLT'),
        algo_number = cms.vint32(6, 7, 8),
        clusterFilter = cms.string('ClusterFilterByAlgoAndSize'),
        iteration_label = cms.string('ticlSimTrackstersHLT'),
        max_cluster_size = cms.int32(9999),
        max_layerId = cms.int32(9999),
        mightGet = cms.optional.untracked.vstring,
        min_cluster_size = cms.int32(0),
        min_layerId = cms.int32(0)
    )

    process.ticlSimTrackstersHLT = cms.EDProducer("SimTrackstersProducer",
        MtdSimTracksters = cms.InputTag("mix","MergedMtdTruthST"),
        caloparticles = cms.InputTag("mix","MergedCaloTruth"),
        computeLocalTime = cms.bool(True),
        cutTk = cms.string('1.48 < abs(eta) < 3.0 && pt > 1. && quality("highPurity") && hitPattern().numberOfLostHits("MISSING_OUTER_HITS") < 5'),
        detector = cms.string('HGCAL'),
        filtered_mask = cms.InputTag("filteredLayerClustersSimTrackstersHLT","ticlSimTrackstersHLT"),
        fractionCut = cms.double(0),
        layerClusterCaloParticleAssociator = cms.InputTag("layerClusterCaloParticleAssociationProducerHLT"),
        layerClusterSimClusterAssociator = cms.InputTag("layerClusterSimClusterAssociationProducerHLT"),
        layer_clusters = cms.InputTag("hgcalMergeLayerClusters", '', 'HLT'),
        mightGet = cms.optional.untracked.vstring,
        qualityCutTrack = cms.double(0.75),
        recoTracks = cms.InputTag("generalTracks", '', 'HLT'),
        simTrackToTPMap = cms.InputTag("simHitTPAssocProducer","simTrackToTP"),
        simVertices = cms.InputTag("g4SimHits"),
        simclusters = cms.InputTag("mix","MergedCaloTruth"),
        time_layerclusters = cms.InputTag("hgcalMergeLayerClusters","timeLayerCluster", "HLT"),
        tpToTrack = cms.InputTag("trackingParticleRecoTrackAsssociationHLT"),
        trackingParticles = cms.InputTag("mix","MergedTrackTruth")
    )

    process.recHitMapProducerHLT = cms.EDProducer('RecHitMapProducer',
        EEInput = cms.InputTag('HGCalRecHit', 'HGCEERecHits', 'HLT'),
        FHInput = cms.InputTag('HGCalRecHit', 'HGCHEFRecHits', 'HLT'),
        BHInput = cms.InputTag('HGCalRecHit', 'HGCHEBRecHits', 'HLT'),
        EBInput = cms.InputTag('particleFlowRecHitECAL', 'Cleaned', 'HLT'),
        HBInput = cms.InputTag('particleFlowRecHitHBHE', 'Cleaned', 'HLT'),
        HOInput = cms.InputTag('particleFlowRecHitHO', 'Cleaned', 'HLT'),
        hgcalOnly = cms.bool(True),
        mightGet = cms.optional.untracked.vstring
    )

    process.tracksterSimTracksterAssociationLinkingHLT = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
        associator = cms.InputTag("simTracksterHitLCAssociatorByEnergyScoreProducer"),
        label_cp = cms.InputTag("mix","MergedCaloTruth"),
        label_lcl = cms.InputTag("hgcalMergeLayerClusters", '', 'HLT'),
        label_scl = cms.InputTag("mix","MergedCaloTruth"),
        label_simTst = cms.InputTag("ticlSimTrackstersHLT","fromCPs"),
        label_tst = cms.InputTag('hltTiclCandidate')
    )

    process.tracksterSimTracksterAssociationPRHLT = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
        associator = cms.InputTag("simTracksterHitLCAssociatorByEnergyScoreProducer"),
        label_cp = cms.InputTag("mix","MergedCaloTruth"),
        label_lcl = cms.InputTag("hgcalMergeLayerClusters", "", "HLT"),
        label_scl = cms.InputTag("mix","MergedCaloTruth"),
        label_simTst = cms.InputTag("ticlSimTrackstersHLT"),
        label_tst = cms.InputTag("hltTiclCandidate")
    )

    process.tracksterSimTracksterAssociationLinkingbyCLUE3DHLT = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
        associator = cms.InputTag("simTracksterHitLCAssociatorByEnergyScoreProducer"),
        label_cp = cms.InputTag("mix","MergedCaloTruth"),
        label_lcl = cms.InputTag("hgcalMergeLayerClusters", "", "HLT"),
        label_scl = cms.InputTag("mix","MergedCaloTruth"),
        label_simTst = cms.InputTag("ticlSimTrackstersHLT", "fromCPs"),
        label_tst = cms.InputTag("ticlTrackstersCLUE3DHigh", "", "HLT")
    )

    process.tracksterSimTracksterAssociationPRbyCLUE3DHLT = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
        associator = cms.InputTag("simTracksterHitLCAssociatorByEnergyScoreProducer"),
        label_cp = cms.InputTag("mix","MergedCaloTruth"),
        label_lcl = cms.InputTag("hgcalMergeLayerClusters", "", "HLT"),
        label_scl = cms.InputTag("mix","MergedCaloTruth"),
        label_simTst = cms.InputTag("ticlSimTrackstersHLT"),
        label_tst = cms.InputTag("ticlTrackstersCLUE3DHigh", "", "HLT")
    )

    process.tracksterSimTracksterAssociationLinkingPUHLT = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
        associator = cms.InputTag("simTracksterHitLCAssociatorByEnergyScoreProducer"),
        label_cp = cms.InputTag("mix","MergedCaloTruth"),
        label_lcl = cms.InputTag("hgcalMergeLayerClusters", "", "HLT"),
        label_scl = cms.InputTag("mix","MergedCaloTruth"),
        label_simTst = cms.InputTag("ticlSimTrackstersHLT","PU"),
        label_tst = cms.InputTag("hltTiclCandidate")
    )

    process.tracksterSimTracksterAssociationPRPUHLT = cms.EDProducer("TSToSimTSHitLCAssociatorEDProducer",
        associator = cms.InputTag("simTracksterHitLCAssociatorByEnergyScoreProducer"),
        label_cp = cms.InputTag("mix","MergedCaloTruth"),
        label_lcl = cms.InputTag("hgcalMergeLayerClusters", "", "HLT"),
        label_scl = cms.InputTag("mix","MergedCaloTruth"),
        label_simTst = cms.InputTag("ticlSimTrackstersHLT","PU"),
        label_tst = cms.InputTag("hltTiclCandidate")
    )

    process.ticlDumperHLT = ticlDumper.clone(
       trackstersclue3d = cms.InputTag('ticlTrackstersCLUE3DHigh', '', 'HLT'),
       trackstersInCand = cms.InputTag('ticlTrackstersCLUE3DHigh', '', 'HLT'),
       layerClusters = cms.InputTag('hgcalMergeLayerClusters', '', 'HLT'),
       layer_clustersTime = cms.InputTag('hgcalMergeLayerClusters', 'timeLayerCluster', 'HLT'),
       ticlcandidates = cms.InputTag('ticlTrackstersMerge', '', 'HLT'),
       tracks = cms.InputTag('generalTracks', '', 'HLT'),
       #tracksTime = cms.InputTag('tofPID', 't0'),
       #tracksTimeQual = cms.InputTag('mtdTrackQualityMVA', 'mtdQualMVA'),
       #tracksTimeErr = cms.InputTag('tofPID', 'sigmat0'),
       #tracksBeta = cms.InputTag('trackExtenderWithMTD', 'generalTrackBeta'),
       #tracksTimeMtd = cms.InputTag('trackExtenderWithMTD', 'generalTracktmtd'),
       #tracksTimeMtdErr = cms.InputTag('trackExtenderWithMTD', 'generalTracksigmatmtd'),
       #tracksPosMtd = cms.InputTag('trackExtenderWithMTD', 'generalTrackmtdpos'),
       trackstersmerged = cms.InputTag('ticlTrackstersMerge', '', 'HLT'),
       muons = cms.InputTag('hltPhase2L3Muons'),
       simtrackstersSC = cms.InputTag('ticlSimTrackstersHLT'),
       simtrackstersCP = cms.InputTag('ticlSimTrackstersHLT', 'fromCPs'),
       simtrackstersPU = cms.InputTag('ticlSimTrackstersHLT', 'PU'),
       simTICLCandidates = cms.InputTag('ticlSimTrackstersHLT'),
       recoToSimAssociatorSC = cms.InputTag('tracksterSimTracksterAssociationPRbyCLUE3DHLT', 'recoToSim'),
       simToRecoAssociatorSC = cms.InputTag('tracksterSimTracksterAssociationPRbyCLUE3DHLT', 'simToReco'),
       recoToSimAssociatorCP = cms.InputTag('tracksterSimTracksterAssociationLinkingbyCLUE3DHLT', 'recoToSim'),
       simToRecoAssociatorCP = cms.InputTag('tracksterSimTracksterAssociationLinkingbyCLUE3DHLT', 'simToReco'),
       MergerecoToSimAssociatorSC = cms.InputTag('tracksterSimTracksterAssociationPRHLT', 'recoToSim'),
       MergesimToRecoAssociatorSC = cms.InputTag('tracksterSimTracksterAssociationPRHLT', 'simToReco'),
       MergerecoToSimAssociatorCP = cms.InputTag('tracksterSimTracksterAssociationLinkingHLT', 'recoToSim'),
       MergesimToRecoAssociatorCP = cms.InputTag('tracksterSimTracksterAssociationLinkingHLT', 'simToReco'),
       MergerecoToSimAssociatorPU = cms.InputTag('tracksterSimTracksterAssociationLinkingPUHLT', 'recoToSim'),
       MergesimToRecoAssociatorPU = cms.InputTag('tracksterSimTracksterAssociationLinkingPUHLT', 'simToReco'),
       simclusters = cms.InputTag('mix', 'MergedCaloTruth'),
       caloparticles = cms.InputTag('mix', 'MergedCaloTruth'),
       detector = cms.string('HGCAL'),
       propagator = cms.string('PropagatorWithMaterial'),
       useMTDtiming = cms.bool(False),
       saveLCs = cms.bool(True),
       saveCLUE3DTracksters = cms.bool(True),
       saveTrackstersMerged = cms.bool(True),
       saveSimTrackstersSC = cms.bool(True),
       saveSimTrackstersCP = cms.bool(True),
       saveTICLCandidate = cms.bool(True),
       saveSimTICLCandidate = cms.bool(True),
       saveTracks = cms.bool(True),
       saveAssociations = cms.bool(True),
    )

    from Configuration.ProcessModifiers.ticl_v5_cff import ticl_v5
    ticl_v5.toModify(process.ticlDumperHLT,
                     # trackstersclue3d = cms.InputTag('mergedTrackstersProducer'), # For future separate iterations
                     trackstersclue3d = cms.InputTag('ticlTrackstersCLUE3DHigh', '', 'HLT'),
                     ticlcandidates = cms.InputTag('hltTiclCandidate'),
                     trackstersmerged = cms.InputTag('hltTiclCandidate'),
                     trackstersInCand = cms.InputTag('hltTiclCandidate'))

    process.TFileService = cms.Service("TFileService",
                                       fileName=cms.string("histoHLT.root")
                                       )
    process.FEVTDEBUGHLToutput_step = cms.EndPath(
        process.FEVTDEBUGHLToutput+process.recHitMapProducerHLT+process.lcAssocByEnergyScoreProducerHLT+process.layerClusterCaloParticleAssociationProducerHLT+process.scAssocByEnergyScoreProducerHLT+process.layerClusterSimClusterAssociationProducerHLT+process.tpClusterProducerHLT+process.quickTrackAssociatorByHitsHLT+process.trackingParticleRecoTrackAsssociationHLT+process.filteredLayerClustersSimTrackstersHLT+process.ticlSimTrackstersHLT+process.tracksterSimTracksterAssociationLinkingHLT+process.tracksterSimTracksterAssociationPRHLT+process.tracksterSimTracksterAssociationLinkingbyCLUE3DHLT+process.tracksterSimTracksterAssociationPRbyCLUE3DHLT+process.tracksterSimTracksterAssociationLinkingPUHLT+process.tracksterSimTracksterAssociationPRPUHLT+process.ticlDumperHLT)
    return process
