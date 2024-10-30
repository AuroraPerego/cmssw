import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import *
from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import *

from RecoLocalCalo.HGCalRecProducers.recHitMapProducer_cfi import recHitMapProducer

# patch particle flow clusters for HGC into local reco sequence
# (for now until global reco is going with some sort of clustering)
from RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGC_cfi import *
from RecoParticleFlow.PFClusterProducer.particleFlowClusterHGC_cfi import *
from RecoLocalCalo.HGCalRecProducers.hgcalMultiClusters_cfi import *
from RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cff import hgcalLayerClustersHFNose, hgcalLayerClustersEE, hgcalLayerClustersHSi, hgcalLayerClustersHSci, hgcalMergeLayerClusters

from RecoLocalCalo.HGCalRecProducers.hgcalHeterogeneousModules_cfi import hgcalSoARecHits, hgcalSoARecHitsLayerClusters, hgcalSoALayerClusters, hgcalLayerClustersFromSoA

hgcalLocalRecoTask1 = cms.Task( HGCalUncalibRecHit,
                                       HGCalRecHit,
                                       recHitMapProducer)

hgcalLocalRecoTask2 = cms.Task( hgcalLayerClustersHSi,
                                       hgcalLayerClustersHSci,
                                       hgcalMergeLayerClusters,
                                       hgcalMultiClusters,
                                       particleFlowRecHitHGC,
                                       particleFlowClusterHGCal )

hgcalLocalRecoTask = cms.Task(hgcalLocalRecoTask1,
                              hgcalLayerClustersEE,
                              hgcalLocalRecoTask2)

_heterogeneous_hgcalLocalRecoTask = cms.Task(hgcalLocalRecoTask1,
                                         hgcalSoARecHits,
                                         hgcalSoARecHitsLayerClusters,
                                         hgcalSoALayerClusters,
                                         hgcalLayerClustersFromSoA,
                                         hgcalLocalRecoTask2)

_hfnose_hgcalLocalRecoTask = hgcalLocalRecoTask.copy()
_hfnose_hgcalLocalRecoTask.add(hgcalLayerClustersHFNose)

from Configuration.ProcessModifiers.alpaka_cff import alpaka
alpaka.toReplaceWith(hgcalLocalRecoTask, _heterogeneous_hgcalLocalRecoTask)
alpaka.toModify(hgcalMergeLayerClusters,
         layerClustersEE = cms.InputTag("hgcalLayerClustersFromSoA",),
         time_layerclustersEE = cms.InputTag("hgcalLayerClustersFromSoA", "timeLayerCluster"))

from Configuration.Eras.Modifier_phase2_hfnose_cff import phase2_hfnose
phase2_hfnose.toReplaceWith(
    hgcalLocalRecoTask, _hfnose_hgcalLocalRecoTask )

hgcalLocalRecoSequence = cms.Sequence(hgcalLocalRecoTask)
