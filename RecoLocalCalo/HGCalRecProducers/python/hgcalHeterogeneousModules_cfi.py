import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.HGCalRecProducers.hgCalSoARecHitsProducer_cfi import hgCalSoARecHitsProducer
from RecoLocalCalo.HGCalRecProducers.hgCalSoARecHitsLayerClustersProducer_cfi import hgCalSoARecHitsLayerClustersProducer
from RecoLocalCalo.HGCalRecProducers.hgCalSoALayerClustersProducer_cfi import hgCalSoALayerClustersProducer
from RecoLocalCalo.HGCalRecProducers.hgCalLayerClustersFromSoAProducer_cfi import hgCalLayerClustersFromSoAProducer

from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import HGCalUncalibRecHit
from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import HGCalRecHit
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import HGCAL_noises

hgcalSoARecHits = hgCalSoARecHitsProducer.clone(detector = cms.string('EE'),
                                                fcPerMip = HGCalUncalibRecHit.HGCEEConfig.fCPerMIP.value() + HGCalUncalibRecHit.HGCHEFConfig.fCPerMIP.value(),
                                                thicknessCorrection = HGCalRecHit.thicknessCorrection.value(),
                                                noises = HGCAL_noises.values.value() + HGCAL_noises.values.value(),
                                                dEdXweights = HGCalRecHit.layerWeights.value())
hgcalSoARecHitsLayerClusters = hgCalSoARecHitsLayerClustersProducer.clone(
                                                hgcalRecHitsSoA = cms.InputTag("hgcalSoARecHits"))
hgcalSoALayerClusters = hgCalSoALayerClustersProducer.clone(
                                                hgcalRecHitsLayerClustersSoA = cms.InputTag("hgcalSoARecHitsLayerClusters"),
                                                hgcalRecHitsSoA = cms.InputTag("hgcalSoARecHits"))
hgcalLayerClustersFromSoA = hgCalLayerClustersFromSoAProducer.clone(
                                                src = cms.InputTag("hgcalSoALayerClusters"),
                                                hgcalRecHitsLayerClustersSoA = cms.InputTag("hgcalSoARecHitsLayerClusters"),
                                                hgcalRecHitsSoA = cms.InputTag("hgcalSoARecHits"))
