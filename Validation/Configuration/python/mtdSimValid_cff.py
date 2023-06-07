import FWCore.ParameterSet.Config as cms

# MTD validation sequences
from Validation.MtdValidation.btlSimHitsValid_cfi import btlSimHitsValid
from Validation.MtdValidation.btlDigiHitsValid_cfi import btlDigiHitsValid
from Validation.MtdValidation.btlLocalRecoValid_cfi import btlLocalRecoValid
from Validation.MtdValidation.etlLocalRecoValid_cfi import etlLocalRecoValid
from Validation.MtdValidation.etlSimHitsValid_cfi import etlSimHitsValid
from Validation.MtdValidation.etlDigiHitsValid_cfi import etlDigiHitsValid
from Validation.MtdValidation.mtdTracksValid_cfi import mtdTracksValid
from Validation.MtdValidation.vertices4DValid_cfi import vertices4DValid
from Validation.MtdValidation.mtdEleIsoAnalyzer_cfi import mtdEleIsoAnalyzer
from Validation.MtdValidation.mtdEleIsoNtupler_cfi import mtdEleIsoNtupler
#from SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi import trackingParticleRecoTrackAsssociation
#from SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi import quickTrackAssociatorByHits
#
#quickTrackAssociator_lowCut = quickTrackAssociatorByHits.clone(
#    Cut_RecoToSim = cms.double(0.01),
#)
#
#gsfAssociator = trackingParticleRecoTrackAsssociation.clone(
#    associator = ("quickTrackAssociator_lowCut"),
#)

mtdSimValid  = cms.Sequence(btlSimHitsValid  + etlSimHitsValid )
mtdDigiValid = cms.Sequence(btlDigiHitsValid + etlDigiHitsValid)
mtdRecoValid = cms.Sequence(btlLocalRecoValid  + etlLocalRecoValid + mtdTracksValid + vertices4DValid + mtdEleIsoAnalyzer + mtdEleIsoNtupler)

