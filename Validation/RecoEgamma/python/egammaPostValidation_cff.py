import FWCore.ParameterSet.Config as cms

from Validation.RecoEgamma.photonPostProcessor_cff import *
from Validation.RecoEgamma.electronPostValidationSequence_cff import *
from Validation.RecoEgamma.egammaPostProcessor_cfi import *

egammaPostValidation = cms.Sequence(photonPostProcessor+electronPostValidationSequence+postProcessorEGamma)
