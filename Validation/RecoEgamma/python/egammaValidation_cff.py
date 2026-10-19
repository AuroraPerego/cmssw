import FWCore.ParameterSet.Config as cms

from Validation.RecoEgamma.electronValidationSequence_cff import *
from Validation.RecoEgamma.photonValidationSequence_cff import *
from Validation.RecoEgamma.egammaValidator_cff import egammaValidator

photonValidation.isRunCentrally = True
pfPhotonValidation.isRunCentrally = True
#oldpfPhotonValidation.isRunCentrally = True
tkConversionValidation.isRunCentrally = True

egammaValidatorSequence = cms.Sequence(egammaValidator)

egammaValidation = cms.Sequence(electronValidationSequence+photonValidationSequence)
