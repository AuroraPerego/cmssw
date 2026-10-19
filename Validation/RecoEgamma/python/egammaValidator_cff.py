import FWCore.ParameterSet.Config as cms

from DQMServices.Core.DQM_cfg import *
from Validation.RecoEgamma.egammaValidator_cfi import egammaValidator as _egammaValidator

egammaValidator = _egammaValidator.clone()
