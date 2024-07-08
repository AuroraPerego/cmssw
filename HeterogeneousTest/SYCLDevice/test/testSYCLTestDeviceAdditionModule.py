import FWCore.ParameterSet.Config as cms

process = cms.Process('TestSYCLTestDeviceAdditionModule')
process.load('HeterogeneousCore.SYCLCore.ProcessAcceleratorSYCL_cfi')

process.source = cms.Source('EmptySource')

process.syclTestDeviceAdditionModule = cms.EDAnalyzer('SYCLTestDeviceAdditionModule',
    size = cms.uint32( 1024*1024 )
)

process.path = cms.Path(process.syclTestDeviceAdditionModule)

process.maxEvents.input = 1
