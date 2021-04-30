import FWCore.ParameterSet.Config as cms

process = cms.Process("Muontry")

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 2000
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)

process.demo = cms.EDAnalyzer('Muontry',
                          MuonTag    = cms.InputTag("slimmedMuons"),
        VertexTag  = cms.InputTag("offlineSlimmedPrimaryVertices"),
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("muontry_output.root"),
                                   closeFileFast = cms.untracked.bool(False)
)

process.p = cms.Path(process.demo)
