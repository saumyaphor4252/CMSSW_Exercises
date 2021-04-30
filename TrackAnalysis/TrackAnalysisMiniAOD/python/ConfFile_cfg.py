import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )
inputFilesMiniAOD = cms.untracked.vstring('/store/data/Run2017B/SingleElectron/MINIAOD/31Mar2018-v1/90000/FC89D712-AF37-E811-AD13-008CFAC93F84.root')
inputFiles = inputFilesMiniAOD
process.source = cms.Source("PoolSource", fileNames = inputFiles)

process.demo = cms.EDAnalyzer('TrackAnalysis',
tracks = cms.InputTag("packedPFCandidates")
)

process.TFileService = cms.Service("TFileService",
fileName = cms.string('outputfile6b.root')
)

process.p = cms.Path(process.demo)
