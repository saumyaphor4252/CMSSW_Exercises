import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
     'root://cms-xrd-global.cern.ch//store/data/Run2017B/SingleElectron/MINIAOD/31Mar2018-v1/90000/FC89D712-AF37-E811-AD13-008CFAC93F84.root'
    )
)

process.demo = cms.EDAnalyzer('patTrackAnalyzer',
pfs = cms.InputTag("packedPFCandidates")
)

process.TFileService = cms.Service("TFileService",
fileName = cms.string('outputfile5.root')
)


process.p = cms.Path(process.demo)
