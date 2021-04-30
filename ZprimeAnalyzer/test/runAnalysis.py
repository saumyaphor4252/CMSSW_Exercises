import FWCore.ParameterSet.Config as cms

process = cms.Process("ZprimeAnalyser")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'root://cms-xrd-global.cern.ch//store/data/Run2017B/SingleElectron/MINIAOD/31Mar2018-v1/90000/FC89D712-AF37-E811-AD13-008CFAC93F84.root'
    ) 
)

process.tree = cms.EDAnalyzer("ZprimeAnalyser",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    elecs = cms.InputTag("slimmedElectrons"),
    jets = cms.InputTag("slimmedJets"),
    mets = cms.InputTag("slimmedMETs"),
    mcLabel = cms.InputTag("prunedGenParticles"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("out_tree.root"
))


process.p = cms.Path(process.tree)

