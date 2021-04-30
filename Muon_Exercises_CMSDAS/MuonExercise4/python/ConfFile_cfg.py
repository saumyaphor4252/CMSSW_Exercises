import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 2000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/s/ssaumya/DataFiles/dymm.root'))

process.demo = cms.EDAnalyzer('MuonExercise4',
#vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#muons=cms.InputTag("slimmedMuons"),
#pruned=cms.InputTag("prunedGenParticles"),
#packed=cms.InputTag("packedGenParticles"),
#bits=cms.InputTag("TriggerResults"),
#objects=cms.InputTag("selectedPatTrigger"),
#prescales=cms.InputTag("patTrigger")
# ORDER NOT WORKING
     vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
     muons = cms.InputTag("slimmedMuons"),
     bits = cms.InputTag("TriggerResults","","HLT"),
     prescales = cms.InputTag("patTrigger"),
     objects = cms.InputTag("selectedPatTrigger"),
     packed = cms.InputTag("packedGenParticles"),
     pruned = cms.InputTag("prunedGenParticles")

)

process.TFileService = cms.Service("TFileService",fileName = cms.string('histos4.root'))

process.p = cms.Path(process.demo)
