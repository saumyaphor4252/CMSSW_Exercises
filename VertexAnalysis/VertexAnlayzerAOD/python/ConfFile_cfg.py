import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:/afs/cern.ch/cms/Tutorials/PAT_tutorial_Summer12/AOD_DoubleMu_Run195013_numEvent1000.root'
    )
)

process.demo = cms.EDAnalyzer('VertexAnalyzerAOD',
vertices = cms.InputTag("offlinePrimaryVertices"),
tracks = cms.InputTag("generalTracks"))

process.TFileService = cms.Service("TFileService",
fileName = cms.string('OutputfileVertex_c22.root')
)


process.p = cms.Path(process.demo)
