import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 2000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:/afs/cern.ch/user/s/ssaumya/DataFiles/dymm.root' )
#     'root://cmseos.fnal.gov///store/user/cmsdas/2019/short_exercises/Muons/Samples/dymm.root')
)

process.demo = cms.EDAnalyzer('MuonExercise1')

process.TFileService = cms.Service("TFileService",fileName = cms.string('histos1_final.root'))

process.p = cms.Path(process.demo)
