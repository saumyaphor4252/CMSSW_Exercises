import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 2000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(
#'file:/afs/cern.ch/user/s/ssaumya/DataFiles/dymm.root'))
'root://cmseos.fnal.gov///store/user/hats/2018/Muon/Samples/dymm.root'))
#'root://cmseos.fnal.gov///store/user/hats/2018/Muon/Samples/dymm_1.root'))
#'root://cmseos.fnal.gov///store/user/hats/2017/muon/dymm.root'))
#'root://cmseos.fnal.gov///store/user/cmsdas/2016/SHORT_EXERCISES/Muons/dymm.root'))

process.demo = cms.EDAnalyzer('MuonExercise2')

process.TFileService = cms.Service("TFileService",fileName = cms.string('histos2try22.root'))

process.p = cms.Path(process.demo)
