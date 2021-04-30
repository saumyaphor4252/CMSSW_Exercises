import FWCore.ParameterSet.Config as cms

process = cms.Process("MUONEX4TAGPROBE")

process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
                                'root://cmseos.fnal.gov///store/user/hats/2018/Muon/Samples/dymm.root',
                                'root://cmseos.fnal.gov///store/user/hats/2018/Muon/Samples/dymm_1.root'
                            )
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 2000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.muonexercise4tp = cms.EDAnalyzer("MuonExercise4TagProbe",
                                         TagTriggerList = cms.vstring(
                                             "HLT_IsoMu24",
                                             "HLT_IsoTkMu24"
                                         ),
                                         TriggerList = cms.vstring(
                                             "HLT_IsoMu24",
                                             "HLT_IsoTkMu24"
                                         )
)

process.thepath = cms.Path(process.muonexercise4tp)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("muonExercise4TagProbe_output.root"),
                                   closeFileFast = cms.untracked.bool(False)
)
