import FWCore.ParameterSet.Config as cms

process = cms.Process("MUONEX4TAGPROBE")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
                                'root://cmseos.fnal.gov///store/user/cmsdas/2018/short_exercises/Muons/Samples/SingleMuon_Run2016G_18Apr2017_10.root',
                                'root://cmseos.fnal.gov///store/user/cmsdas/2018/short_exercises/Muons/Samples/SingleMuon_Run2016G_18Apr2017_3.root',
                                'root://cmseos.fnal.gov///store/user/cmsdas/2018/short_exercises/Muons/Samples/SingleMuon_Run2016G_18Apr2017_4.root',
                                'root://cmseos.fnal.gov///store/user/cmsdas/2018/short_exercises/Muons/Samples/SingleMuon_Run2016G_18Apr2017_9.root',
                            )
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016LegacyRepro_v4', '')

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
                                   fileName = cms.string("muonExercise4TagProbe_data_output.root"),
                                   closeFileFast = cms.untracked.bool(False)
)
