import FWCore.ParameterSet.Config as cms

process = cms.Process("MuonExercise3b")

#sample='DY'
#sample='TT'
sample='QCD'

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring()
)

if sample == 'DY':
    process.source.fileNames.append('root://cmseos.fnal.gov///store/user/cmsdas/2019/short_exercises/Muons/Samples/dymm.root')
elif sample == 'TT':
    process.source.fileNames.append('root://cmseos.fnal.gov///store/user/cmsdas/2019/short_exercises/Muons/Samples/TTJets_DiLept_MadGraph_LO.root')
#    process.source.fileNames.append('root://cmseos.fnal.gov///store/user/cmsdas/2019/short_exercises/Muons/Samples/dymm.root')
elif sample == 'QCD':
#   process.source.fileNames.append('root://cmseos.fnal.gov///store/user/cmsdas/2019/short_exercises/Muons/Samples/QCD_Pt-20toInf_MuEnrichedPt15_Pythia.root')
     process.source.fileNames.append('root://cmseos.fnal.gov///store/user/cmsdas/2016/SHORT_EXERCISES/Muons/qcd.root')
else:
    print " *** ERROR: Please select a sample (DY, TT, QCD)! ***"
    exit(1)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 2000
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)

process.demo = cms.EDAnalyzer('MuonExercise3b',
                          MuonTag    = cms.InputTag("slimmedMuons"),
                          GenPartTag = cms.InputTag("prunedGenParticles"),
                          VertexTag  = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          IsQCD      = cms.untracked.bool(True)
)

#f sample == 'QCD':
#   process.demo = True

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("muonExercise3_output_"+sample+"2016.root"),
                                   closeFileFast = cms.untracked.bool(False)
)

process.p = cms.Path(process.demo)
