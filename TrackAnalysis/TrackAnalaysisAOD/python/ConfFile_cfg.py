import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root'
    )
)

process.demo = cms.EDAnalyzer('TrackAnalyzer',
tracks = cms.untracked.InputTag("generalTracks")
)

process.TFileService = cms.Service("TFileService",
fileName = cms.string('outputfile_Track.root')
)


process.p = cms.Path(process.demo)
