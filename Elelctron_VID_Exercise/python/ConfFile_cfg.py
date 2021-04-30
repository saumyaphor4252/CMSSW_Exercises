import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

inputFilesMiniAOD = cms.untracked.vstring('/store/data/Run2017B/SingleElectron/MINIAOD/31Mar2018-v1/90000/FC89D712-AF37-E811-AD13-008CFAC93F84.root')

# Set up input/output depending on the format
# You can list here either AOD or miniAOD files, but not both types mixed
#
useAOD = False
if useAOD == True :
    inputFiles = inputFilesAOD
    outputFile = "electron_ntuple.root"
    print("AOD input files are used")
else :
    inputFiles = inputFilesMiniAOD
    outputFile = "electron_ntuple_mini.root"
    print("MiniAOD input files are used")

process.source = cms.Source("PoolSource", fileNames = inputFiles )

#
# Set up electron ID (VID framework)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#
# Configure the module
#
process.demo = cms.EDAnalyzer('DemoAnalyzer',
 electrons = cms.InputTag("slimmedElectrons"),

 # ID decisions (common to all formats)
 #
 # all IDs listed below are available given the content of "my_id_modules" defined above.
 # only one is exercised for this example.
 #

 #eleIdMapLoose = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose"),
 #eleIdMapMedium = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium"),
 #eleIdMapTight = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight"),
 eleIdMapLoose = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-loose"),
 eleIdMapMedium = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-medium"),
 eleIdMapTight = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-tight")
)



process.TFileService = cms.Service("TFileService",
fileName = cms.string(outputFile)
)

process.p = cms.Path(process.egmGsfElectronIDSequence*process.demo)
