import FWCore.ParameterSet.Config as cms


process = cms.Process("TRIGGERSELECTOR")

process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(
                                'file:/afs/cern.ch/work/b/bkailasa/CMSSW/slc7/CMSSW_10_2_16_UL2/miniAOD_for_testing/DarkSusy-RunIIAutumn18DR_step3-MiniAOD.root'
                            )
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.TriggerSelector = cms.EDFilter("TriggerSelector",
                                  TriggerList = cms.vstring("HLT_Mu9","HLT_Quad")
)

process.thepath = cms.Path(process.TriggerSelector)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("TriggerSelector.root"),
                                   closeFileFast = cms.untracked.bool(False)
)
