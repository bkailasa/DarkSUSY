#importing CMS-specific Python classes and functions
import FWCore.ParameterSet.Config as cms

process = cms.Process("fastJetNtupler")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.MessageLogger.cerr.threshold = 'INFO'

process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit = cms.untracked.int32(0))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False),SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(

'file:/lustrehome/mmaggi/DarkMuonJet/MNeu60_MHd30_Mad14.00_Mnd0.1/003/DarkSusy-RunIIAutumn18DR_step3-MiniAOD.root'
))

#process.demo = cms.EDAnalyzer('fastJetNtupler',jetpat = cms.untracked.InputTag("slimmedJets"))

process.demo1 = cms.EDAnalyzer('fastJetNtupler',
                                muonTag = cms.untracked.InputTag("slimmedMuons"),
                                metTag = cms.untracked.InputTag("slimmedMETsPuppi")
                                )

process.TFileService = cms.Service("TFileService",fileName = cms.string('MNeu60_MHd30_Mad14.00_Mnd0.1.root'))


process.p = cms.Path(process.demo1)
