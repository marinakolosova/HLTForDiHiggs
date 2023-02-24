import FWCore.ParameterSet.VarParsing as VarParsing
import FWCore.ParameterSet.Config as cms
import sys

# ============================================
#             Setup the process
# ============================================
from hltMC import *
#from hltMC_08Feb2023 import *
#from hltMC_Feb03_RAWAOD import *
#from hltMC_Feb03 import *
#from hlt_Data import *

# ===========================================
#             Output Files                   
# ===========================================
process.TFileService = cms.Service ('TFileService',
                                    fileName = cms.string(
                                        'efficiency.root',
                                    )
                                )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.HLTForDiHiggs = cms.EDAnalyzer("HLTForDiHiggsAnalyzer",

    ### Trigger results
    hltresults  = cms.InputTag('TriggerResults::HLTX'),
    trigSummary = cms.InputTag('hltTriggerSummaryAOD','','HLTX'),
    l1results   = cms.InputTag("hltGtStage2Digis", "", "HLTX"),
                                       
    ### Reco Objects
    Rho          = cms.InputTag('fixedGridRhoFastjetAll'),                                       
    pileupInfo   = cms.InputTag('addPileupInfo'),                                       
    vertex       = cms.InputTag("offlinePrimaryVertices"),
    muons        = cms.InputTag('muons'),
    electrons    = cms.InputTag("gedGsfElectrons"),
    jets         = cms.InputTag("ak4PFJetsPuppi"),
    bjets        = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    jetCorrector = cms.InputTag("ak4PFCHSL1FastL2L3Corrector"),                                       
    fatjets      = cms.InputTag("ak8PFJetsPuppi"),
    genparticles = cms.InputTag("genParticles"),
    genjets      = cms.InputTag("ak4GenJets"),
    genfatjets   = cms.InputTag("ak8GenJets"),
    isMC         = cms.bool(True),
    isSignal     = cms.bool(True),

    ### Output file
    outRoot  = cms.untracked.string('efficiency.root')
)                           

process.p = cms.FinalPath(process.HLTForDiHiggs)
process.schedule.append(process.p)
# ===========================================
#            Process Summary
# ===========================================
#process.output = cms.OutputModule("PoolOutputModule",
#                                  outputCommands = cms.untracked.vstring(
#                                  "keep *",
#                                  ),
#                                  fileName = cms.untracked.string("CMSSW.root")
#                                  )
#process.ep = cms.EndPath(process.output)
