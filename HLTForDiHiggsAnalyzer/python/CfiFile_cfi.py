import FWCore.ParameterSet.Config as cms

HLTForDiHiggsAnalyzer = cms.EDAnalyzer('HLTForDiHiggsAnalyzer',

    ### Trigger results
    hltresults  = cms.InputTag('TriggerResults::HLTX'),
    trigSummary = cms.InputTag('hltTriggerSummaryAOD','','HLTX'),
    #l1results       = cms.InputTag("hltGtStage2Digis", "", "HLTX"),                      

    ### Output file
    outRoot  = cms.untracked.string('efficiency.root')
)
