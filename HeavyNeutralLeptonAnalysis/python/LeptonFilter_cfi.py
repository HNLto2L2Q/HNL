import FWCore.ParameterSet.Config as cms

LeptonsFilter = cms.EDFilter('LeptonFilter',
                             triggerResultSrc      = cms.InputTag("TriggerResults","","HLT"),
                             muonSrc               = cms.InputTag("slimmedMuons"),
                             electronSrc           = cms.InputTag("slimmedElectrons"),
                             MinimalNumberOfMuons  = cms.uint32(2),
                             MinimalNumberOfElectrons = cms.uint32(2)
)
