import FWCore.ParameterSet.Config as cms

displacedSVAssociationIVF = cms.EDProducer("DisplacedSVAssociator", 
                                           #pfJets          = cms.untracked.InputTag("ak4PFJetsCHS","",""), 
                                           #pfJets          = cms.untracked.InputTag("ak8PFJetsCHS","",""),
                                           muonTag                       = cms.untracked.InputTag("slimmedMuons","",""),
                                           secondaryVertices = cms.untracked.InputTag("displacedInclusiveSecondaryVertices"),
                                           primaryVertices   = cms.untracked.InputTag("offlinePrimaryVertices"),
                                           outputLabel       = cms.untracked.string('displacedIVFAssoc'),
                                           muonPtCut          = cms.untracked.double(20.0), #5
                                           algoName          = cms.untracked.string("oneOverDeltaR"),
                                           debug             = cms.untracked.int32(1)
                                           #jets = cms.untracked.InputTag("ak4PFJetsCHSL1FastL2L3","","")
)
