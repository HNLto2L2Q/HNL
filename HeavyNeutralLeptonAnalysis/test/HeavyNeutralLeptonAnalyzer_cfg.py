from os import path as path
import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD#MET correction

debugLevel    = -1 

isMC_         = True
isMCSignal_   = False
hasLHE_       = True #Only for MC with Matrix Element generators

#jec_tag_DATA  = 'JetCorrectorParametersCollection_Summer16_23Sep2016AllV4_DATA_AK4PFchs'
#jec_tag_MC    = 'JetCorrectorParametersCollection_Summer16_23Sep2016V4_MC_AK4PFchs'
#jec_file_DATA = 'sqlite_file:Summer16_23Sep2016AllV4_DATA.db'
#jec_file_MC   = 'sqlite_file:Summer16_23Sep2016V4_MC.db'
jec_tag       = jec_tag_MC if isMC_ else jec_tag_DATA
jec_file      = jec_file_MC if isMC_ else jec_file_DATA
algorithm     = "AK4PFchs"

#GT_MC   = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
#GT_DATA = '80X_dataRun2_2016SeptRepro_v7'
GT_MC = '94X_mc2017_realistic_v10'#94X_mc2017_realistic_v14
GT_DATA = '92X_dataRun2_2017Repro_v4'#94X_dataRun2_v6

GT      =  GT_MC if isMC_ else GT_DATA


process = cms.Process("AnalysisProc")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

#b-tagging
process.load('RecoBTag/Configuration/RecoBTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1

import FWCore.PythonUtilities.LumiList as LumiList
LumiList.LumiList().getVLuminosityBlockRange()

#from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, GT)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
)
process.source = cms.Source("PoolSource", 
                            fileNames =  cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root'
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root'
#'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/RECOSIMstep_94X_mc2017_realistic_v10-v1/50000/EE9CC3E0-0DED-E711-BCAC-00E081CB560C.root'
'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/80000/4E597432-24BE-E611-ACBB-00266CFFBFC0.root'
#cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/80000/4E597432-24BE-E611-ACBB-00266CFFBFC0.root'
#'root://cms-xrd-global.cern.ch//store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/1110000/72446D9C-D89C-E611-9060-002590A3C984.root'
#'file:/pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018/displaced/HeavyNeutrino_lljj_M-2_V-0.00316227766017_mu_massiveAndCKM_LO/heavyNeutrino_40.root'
#'file:/pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018/displaced/HeavyNeutrino_lljj_M-8_V-0.004472135955_mu_massiveAndCKM_LO/heavyNeutrino_96.root'
#'file:/pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-1_V-0.00836660026534_e_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_96.root'
))
process.TFileService = cms.Service("TFileService", fileName = cms.string("Analysis_output.root"))
process.load('HNL.HeavyNeutralLeptonAnalysis.MetaNtuplizer_cfi')
process.metaTree.isMC = isMC_
process.metaTree.weightsSrc = cms.InputTag('externalLHEProducer')
process.metaTree.globalTag = GT
process.metaTree.args = cms.string('USELESS') #FILL ME!
process.metaTree.hasLHE = cms.bool(hasLHE_ and isMC_)

process.load('HNL.DisplacedAdaptiveVertexFinder.displacedInclusiveVertexing_cff')

from HNL.HeavyNeutralLeptonAnalysis.ele_Sequence_cff import addElectronSequence

addElectronSequence(process)

process.load("CondCore.CondDB.CondDB_cfi")

process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")

process.jetCorrFactors = process.updatedPatJetCorrFactors.clone(
    src = cms.InputTag("slimmedJets"),
    levels = ['L1FastJet', 
              'L2Relative', 
              'L3Absolute',
              'L2L3Residual'],
    payload = 'AK4PFchs' ) 

process.slimmedJetsJEC = process.updatedPatJets.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jetCorrFactors"))
    labelName = "NewJEC"
    )

runMetCorAndUncFromMiniAOD(process, isData = !isMC_is, jetCollUnskimmed = NewJEC, postfix="NewJEC")

#process.ak4PFJetsCHSL1FastL2L3 = process.ak4PFCHSJetsL1.clone(correctors = ['ak4PFL1FastL2L3Corrector'])
#process.correctJets           = cms.Sequence(process.ak4PFL1FastL2L3CorrectorChain * process.ak4PFJetsCHSL1FastL2L3)

process.HeavyNeutralLepton = cms.EDAnalyzer('HeavyNeutralLeptonAnalysis',
                                            debugLevel            = cms.int32(debugLevel),
                                            isMC                  = cms.bool(isMC_),
                                            isMCSignal            = cms.bool(isMCSignal_),
                                            vtxSrc                = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                            rho                   = cms.InputTag("fixedGridRhoFastjetAll"),#cambiato tag?
                                            muonSrc               = cms.InputTag("slimmedMuons"),
                                            electronSrc           = cms.InputTag("slimmedElectrons"),
                                            recHitCollectionEBSrc = cms.InputTag("reducedEgamma","reducedEBRecHits"),
                                            recHitCollectionEESrc = cms.InputTag("reducedEgamma","reducedEERecHits"),
                                            tauSrc                = cms.InputTag("slimmedTaus"),
                                            packCandSrc           = cms.InputTag("packedPFCandidates"),
                                            jetSrc                = cms.InputTag("slimmedJetsJEC"),
                                            pfMETSrc              = cms.InputTag("slimmedMETs"),
                                            triggerResultSrc      = cms.InputTag("TriggerResults","","HLT"),
                                            metFilterResultSrc    = cms.InputTag("TriggerResults","","PAT"),
                                            genParticleSrc        = cms.InputTag("prunedGenParticles"),
                                            genEventInfoProduct   = cms.InputTag("generator"),
                                            PUInfo                = cms.InputTag("slimmedAddPileupInfo"),
                                            lheEventProducts      = cms.InputTag("externalLHEProducer"),
                                            SecondaryVertices     = cms.InputTag("displacedInclusiveSecondaryVertices"), 
                                            bDiscbb = cms.vstring("pfDeepCSVJetTags:probb"),
                                            bDiscbbb = cms.vstring("pfDeepCSVJetTags:probbb"),
                                            bDiscbc = cms.vstring("pfDeepCSVJetTags:probc"),
                                            electronsMva          = cms.InputTag("mvaEleID-Fall17-noIso-V1-wp90"),
                                            electronsVeto  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto"),
                                            electronsLoose = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose"),
                                            electronsMedium= cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium"),
                                            electronsTight = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight")
                                            )

#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string('test6.root')
#)

process.p = cms.Path(
    process.btagging
    *process.metaTree
    *process.displacedInclusiveVertexing 
    *process.ele_Sequence
    *process.jetCorrFactors
    *process.slimmedJetsJEC
    *process.HeavyNeutralLepton
    )

#process.e = cms.EndPath(process.out)
