from os import *
from os import path as path
import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
from RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi import *
from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
from RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi import *
from HNL.HeavyNeutralLeptonAnalysis.ele_Sequence_cff import addElectronSequence
from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import *

import sys

import re
import importlib
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
#options.register('Flag', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Apply trigger matching for signal objects. Default: True")
options.register('isMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Running on MC samples")
options.parseArguments()

#hasLHE_ = options.Flag

hasLHE_ = True

debugLevel    = -1 

isMC_          = options.isMC
isMCSignal_    = False
#hasLHE_       = False #Only for MC with Matrix Element generators

algorithm     = "AK4PFchs"

GT_MC = '94X_mc2017_realistic_v10'#94X_mc2017_realistic_v14
GT_DATA = '92X_dataRun2_2017Repro_v4'#94X_dataRun2_v6

GT      =  GT_MC if isMC_ else GT_DATA

#system('ls -ltr')

process = cms.Process("AnalysisProc")
process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('HNL.HeavyNeutralLeptonAnalysis.LeptonFilter_cfi')

#process.LeptonsFilter.MinimalNumberOfMuons = cms.untracked.int32(2)
#process.LeptonsFilter.MinimalNumberOfElectrons = cms.untracked.int32(2)

#b-tagging
process.load('RecoBTag/Configuration/RecoBTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1

import FWCore.PythonUtilities.LumiList as LumiList
LumiList.LumiList().getVLuminosityBlockRange()

#from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, GT)

#input_file_list = options.inputFile
#myfilelist      = cms.untracked.vstring()
##if input_file_list != None:
#list_from_input_list = open("list.txt", "r")
#lines = list_from_input_list.readlines()
#stripped_lines = map(lambda x: x.rstrip("\n"), lines)
#for line in stripped_lines:
#    myfilelist.extend([line])
    
#print "my file " , myfilelist
    
#-------------------------------------------------------data section

#data 2017 lumi 41.86 /fb

#if(isMC_ == False):
#    if len('/afs/cern.ch/user/a/atalierc/CMSSW_9_4_13_patch4/src/HNL/data/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt') > 0:
#        import FWCore.PythonUtilities.LumiList as LumiList
#        process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt').getVLuminosityBlockRange()
#        import FWCore.PythonUtilities.LumiList as LumiList
#        import FWCore.ParameterSet.Types as CfgTypes
#        process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#        JSONfile = '/afs/cern.ch/user/a/atalierc/CMSSW_9_4_13_patch4/src/HNL/data/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
#        myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
#        process.source.lumisToProcess.extend(myLumis)

        
#-------------------------------------------------------------------
    
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    )
    
    

process.source = cms.Source("PoolSource", 
                            fileNames =  cms.untracked.vstring(

#'root://xrootd-cms.infn.it//store/data/Run2017B/SingleElectron/MINIAOD/31Mar2018-v1/30000/04B05308-0038-E811-99AB-008CFAC94314.root')
'root://xrootd-cms.infn.it//store/data/Run2017B/SingleMuon/MINIAOD/PromptReco-v2/000/299/329/00000/D6E915C7-3E6D-E711-8384-02163E014126.root')
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/100000/D4B7750A-4D94-E811-B78B-842B2B1\
#81788.root') 
#'file:heavyNeutrino_1.root')
                           #fileNames = myfilelist
#'file:/afs/cern.ch/user/a/atalierc/Merged_3GeVgood.root'
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/100000/D4B7750A-4D94-E811-B78B-842B2B181788.root'
#'file:/afs/cern.ch/user/a/atalierc/public/HIG-RunIIFall17MiniAODv2-00666_99.root'
#'file:/afs/cern.ch/user/a/atalierc/CMSSW_9_4_13_patch4/src/HNL/heavyNeutrino_150.root'
#'file:/afs/cern.ch/user/a/atalierc/CMSSW_9_4_13/src/HNL/HeavyNeutralLeptonAnalysis/test/Signal-RunIIFall17MiniAODv2-00666.root'
#'file:/afs/cern.ch/user/a/atalierc/Signal-RunIIFall17MiniAODv2-00666_38.root'
#'file:/afs/cern.ch/user/a/atalierc/CMSSW_9_4_13/src/Signal_300GeV.root'
#'file:/afs/cern.ch/user/a/atalierc/HIG-RunIIFall17MiniAODv2-00666_99.root' 
#'file:/afs/cern.ch/user/a/atalierc/CMSSW_9_4_13/src/HNL/HeavyNeutralLeptonAnalysis/test/04C8B197-4042-E811-BD46-FA163E81B685.root'
#'/store/mc/RunIISpring16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v4/40000/04E3024A-EF2B-E611-9794-02163E013F44.root'#2016 sample
#'file:/afs/cern.ch/user/a/atalierc/CMSSW_9_4_10/src/HNL/HeavyNeutralLeptonAnalysis/test/HIG-RunIIFall17MiniAODv2-00666.root'
#####################################################'file:/afs/cern.ch/user/a/atalierc/public/HIG-RunIIFall17MiniAODv2-00666_99.root' 
#'file:/afs/cern.ch/user/a/atalierc/HIG-RunIIFall17MiniAODv2-00666_324.root'
#'file:/afs/cern.ch/user/a/atalierc/HIG-RunIIFall17MiniAODv2-00666_325.root'
#'file:/afs/cern.ch/user/a/atalierc/HIG-RunIIFall17MiniAODv2-00666_327.root'
#'file:/afs/cern.ch/user/a/atalierc/HIG-RunIIFall17MiniAODv2-00666_329.root'
#'file:/afs/cern.ch/user/a/atalierc/HIG-RunIIFall17MiniAODv2-00666_445.root'
#'file:/afs/cern.ch/user/a/atalierc/HIG-RunIIFall17MiniAODv2-00666_558.root'
#'file:/afs/cern.ch/user/a/atalierc/HIG-RunIIFall17MiniAODv2-00666_674.root'
#'file:/afs/cern.ch/user/a/atalierc/HIG-RunIIFall17MiniAODv2-00666_676.root'
#'file:/afs/cern.ch/user/a/atalierc/HIG-RunIIFall17MiniAODv2-00666_678.root'
#'file:/afs/cern.ch/user/a/atalierc/HIG-RunIIFall17MiniAODv2-00666_768.root'
#'file:/afs/cern.ch/user/a/atalierc/HIG-RunIIFall17MiniAODv2-00666_769.root'
#'file:/afs/cern.ch/user/a/atalierc/CMSSW_9_4_13/src/HNL/HeavyNeutralLeptonAnalysis/test/04C8B197-4042-E811-BD46-FA163E81B685.root'
#'/store/mc/RunIISpring16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v4/40000/04E3024A-EF2B-E611-9794-02163E013F44.root'#2016 sample
#'file:/afs/cern.ch/user/a/atalierc/CMSSW_9_4_10/src/HNL/HeavyNeutralLeptonAnalysis/test/HIG-RunIIFall17MiniAODv2-00666.root'
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/50000/FE8D896F-386C-E811-AAAB-001E6779264E.root') 
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/20000/3279EE6B-108C-E811-804C-F01FAFD8EA6A.root' 
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/100000/D8FD945E-5588-E811-A866-D8D385FF33B9.root'
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/120000/96922A9A-B5B8-E811-986B-02163E017F81.root'
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_BGenFilter_Wpt-200toInf_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/30000/FEAFC6E4-ED82-E811-8398-0025904CF766.root'
#store/mc/RunIIFall17MiniAOD/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/70000/FEFA6784-D0F6-E711-A31A-008CFAC93ECC.root'
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root'
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root',
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Fall17/displaced/HeavyNeutrino_lljj_M-5_V-0.00836660026534_mu_massiveAndCKM_LO/heavyNeutrino_1.root'
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-50_Zpt-150toInf_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/RECOSIMstep_94X_mc2017_realistic_v10-v1/50000/EE9CC3E0-0DED-E711-BCAC-00E081CB560C.root'
#root://xrootd-cms.infn.it//WZTo3LNu_0Jets_MLL-4to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
#'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/80000/4E597432-24BE-E611-ACBB-00266CFFBFC0.root'
#cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/80000/4E597432-24BE-E611-ACBB-00266CFFBFC0.root'
#'root://cms-xrd-global.cern.ch//store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/1110000/72446D9C-D89C-E611-9060-002590A3C984.root'
#'file:/pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018/displaced/HeavyNeutrino_lljj_M-2_V-0.00316227766017_mu_massiveAndCKM_LO/heavyNeutrino_40.root'
#'file:/pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018/displaced/HeavyNeutrino_lljj_M-8_V-0.004472135955_mu_massiveAndCKM_LO/heavyNeutrino_96.root'
#'file:/pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-1_V-0.00836660026534_e_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_96.root'
)#)

process.TFileService = cms.Service("TFileService", fileName = cms.string("Analysis_output_prova.root"))
process.load('HNL.HeavyNeutralLeptonAnalysis.MetaNtuplizer_cfi')
process.metaTree.isMC = isMC_
process.metaTree.weightsSrc = cms.InputTag('externalLHEProducer')
process.metaTree.globalTag = GT
process.metaTree.args = cms.string('USELESS') #FILL ME!
process.metaTree.hasLHE = cms.bool(hasLHE_ and isMC_)

process.load('HNL.DisplacedAdaptiveVertexFinder.displacedInclusiveVertexing_cff')

addElectronSequence(process)

process.load("CondCore.CondDB.CondDB_cfi")

process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")

process.jetCorrFactors = process.updatedPatJetCorrFactors.clone(
    src = cms.InputTag("slimmedJets"),
    levels = ['L1FastJet', 
              'L2Relative', 
              'L3Absolute',
              'L2L3Residual'],
    payload = 'AK4PFchs') 

process.slimmedJetsJEC = process.updatedPatJets.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jetCorrFactors"))
    )

runMetCorAndUncFromMiniAOD(process, isData = not(isMC_), jetCollUnskimmed = "slimmedJetsJEC", postfix="NewJEC")

from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat, switchOnVIDElectronIdProducer, setupAllVIDIdsInModule, setupVIDElectronSelection
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

process.egmGsfElectronIDs.physicsObjectSrc = "slimmedElectrons"

id_modules = [
    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff',
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff'
    ]

for mod in id_modules:
    setupAllVIDIdsInModule(process, mod, setupVIDElectronSelection)

    setupEgammaPostRecoSeq(process,
                           runVID=True, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                           era='2017-Nov17ReReco')

#process.Filter = cms.EDFilter('LeptonFilter',
#                              MinimalNumberOfMuons = cms.untracked.int32(2),
#                              MinimalNumberOfElectrons = cms.untracked.int32(2)
#    )
    
#MET correction for systematics
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

runMetCorAndUncFromMiniAOD (
    process,
    isData = True, # false for MC
    fixEE2017 = True,
    fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
    postfix = "ModifiedMET"
    )

if isMC_:
    process.jetSmearing = cms.EDProducer('SmearedPATJetProducer',
                                         src          = cms.InputTag('slimmedJetsJEC'),
                                         enabled      = cms.bool(True),
                                         rho          = cms.InputTag("fixedGridRhoFastjetAll"),
                                         algo         = cms.string('AK4PFchs'),
                                         algopt       = cms.string('AK4PFchs_pt'),
                                         genJets      = cms.InputTag('slimmedGenJets'),
                                         dRMax        = cms.double(0.2),
                                         dPtMaxFactor = cms.double(3),
                                         debug        = cms.untracked.bool(False),
                                         variation    = cms.int32(0),
                                         )
    process.jetSmearingUp   = process.jetSmearing.clone(variation    = cms.int32(1))
    process.jetSmearingDown = process.jetSmearing.clone(variation    = cms.int32(-1))

process.prefiringweight = cms.EDProducer("L1ECALPrefiringWeightProducer",
                                         ThePhotons = cms.InputTag("slimmedPhotons"),
                                         TheJets = cms.InputTag("slimmedJets"),
                                         L1Maps = cms.string("${CMSSW_BASE}/src/L1Prefiring/EventWeightProducer/data/L1PrefiringMaps_new.root"), # update this line with the location of this file
                                         DataEra = cms.string("2017BtoF"), #Use 2016BtoH for 2016
                                         UseJetEMPt = cms.bool(False), #can be set to true to use jet prefiring maps parametrized vs pt(em) instead of pt
                                         PrefiringRateSystematicUncty = cms.double(0.2) #Minimum relative prefiring uncty per object
                                         )



process.HeavyNeutralLepton = cms.EDAnalyzer('HeavyNeutralLeptonAnalysis',#HeavyNeutralLeptonAnalysis
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
                                            pfMETSrc              = cms.InputTag("slimmedMETsModifiedMET"),
                                            triggerResultSrc      = cms.InputTag("TriggerResults","","HLT"),
                                            metFilterResultSrc    = cms.InputTag("TriggerResults","","PAT"),
                                            genParticleSrc        = cms.InputTag("prunedGenParticles"),
                                            genEventInfoProduct   = cms.InputTag("generator"),
                                            PUInfo                = cms.InputTag("slimmedAddPileupInfo"),
                                            lheEventProducts      = cms.InputTag("externalLHEProducer"),
                                            SecondaryVertices     = cms.InputTag("displacedInclusiveSecondaryVertices"), 
                                            bDiscbb               = cms.vstring("pfDeepCSVJetTags:probb"),
                                            bDiscbbb              = cms.vstring("pfDeepCSVJetTags:probbb"),
                                            bDiscbc               = cms.vstring("pfDeepCSVJetTags:probc"),
                                            #electronsMva          = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto"),#egmGsfElectronIDs:mvaEleID-Fall17-iso-V1-wp90 #egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto mvaEleID-Fall17-noIso-V1-wp90
                                            electronsVeto  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto"),
                                            electronsLoose = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose"),
                                            electronsMedium= cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium"),
                                            electronsTight = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight"),
                                            jetsSmeared     = cms.InputTag("jetSmearing" if isMC_ else "slimmedJetsJEC"),
                                            jetsSmearedUp   = cms.InputTag("jetSmearingUp" if isMC_ else "slimmedJetsJEC"),
                                            jetsSmearedDown = cms.InputTag("jetSmearingDown" if isMC_ else "slimmedJetsJEC")
                                            )

process.MessageLogger = cms.Service("MessageLogger",
                                suppressWarning= cms.untracked.vstring('displacedInclusiveVertexFinder')
)
if (isMC_):
    process.p = cms.Path(
        process.metaTree
        *process.LeptonsFilter
        *process.egmGsfElectronIDSequence
        *process.prefiringweight
        *process.fullPatMetSequenceModifiedMET
        *process.electronMVAValueMapProducer
        *process.btagging
        *process.displacedInclusiveVertexing
        *process.ele_Sequence
        *process.jetCorrFactors
        *process.slimmedJetsJEC
        *process.jetSmearing
        *process.jetSmearingUp
        *process.jetSmearingDown
        *process.HeavyNeutralLepton
        )
else:
    process.p = cms.Path(
        process.LeptonsFilter
#        *process.egmGsfElectronIDSequence
#        *process.prefiringweight
#        *process.fullPatMetSequenceModifiedMET
#        *process.electronMVAValueMapProducer
#        *process.btagging
#        *process.displacedInclusiveVertexing
#        *process.ele_Sequence
#        *process.jetCorrFactors
#        *process.slimmedJetsJEC
#        *process.HeavyNeutralLepton
        *process.egmGsfElectronIDSequence
        *process.jetCorrFactors
        *process.slimmedJetsJEC
        *process.prefiringweight
        *process.electronMVAValueMapProducer
        *process.btagging
        *process.displacedInclusiveVertexing
        *process.HeavyNeutralLepton 


       )
