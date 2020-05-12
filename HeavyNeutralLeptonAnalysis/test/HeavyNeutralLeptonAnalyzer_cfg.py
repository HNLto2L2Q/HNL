from os import path as path
import FWCore.ParameterSet.Config as cms

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
from RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi import *
from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
from RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi import *

#parser
#import argparse
#from input_crab import dataset_files
hasLHE_ = True

#parser = argparse.ArgumentParser()
#parser.add_argument("PythonOn", type=int, help="hasLHE congif ON when 1")
#args = parser.parse_args()
#pyEnable = args.PythonOn

#for dataset, infos in dataset_files.items():
#    if infos[2] == 0:
#        hasLHE_ = False
#    else:
#        hasLHE_ = True

#if pyEnable == 1:
#    hasLHE_ == True
#else:
#    hasLHE_ == False

#print hasLHE_

debugLevel    = -1

isMC_         = True
isMCSignal_   = False
#hasLHE_       = False #Only for MC with Matrix Element generators

if isMCSignal_ : print "this the signal please note this muon signals or electron signal "


algorithm     = "AK4PFchs"

GT_MC = '94X_mcRun2_asymptotic_v3'#94X_mc2017_realistic_v14
GT_DATA = '94X_dataRun2_v10'#94X_dataRun2_v6

GT      =  GT_MC if isMC_ else GT_DATA


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

process.MessageLogger.cerr.FwkReport.reportEvery = 500

import FWCore.PythonUtilities.LumiList as LumiList
LumiList.LumiList().getVLuminosityBlockRange()

#from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, GT)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
                            fileNames =  cms.untracked.vstring(
#'file:/afs/cern.ch/user/a/atalierc/CMSSW_9_4_10/src/HNL/HeavyNeutralLeptonAnalysis/test/HIG-RunIIFall17MiniAODv2-00666.root'

#'/store/mc/RunIISpring16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v4/40000/04E3024A-EF2B-E611-9794-02163E013F44.root'#2016 sample
#'file:/afs/cern.ch/user/a/atalierc/CMSSW_9_4_10/src/HNL/HeavyNeutralLeptonAnalysis/test/HIG-RunIIFall17MiniAODv2-00666.root'
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/50000/FE8D896F-386C-E811-AAAB-001E6779264E.root'
#'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv3/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/100000/4C8087A0-F5E9-E811-829E-001E6724807F.root'
#'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv3/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/280000/42B83DB5-9C21-E911-AE73-0CC47AC52A94.root'
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/20000/3279EE6B-108C-E811-804C-F01FAFD8EA6A.root'
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/100000/D8FD945E-5588-E811-A866-D8D385FF33B9.root'
#'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv3/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/100000/4C8087A0-F5E9-E811-829E-001E6724807F.root'
#'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv3/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/00000/06D08F8C-FBC2-E811-AC7B-0090FAA57460.root'
#'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv3/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/280000/009D1E1A-111F-E911-AD7E-6CC2173DABE0.root'
#'root://xrootd-cms.infn.it//store//mc/RunIISummer16MiniAODv3/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/280000/1A569771-871E-E911-AE14-0CC47AC17550.root'
'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv3/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/00000/2EBE1EFA-CFC2-E811-96B1-0090FAA57FA4.root'
#'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv3/WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/80000/08FECAA0-2CEF-E811-B07D-008CFAC93C10.root'
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/120000/96922A9A-B5B8-E811-986B-02163E017F81.root'
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAODv2/QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/70000/F0B31586-DE42-E811-9BF1-0242AC1C0500.root'
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
#'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/80000/4E597432-24B-E611-ACBB-00266CFFBFC0.root'
#cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/80000/4E597432-24BE-E611-ACBB-00266CFFBFC0.root'
#'root://cms-xrd-global.cern.ch//store/data/Run2016G/SingleMuon/MINIAOD/23Sep2016-v1/1110000/72446D9C-D89C-E611-9060-002590A3C984.root'
#'file:/pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018/displaced/HeavyNeutrino_lljj_M-2_V-0.00316227766017_mu_massiveAndCKM_LO/heavyNeutrino_40.root'
#'file:/pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018/displaced/HeavyNeutrino_lljj_M-8_V-0.004472135955_mu_massiveAndCKM_LO/heavyNeutrino_96.root'
#'file:/pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17/displaced/HeavyNeutrino_lljj_M-1_V-0.00836660026534_e_onshell_pre2017_leptonFirst_NLO/heavyNeutrino_96.root'
#'file:heavyNeutrino_10.root'
#'file:heavyNeutrino_150.root'
#'file:/afs/cern.ch/user/a/atalierc/public/Signal-RunIIFall17MiniAODv2-00666.root'
#'root://cms-xrd-global.cern.ch//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018/displaced/HeavyNeutrino_lljj_M-2_V-0.022360679775_mu_massiveAndCKM_LO/heavyNeutrino_55.root'
))
process.TFileService = cms.Service("TFileService", fileName = cms.string("Analysis_output.root"))
process.load('HNL.HeavyNeutralLeptonAnalysis.MetaNtuplizer_cfi')
process.metaTree.isMC = isMC_
process.metaTree.weightsSrc = cms.InputTag('externalLHEProducer')
process.metaTree.globalTag = GT
process.metaTree.args = cms.string('USELESS') #FILL ME!
process.metaTree.hasLHE = cms.bool(hasLHE_ and isMC_)

process.load('HNL.DisplacedAdaptiveVertexFinder.displacedInclusiveVertexing_cff')
process.load("CondCore.CondDB.CondDB_cfi")
process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")

process.MessageLogger = cms.Service("MessageLogger",
                                suppressWarning= cms.untracked.vstring('displacedInclusiveVertexFinder')
)

if isMC_: jetCorrectorLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
else : jetCorrectorLevels = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']


process.jetCorrFactors = process.updatedPatJetCorrFactors.clone(
    src = cms.InputTag("slimmedJets"),
    levels = jetCorrectorLevels,
    payload = 'AK4PFchs')

process.slimmedJetsJEC = process.updatedPatJets.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jetCorrFactors"))
    )


from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat, switchOnVIDElectronIdProducer, setupAllVIDIdsInModule, setupVIDElectronSelection
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

process.egmGsfElectronIDs.physicsObjectSrc = "slimmedElectrons"

id_modules = [
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff'
    ]

for mod in id_modules:
    setupAllVIDIdsInModule(process, mod, setupVIDElectronSelection)

    setupEgammaPostRecoSeq(process,
                           runVID=False, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                           era='2016-Legacy')

#process.options = cms.untracked.PSet(

#)
#process.options.numberOfThreads=cms.untracked.uint32(2)
#process.options.numberOfStreams=cms.untracked.uint32(0)

#process.Filter = cms.EDFilter('LeptonFilter',
#                              MinimalNumberOfMuons = cms.untracked.int32(2),
#                              MinimalNumberOfElectrons = cms.untracked.int32(2)
#    )

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


#prefiring correction for ECAL in the forward region
process.prefiringweight = cms.EDProducer("L1ECALPrefiringWeightProducer",
                                         ThePhotons = cms.InputTag("slimmedPhotons"),
                                         TheJets = cms.InputTag("slimmedJets"),
                                         L1Maps = cms.string("${CMSSW_BASE}/src/L1Prefiring/EventWeightProducer/data/L1PrefiringMaps_new.root"), # update this line with the location of this file
                                         DataEra = cms.string("2016BtoH"), #Use 2017BtoF for 2017
                                         UseJetEMPt = cms.bool(True), #can be set to true to use jet prefiring maps parametrized vs pt(em) instead of pt
                                         PrefiringRateSystematicUncty = cms.double(0.2) #Minimum relative prefiring uncty per object
                                         )

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
                                            electronsVeto  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
                                            electronsLoose = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
                                            electronsMedium= cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
                                            electronsTight = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
                                            jetsSmeared     = cms.InputTag("jetSmearing" if isMC_ else "slimmedJetsJEC"),
                                            jetsSmearedUp   = cms.InputTag("jetSmearingUp" if isMC_ else "slimmedJetsJEC"),
                                            jetsSmearedDown = cms.InputTag("jetSmearingDown" if isMC_ else "slimmedJetsJEC"),
                                            electronsEffectiveAreas = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt')

                                        )
triggers = [
    'HLT_IsoMu24_v*',
    'HLT_IsoTkMu24_v*',
    'HLT_IsoTkMu27_v*',
    'HLT_IsoMu27_v*',
    'HLT_Ele27_WPTight_Gsf_v*',

]
process.TriggerSelection = cms.EDFilter( "TriggerResultsFilter",
                                         triggerConditions = cms.vstring(*triggers),
                                         hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                         l1tResults = cms.InputTag( "" ),
                                         l1tIgnoreMask = cms.bool( False ),
                                         l1techIgnorePrescales = cms.bool( False ),
                                         daqPartitions = cms.uint32( 1 ),
                                         throw = cms.bool( True )
                                     )


if isMC_:
        process.p = cms.Path(
        process.metaTree
        *process.TriggerSelection
        *process.LeptonsFilter
        *process.egmGsfElectronIDSequence
        *process.jetCorrFactors
        *process.slimmedJetsJEC
        *process.prefiringweight
        #*process.fullPatMetSequenceModifiedMET
        *process.electronMVAValueMapProducer
        *process.btagging
        *process.displacedInclusiveVertexing
        *process.jetSmearing
        *process.jetSmearingUp
        *process.jetSmearingDown
        *process.HeavyNeutralLepton
    )
else :
    process.p = cms.Path(
        process.metaTree
        *process.LeptonsFilter
        *process.egmGsfElectronIDSequence
        *process.jetCorrFactors
        *process.slimmedJetsJEC
        *process.prefiringweight
        *process.electronMVAValueMapProducer
        *process.btagging
        *process.displacedInclusiveVertexing
        *process.HeavyNeutralLepton
    )
