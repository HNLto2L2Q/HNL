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

obj = VarParsing ('analysis')
# outputFile already defined in VarParsing
#obj.register ('outputFile',
#              'Analysis_output',
#               VarParsing.multiplicity.singleton,
#               VarParsing.varType.string,
#               info="input file name")
obj.register ('inputFile',
	          '',
              VarParsing.multiplicity.singleton,
              VarParsing.varType.string,
	          info="input file name")
obj.register ('newIVF',
               False,
               VarParsing.multiplicity.singleton,
               VarParsing.varType.bool,
               info="switch between the different version of IVF")
# get and parse the command line arguments
obj.parseArguments()

hasLHE_ = True

debugLevel    = -1

isMC_         = True
isMCSignal_    = True
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

process.MessageLogger.cerr.FwkReport.reportEvery = 1

import FWCore.PythonUtilities.LumiList as LumiList
LumiList.LumiList().getVLuminosityBlockRange()
#LumiList.LumiList(filename = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt').getVLuminosityBlockRange()

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
    input = cms.untracked.int32(-1)
    )


process.source = cms.Source("PoolSource",
                            fileNames =  cms.untracked.vstring('file:'+obj.inputFile)#'file:/pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-3_V-0.00244948974278_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_1.root')
                            )


process.TFileService = cms.Service("TFileService", fileName = cms.string(obj.outputFile))
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

#bad channels correction TO BE CHECK
process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313])


process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist,
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
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
                                            pfMETSrc              = cms.InputTag("slimmedMETsModifiedMET" if isMC_ else "slimmedMETs"),
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

process.Timing = cms.Service("Timing",
  summaryOnly = cms.untracked.bool(True),
  useJobReport = cms.untracked.bool(False)
)

print 'Getting HNL.DisplacedAdaptiveVertexFinder.displacedInclusiveVertexing PSet'
print ' Default useObjectForSeeding.value(): %s' % process.displacedInclusiveVertexFinder.useObjectForSeeding.value()
process.displacedInclusiveVertexFinder.useObjectForSeeding.setValue(obj.newIVF)
print ' New useObjectForSeeding.value(): %s' % process.displacedInclusiveVertexFinder.useObjectForSeeding.value()

if (isMC_):
    process.p = cms.Path(
       # process.metaTree
       # process.LeptonsFilter
        process.ecalBadCalibReducedMINIAODFilter
        *process.egmGsfElectronIDSequence
        *process.prefiringweight
        *process.fullPatMetSequenceModifiedMET
        *process.electronMVAValueMapProducer
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
        #process.LeptonsFilter
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
        process.ecalBadCalibReducedMINIAODFilter
        *process.egmGsfElectronIDSequence
        *process.jetCorrFactors
        *process.slimmedJetsJEC
        *process.prefiringweight
        *process.electronMVAValueMapProducer
        *process.displacedInclusiveVertexing
        *process.HeavyNeutralLepton
       )
