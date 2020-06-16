from os import *
from os import path as path
import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
from RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi import *
from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
from RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi import *
from HNL.HeavyNeutralLeptonAnalysis.ele_Sequence_cff import addElectronSequence


import sys

import re
import importlib
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register(
    'period', 'A',
    VarParsing.multiplicity.singleton, VarParsing.varType.string,
    "Period: default A"
)

options.register(
    'isMC', False,
    VarParsing.multiplicity.singleton, VarParsing.varType.bool,
    "isMC: default False"
)

options.register(
    'isMCSignal', False,
    VarParsing.multiplicity.singleton, VarParsing.varType.bool,
    "isMCSignal: default False"
)

options.register(
    'prescale', -1,
    VarParsing.multiplicity.singleton, VarParsing.varType.int,
    "prescale to split sample"
)

options.register(
    'offset', 0,
    VarParsing.multiplicity.singleton, VarParsing.varType.int,
    "offset to split sample"
)

options.register(
    'hasLHE', True,
    VarParsing.multiplicity.singleton, VarParsing.varType.bool,
    "The sample has LHE products"
)

options.register ('inputFile',
	          '',
              VarParsing.multiplicity.singleton,
              VarParsing.varType.string,
	          info="input file name")

options.register ('newIVF',
               False,
               VarParsing.multiplicity.singleton,
               VarParsing.varType.bool,
               info="switch between the different version of IVF")

options.parseArguments()
#hasLHE_ = options.Flag

period_ = options.period
hasLHE_ = options.hasLHE
debugLevel    = -1
isMC_         = options.isMC
isMCSignal_    = options.isMCSignal
prescale_ = options.prescale
offset_ = options.offset
GT_MC = '102X_upgrade2018_realistic_v18'#94X_mc2017_realistic_v14
edmOut = False


if period_ in 'ABC':
    GT_DATA = '102X_dataRun2_Sep2018ABC_v2'#94X_dataRun2_v6
elif period_ == 'D':
    GT_DATA = '102X_dataRun2_Prompt_v13' #2018D
else:
    raise RuntimeError()

GT      =  GT_MC if isMC_ else GT_DATA

#system('ls -ltr')

process = cms.Process("AnalysisProc")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('HNL.HeavyNeutralLeptonAnalysis.LeptonFilter_cfi')
process.load("CondCore.CondDB.CondDB_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 1

import FWCore.PythonUtilities.LumiList as LumiList
LumiList.LumiList().getVLuminosityBlockRange()

#from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, GT)


##################### prescale #######################################
if prescale_ >= 0:
  process.prescale = cms.EDFilter("EventPrescaler",
                                  prescale=cms.int32(prescale_),
                                  offset=cms.int32(offset_))

###################### input file for testing ##########################
if len(options.inputFile) == 0:
    process.source = cms.Source("PoolSource",
                                fileNames =  cms.untracked.vstring(
                                '/store/mc/RunIIAutumn18MiniAOD/QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v2/110000/99C13D85-30DE-5349-930E-D1662BE02690.root'
                                )
                               )

else:
    process.source = cms.Source("PoolSource",
                                fileNames =  cms.untracked.vstring('file:'+options.inputFile)#'file:/pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-3_V-0.00244948974278_mu_Dirac_massiveAndCKM_LO/heavyNeutrino_1.root')
                                )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
####################################################################


###################### output file #############################
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile))
################################################################


######################  Meta Ntuplizer #########################
process.load('HNL.HeavyNeutralLeptonAnalysis.MetaNtuplizer_cfi')
process.metaTree.isMC = isMC_
process.metaTree.weightsSrc = cms.InputTag('externalLHEProducer')
process.metaTree.globalTag = GT
process.metaTree.args = cms.string('USELESS') #FILL ME!
process.metaTree.hasLHE = cms.bool(hasLHE_ and isMC_)
################################################################



########################### Displaced IVF ######################
process.load('HNL.DisplacedAdaptiveVertexFinder.displacedInclusiveVertexing_cff')
print 'Getting HNL.DisplacedAdaptiveVertexFinder.displacedInclusiveVertexing PSet'
print ' Default useObjectForSeeding.value(): %s' % process.displacedInclusiveVertexFinder.useObjectForSeeding.value()
process.displacedInclusiveVertexFinder.useObjectForSeeding.setValue(options.newIVF)
print ' New useObjectForSeeding.value(): %s' % process.displacedInclusiveVertexFinder.useObjectForSeeding.value()
################################################################


################################ JETMET################################
# update MET according to new JECs and compute uncertainties
# given that to update the MET we need to update the Jets we also take the updated jets from here
#the module below automatically load the new Jet energy corrections
# from: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD( process,isData = not isMC_)

# To get updated ecalBadCalibReducedMINIAODFilter
# See https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
# Recipe is preliminary, i.e. recommended to check for updates

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
  baddetEcal       = baddetEcallist,
  taggingMode      = cms.bool(True),
  debug            = cms.bool(False)
  )


# we redo the smearing because that in runMetCorAndUncFromMiniAOD has some selection (e.g., remove jets overlapping with muons)
process.jetSmearingSeq = cms.Sequence()
if isMC_:
    process.jetSmearing = process.patSmearedJets.clone(
        src = cms.InputTag("patJetsReapplyJEC")
    )
    process.jetSmearingUp   = process.jetSmearing.clone(variation = cms.int32(101))
    process.jetSmearingDown = process.jetSmearing.clone(variation = cms.int32(-101))
    process.jetSmearingSeq = cms.Sequence(process.jetSmearing + process.jetSmearingUp + process.jetSmearingDown)

################################################################################


#declare producer for ecalBadCalibReducedMINIAODFilter
#https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2

# EGamma recipe: run electronID and photonID  https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#2018_Data_MC
################################Electron ID ###########################
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process, runVID=False, era = '2018-Prompt')
#######################################################################

#https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe#Recipe_details_10_2_X_X_10_or_9   Should not be needed

#https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_9_or_10  Should not be needed


process.HeavyNeutralLepton = cms.EDAnalyzer(
    'HeavyNeutralLeptonAnalysis',#HeavyNeutralLeptonAnalysis
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
    jetSrc                = cms.InputTag("patJetsReapplyJEC"),
    pfMETSrc              = cms.InputTag("slimmedMETs","","AnalysisProc"),
    triggerResultSrc      = cms.InputTag("TriggerResults","","HLT"),
    metFilterResultSrc    = cms.InputTag("TriggerResults","","PAT"),
    genParticleSrc        = cms.InputTag("prunedGenParticles"),
    genEventInfoProduct   = cms.InputTag("generator"),
    PUInfo                = cms.InputTag("slimmedAddPileupInfo"),
    lheEventProducts      = cms.InputTag("externalLHEProducer"),
    SecondaryVertices     = cms.InputTag("displacedInclusiveSecondaryVertices"),
    electronsMVA   = cms.string("ElectronMVAEstimatorRun2Fall17NoIsoV1Values"), # TO CHANGE
    electronsVeto  = cms.string("cutBasedElectronID-Fall17-94X-V2-veto"),
    electronsLoose = cms.string("cutBasedElectronID-Fall17-94X-V2-loose"),
    electronsMedium= cms.string("cutBasedElectronID-Fall17-94X-V2-medium"),
    electronsTight = cms.string("cutBasedElectronID-Fall17-94X-V2-tight"),
    jetsSmeared     = cms.InputTag("jetSmearing"     if isMC_ else "patJetsReapplyJEC"),
    jetsSmearedUp   = cms.InputTag("jetSmearingUp"   if isMC_ else "patJetsReapplyJEC"),
    jetsSmearedDown = cms.InputTag("jetSmearingDown" if isMC_ else "patJetsReapplyJEC"),
    electronsEffectiveAreas = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt')
)

process.MessageLogger = cms.Service(
    "MessageLogger",
    suppressWarning= cms.untracked.vstring('displacedInclusiveVertexFinder')
)

process.Timing = cms.Service("Timing",
  summaryOnly = cms.untracked.bool(True),
  useJobReport = cms.untracked.bool(False)
)

triggers = [
    'HLT_IsoMu24_v*',
    'HLT_IsoMu27_v*',
    'HLT_Mu50_v*',
    'HLT_Ele27_WPTight_Gsf_v*',
    'HLT_Ele32_WPTight_Gsf_v*',
    'HLT_Ele115_CaloIdVT_GsfTrkIdT_v*',
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

process.p = cms.Path()
if prescale_ >= 0:  process.p *= process.prescale

process.p *= process.metaTree*process.TriggerSelection*process.LeptonsFilter*process.egammaPostRecoSeq*process.ecalBadCalibReducedMINIAODFilter*process.fullPatMetSequence*process.displacedInclusiveVertexing*process.jetSmearingSeq*process.HeavyNeutralLepton


# crea un output EDM
if edmOut:
        process.edmOut = cms.OutputModule(
            "PoolOutputModule",
            # use this in case of filter available
            outputCommands = cms.untracked.vstring(
                'keep *',
            ),
            fileName = cms.untracked.string('edmTEST.root')
        )
        process.end = cms.EndPath(
            process.edmOut
        )


def diff(mod1, mod2):
    s1 = mod1.__repr__().split('\n')
    s2 = mod2.__repr__().split('\n')
    m_size = max(
        max(len(i) for i in s1),
        max(len(i) for i in s2),
    )
    form = '%-'+str(m_size)+'s'
    for l1, l2 in zip(s1, s2):
        space = ' ' if l1 == l2 else '|'
        print form % l1, space, l2
