from CRABClient.UserUtilities import config, ClientException, getUsernameFromSiteDB
from os import *

config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'mc_test'
config.section_('Data')
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/%s/mc_test' % (getUsernameFromSiteDB())
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HeavyNeutralLeptonAnalyzer_cfg.py'
#config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader'
#config.JobType.inputFiles = ['L1PrefiringMaps_new.root']
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' 
#'/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#'/HNL_3GeV_2L2Q_DisplacedSignal/atalierc-processed_gen-sim-premix-reco-miniaod-919c80a76a70185609d372d13ecbc645/USER'
#'/store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/100000/D4B7750A-4D94-E811-B78B-842B2B181788.root'
#'/HNL_3GeV_2L2Q_DisplacedSignal/atalierc-processed_gen-sim-premix-reco-miniaod-919c80a76a70185609d372d13ecbc645/USER'
#'/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#'/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
#'/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#'/HNL_2GeV_2L2Q_DisplacedSignal/atalierc-processed_gen-sim-premix-step2-miniaod-919c80a76a70185609d372d13ecbc645/USER'
#'/HNL_2L2Q_DisplacedSignal/atalierc-processed_gen-sim-premix-RECO-MINIAODv2-10a67796329f4238191918f07d0f7633/USER'
config.JobType.outputFiles = ['Analysis_output.root']
config.JobType.maxJobRuntimeMin = 5500
config.JobType.allowUndistributedCMSSW = True
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
#config.Site.storageSite ='T2_BE_IIHE'



