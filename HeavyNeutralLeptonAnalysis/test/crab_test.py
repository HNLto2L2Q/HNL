from CRABClient.UserUtilities import config, ClientException, getUsernameFromSiteDB

config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'prova_2017_signalworong'
config.section_('Data')
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HeavyNeutralLeptonAnalyzer_cfg.py'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader'
config.Data.inputDataset = '/HNL_2GeV_2L2Q_DisplacedSignal/atalierc-processed_gen-sim-premix-step2-miniaod-919c80a76a70185609d372d13ecbc645/USER'
#'/HNL_2L2Q_DisplacedSignal/atalierc-processed_gen-sim-premix-RECO-MINIAODv2-10a67796329f4238191918f07d0f7633/USER'
config.JobType.outputFiles = ['Analysis_output.root']
config.JobType.maxJobRuntimeMin = 5500
config.JobType.allowUndistributedCMSSW = True
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'




