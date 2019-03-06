from CRABClient.UserUtilities import config, ClientException, getUsernameFromSiteDB
from input_crab import dataset_files

config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'HNL_mc_2017'
config.section_('Data')
config.Data.splitting = 'FileBased'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HeavyNeutralLeptonAnalyzer_cfg.py'
config.JobType.maxJobRuntimeMin = 3000
config.JobType.allowUndistributedCMSSW = True
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'

if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    from multiprocessing import Process
        
    for dataset, nfiles in dataset_files.items():        
        
        print dataset, nfiles
        config.Data.inputDataset = dataset
        config.General.requestName = dataset.split('/')[1]
        config.Data.unitsPerJob = nfiles
        crabCommand('submit', config = config)#, dryrun = True)
        
        
        
