from CRABClient.UserUtilities import config, ClientException, getUsernameFromSiteDB
#from input_crab_data import dataset_files
import yaml

config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'HNL_data_2017_corrected'
config.section_('Data')
#config.Data.splitting = 'FileBased'
config.Data.splitting = 'LumiBased'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/%s/HNL_data_2017_corrected' % (getUsernameFromSiteDB())
config.Data.inputDBS = 'global'
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HeavyNeutralLeptonAnalyzer_cfg.py'
config.JobType.maxJobRuntimeMin = 3000
config.JobType.allowUndistributedCMSSW = True
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_IT_Bari'
config.JobType.pyCfgParams = ['isMC=False']

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    from multiprocessing import Process

#    for dataset, infos in dataset_files.items():
    with open('input_data_2017.yml','r') as f:
       doc = yaml.load(f)
       samples = doc['samples'].keys()
       for sample in samples:
         name       = sample
         unitPerJob = doc['samples'][sample]['unit_per_job'] 
         dataset    = doc['samples'][sample]['dsname']
         group      = doc['samples'][sample]['group']
         print dataset, unitPerJob, name
         config.Data.inputDataset = dataset
         config.General.requestName = name
         config.Data.unitsPerJob = int(unitPerJob)
         
         if group == 'data':
           lumi_file = doc['samples'][sample]['lumi_file']
           config.Data.lumiMask = lumi_file
           print '\t '+lumi_file
         
         crabCommand('submit', config = config, dryrun = True)


