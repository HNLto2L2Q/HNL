from CRABClient.UserUtilities import config, ClientException, getUsernameFromCRIC
#from input_crab_data import dataset_files
import yaml

config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'HNL_MC_2018'
config.section_('Data')
config.Data.splitting = 'FileBased'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/%s/HNL_MC_2018' % (getUsernameFromCRIC())
config.Data.inputDBS = 'global'
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HeavyNeutralLeptonAnalyzer_cfg.py'
config.JobType.maxJobRuntimeMin = 3000
config.JobType.allowUndistributedCMSSW = True
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_BE_UCL'

from pdb import set_trace
if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    from multiprocessing import Process

 
#    for dataset, infos in dataset_files.items():
    with open('input_mc_2018.yml','r') as f:
       doc = yaml.load(f)
       samples = doc['samples'].keys()
       for sample in samples:
         if 'QCD' not in sample: continue
         name       = sample
         unitPerJob = doc['samples'][sample]['unit_per_job'] 
         dataset    = doc['samples'][sample]['dsname']
         group      = doc['samples'][sample]['group']
         prescaleValue = int(doc['samples'][sample]['prescale'])
         hasLHE = bool(doc['samples'][sample].get('hasLHE', 1))

         print dataset, unitPerJob, name, prescaleValue
         config.Data.inputDataset = dataset
         config.General.requestName = name
         config.Data.unitsPerJob = int(unitPerJob)
         config.JobType.pyCfgParams = ['isMC=True', 'hasLHE=%s' % hasLHE]
         
         if prescaleValue > 0:
             print 'further splitting...'
             offsetValue = [i for i in range(prescaleValue)]
             for i in offsetValue: 
                 print "offset is: ", i
                 config.Data.inputDataset = dataset
                 config.General.requestName = name + '_'  + str(i)
                 config.Data.unitsPerJob = int(unitPerJob)
                 config.JobType.pyCfgParams = ['isMC=True','prescale=%s'%prescaleValue, 'offset=%s'%i, 'hasLHE=%s' % hasLHE]
                 crabCommand('submit', config = config, dryrun = False)


         else:
             crabCommand('submit', config = config, dryrun = False)


