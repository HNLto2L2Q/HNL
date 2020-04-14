from CRABClient.UserUtilities import config, ClientException, getUsernameFromSiteDB, getUsernameFromCRIC
#from input_crab_data import dataset_files
import os
import yaml
from datetime import datetime

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--ymls", type=str, nargs='*', required=True, help="yml files" )
parser.add_argument("--cfgfile", type=str, required=True, help="path to cfgfile to be used" )
parser.add_argument("--newIVF", action='store_true', default=False, help="Process events with the modify IVF")
parser.add_argument("--dryRun", action='store_true', default=False, help="CRAB submit dryrun")
parser.add_argument("--types", type=str, choices=["data", "background", "signal"], nargs='*',
                               default=["data", "background", "signal"],
                               help="sample types [data, background, signal]" )
parser.add_argument("--storage_site", type=str, required=True, help="Storage site where you are allow to write")
parser.add_argument("--groups", type=str, choices=["DY","DY_MLM","ttbar","WJets","WW", "WZ", "ZZ", "VGamma", "SingleTop"], nargs='*',
                                help="Sample groups to process [DY, DY_MLM, ttbar, WJets, WW, WZ, ZZ, VGamma, SingleTop]" )

now  = datetime.now()
time = now.strftime("%y%m%d")

args         = parser.parse_args()
ymls         = args.ymls
types        = args.types
cfgfile      = os.path.abspath(args.cfgfile)
newIVF       = args.newIVF
dryRun       = args.dryRun
storage_site = args.storage_site
groups       = args.groups

tag = "newIVF_{}".format(time) if newIVF else "usualIVF_{}".format(time)

os.system("voms-proxy-init --voms cms")

config = config()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.section_('Data')
config.Data.publication = False
config.Data.inputDBS = 'global'
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = cfgfile
config.JobType.maxJobRuntimeMin = 3000
config.JobType.allowUndistributedCMSSW = True
config.section_('User')
config.section_('Site')
config.Site.storageSite = storage_site


if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    from multiprocessing import Process

#    for dataset, infos in dataset_files.items():
    for yml in ymls:
        with open(yml,'r') as f:
            doc = yaml.load(f)
            sample_types = doc.keys()
            for type in sample_types:
                if type not in types:
                    continue

                config.General.workArea = 'HNL_{type}_2017_{tag}'.format(type=type, tag=tag)
                config.Data.outLFNDirBase = '/store/user/{username}/HNL_{type}_2017_{tag}'.format(username=getUsernameFromCRIC(),
                                                                                                    type=type, tag=tag)
                config.Data.splitting = 'LumiBased' if (type == 'data') else 'FileBased'

                type_dict = doc[type]
                samples  = type_dict['samples'].keys()
                cfg_args = type_dict['cfg_args']
                if newIVF:
                    cfg_args.append('newIVF=True')

                config.JobType.pyCfgParams = cfg_args

                for sample in samples:
                    name       = sample
                    unitPerJob = type_dict['samples'][sample]['unit_per_job']
                    dataset    = type_dict['samples'][sample]['dsname']
                    group      = type_dict['samples'][sample]['group']

                    if len(groups) > 0 :
                        if ("group" in type_dict['samples'][sample]):
                            if not groups.count(group):
                                continue
                        else:
                            continue

                    print dataset, unitPerJob, name
                    print config.JobType.pyCfgParams
                    print config.General.workArea
                    config.Data.inputDataset = dataset
                    config.General.requestName = name
                    config.Data.unitsPerJob = int(unitPerJob)

                    if type == 'data':
                        lumi_file = type_dict['samples'][sample]['lumi_file']
                        config.Data.lumiMask = lumi_file
                        print '\t '+lumi_file

                    crabCommand('submit', config = config, dryrun = dryRun)
