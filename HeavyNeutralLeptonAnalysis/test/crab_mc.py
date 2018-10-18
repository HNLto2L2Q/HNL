from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'Background/samples/'


config.section_('JobType')
config.JobType.psetName = 'HeavyNeutralLeptonAnalyzer_cfg.py'
config.JobType.pluginName = 'Analysis'
#config.JobType.inputFiles = ['MC_Recent_25ns_2015.root','pileUpData_fromJson.root'] #data files for PileUp reweighting
config.JobType.inputFiles = ['MCpileUp_25ns_Recent2016.root','puData_2016_central.root', 'puData_2016_up.root', 'puData_2016_down.root', 'Summer16_23Sep2016V4_MC.db'] #data files for PileUp reweighting
config.JobType.outputFiles = ['Analysis_output.root']

config.section_('Data')
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/group/phys_exotica/HNL/Background' #% (getUsernameFromSiteDB())                           
config.Data.publishWithGroupName = True
config.Data.publication = False
#config.Data.ignoreLocality = True

config.section_('Site')
#config.Site.storageSite = 'T2_IT_Legnaro'
#config.Site.storageSite = 'T2_US_Nebraska' #write your storage site here                                                                            
config.Site.storageSite = 'T2_BE_IIHE'

#####################################################################################
##################################  IMPORTANT  ######################################
####################  Replace the line runningOnData=cms.bool(True)  ################
#########################  with runningOnData=cms.bool(False) ####################### 
################### in python/Analysis_cfi.py before submitting the job #############
#####################################################################################

if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    from multiprocessing import Process
    
    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)
            
    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

#

config.General.requestName = 'Analysis_ttbar_DiLepton'
config.Data.unitsPerJob = 1
config.Data.inputDataset = '/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()


#
config.General.requestName = 'Analysis_ttbar_DiLepton_ext1'
config.Data.unitsPerJob = 1
config.Data.inputDataset = '/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()



#
config.General.requestName = 'Analysis_ttbar_SingleLeptFromT'
config.Data.unitsPerJob = 1
config.Data.inputDataset = '/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()


#
config.General.requestName = 'Analysis_ttbar_SingleLeptFromT_ext1'
config.Data.unitsPerJob = 1
config.Data.inputDataset = '/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()

config.General.requestName = 'Analysis_ttbar_SingleLeptFromTbar'
config.Data.unitsPerJob = 1 
config.Data.inputDataset = '/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
p = Process(target=submit, args=(config,)) 
p.start() 
p.join()

#
config.General.requestName = 'Analysis_ttbar_SingleLeptFromTbar_ext1'
config.Data.unitsPerJob = 1
config.Data.inputDataset = '/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()


config.General.requestName = 'Analysis_SinglTop_s_channel_leptonDecays'
config.Data.unitsPerJob = 1
config.Data.inputDataset = '/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()


#
config.General.requestName = 'Analysis_SingleTop_tW_top_inclusivedecay'
config.Data.unitsPerJob = 1
config.Data.inputDataset = '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()
#
config.General.requestName = 'Analysis_SingleTop_tW_antitop_inclusivedecay'
config.Data.unitsPerJob = 1
config.Data.inputDataset = '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
p = Process(target=submit, args=(config,))
p.start()
p.join()




