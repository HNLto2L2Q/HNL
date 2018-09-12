import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
def addElectronSequence(process):
  
  process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
  switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
  electronModules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff']
  for idmod in electronModules: setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
  process.ele_Sequence = cms.Sequence(process.egmGsfElectronIDSequence)
                     
