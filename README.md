# HNL
Repository for analyzer Heavy Neutral Leptons in CMS

| Update date | Updates |
| ----------- | ------ |
| 22-3-19 | EGamma recipes |

## Setup working area

Instructions:

```bash
export SCRAM_ARCH=slc6_amd64_gcc630
cmsrel CMSSW_9_4_13
cd CMSSW_9_4_13/src
cmsenv
git cms-init

# EGAMMA
# Energy corrections
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#2017_MiniAOD_V2
# Post RECO tools:
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#2016_2017_Data_MC
git cms-merge-topic cms-egamma:EgammaPostRecoTools #just adds in an extra file to have a setup function to make things easier

# L1 ECAL Prefiring weights
# twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe#Recipe_details_80X_94X
git cms-merge-topic lathomas:L1Prefiring_9_4_9
# Moving L1PrefiringMaps in data folder to make it available on CRAB
pushd L1Prefiring/EventWeightProducer/
mkdir data
cp files/L1PrefiringMaps_new.root data/
popd

# Official Prescription for calculating corrections and uncertainties on Missing Transverse Energy (MET)
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_9_or_10
git cms-merge-topic cms-met:METFixEE2017_949_v2

#PLEASE CHECK THE BAD CHANNEL CORRECTION
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
#they are included in the final cgf but they are 'WORK IN PROGRESS' in the twiki

## Adding the modify IVF
git cms-merge-topic clacaputo:HNL_IVF_Mods

git clone git@github.com:HNLto2L2Q/HNL.git
pushd HNL/
# For 2017
git checkout run_2017
popd

scram b -j 8
```

## Test

```bash
cd src/HNL/HeavyNeutralLeptonAnalysis/test
cmsRun HeavyNeutralLeptonAnalyzer_cfg.py
```
