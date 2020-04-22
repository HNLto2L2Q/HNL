# HNL
Repository for analyzer Heavy Neutral Leptons in CMS

## Setup working area

Instructions:

```bash
cmsrel CMSSW_10_2_15_patch2
cd CMSSW_10_2_15_patch2/src
cmsenv

## Egamma POG
## Ref: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#2018_Data_MC

git cms-init
git cms-merge-topic cms-egamma:EgammaPostRecoTools #just adds in an extra file to have a setup function to make things easier
git cms-merge-topic cms-egamma:PhotonIDValueMapSpeedup1029 #optional but speeds up the photon ID value module so things fun faster
git cms-merge-topic cms-egamma:slava77-btvDictFix_10210 #fixes the Run2018D dictionary issue, see https://github.com/cms-sw/cmssw/issues/26182, may not be necessary for later releases, try it first and see if it works
#now to add the scale and smearing for 2018 (eventually this will not be necessary in later releases but is harmless to do regardless)
git cms-addpkg EgammaAnalysis/ElectronTools
rm EgammaAnalysis/ElectronTools/data -rf
git clone git@github.com:cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data

## MET Unc
## Ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_9_or_10

git cms-merge-topic cms-met:METFixEE2017_949_v2_backport_to_102X


## MET Optional Filter Run2
## Ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM

git cms-addpkg RecoMET/METFilters

## Adding the modify IVF
git cms-addpkg RecoVertex/AdaptiveVertexFinder

#now build everything
scram b -j 8


## HNL2L2Q code

git clone git@github.com:HNLto2L2Q/HNL.git
pushd HNL/
# For 2018
git checkout run_2018
popd

scram b -j 8
```

## Test

```bash
cd src/HNL/HeavyNeutralLeptonAnalysis/test
cmsRun HeavyNeutralLeptonAnalyzer_cfg.py
```
