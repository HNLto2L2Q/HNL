# HNL
Repository to analyzer Heavy Neutral Leptons in CMS

Instructions:

```
cmsrel CMSSW_9_4_12
cd CMSSW_9_4_12/src
cmsenv
git cms-init

# EGAMMA
# ...
git cms-merge-topic cms-egamma:EgammaID_949 #if you want the FallV2 IDs, otherwise skip
git cms-merge-topic cms-egamma:EgammaPostRecoTools_940 #just adds in an extra file to have a setup function to make things easier

# L1 ECAL Prefiring weights
# twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe#Recipe_details_80X_94X
git cms-merge-topic lathomas:L1Prefiring_9_4_9

# Official Prescription for calculating corrections and uncertainties on Missing Transverse Energy (MET)
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_9_or_10
git cms-merge-topic cms-met:METFixEE2017_949_v2

git clone git@github.com:HNLto2L2Q/HNL.git
cd HNL/
# For 2017
git checkout run_2017


scramv1 b -j 8
```
