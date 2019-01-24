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

git clone git@github.com:HNLto2L2Q/HNL.git
cd HNL/
# For 2017
git checkout run_2017


scramv1 b -j 8
```
