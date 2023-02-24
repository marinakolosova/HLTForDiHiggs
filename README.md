# HLTForDiHiggs
Setup the HLT environment from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGlobalHLT#CMSSW_13_0_X_presently_used_for
```
cmsrel CMSSW_13_0_0_pre4
cd CMSSW_13_0_0_pre4/src
cmsenv
git cms-init
scram build -j 4
```
and then clone the repo:
```
git clone git@github.com:marinakolosova/HLTForDiHiggs.git
scram build -j 4
```
