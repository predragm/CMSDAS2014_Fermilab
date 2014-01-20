#! /bin/bash

#CMSSW 
pwd
export OSG_APP=/raid/osgpg/pg/app
export SCRAM_ARCH=slc5_amd64_gcc472
source $OSG_APP/cmssoft/cms/cmsset_default.sh
eval `scramv1 runtime -sh`
#cmsenv

#export SCRAM_ARCH=slc5_amd64_gcc462
#export SCRAM_ARCH=slc5_amd64_gcc472
#source /afs/cern.ch/cms/cmsset_default.sh
#eval `scramv1 runtime -sh`
#cmsenv

# Run hypothesis testing, using nominal value of nuisances and mu for generation
NTOYS=10000  #toys per  job
MH=126.0     # mass of the signal hypothesis
SEED=1234
CARD=StatCard_hzz4l_Spin0M.txt

### FLOATING MU:
text2workspace.py -m $MH $CARD -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=muFloating -o floatMu.root
combine -m ${MH} -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=1 floatMu.root --singlePoint 1 --saveHybridResult --seed $SEED -T $NTOYS -i 1 --clsAcc 0 --fullBToys
