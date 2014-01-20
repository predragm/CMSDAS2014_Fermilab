#! /bin/bash

### Check input: 
if [ $# -gt 0 ]
    then
    echo "Need no arguments"
    exit 1
fi

# mass of the signal hypothesis
MH=126

#make the tree of the test statistics distribution (the macro is under HiggsAnalysis/CombinedLimit/test/plotting)
root -l -q -b higgsCombineTest.HybridNew.mH${MH}.1234.root "${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTree.cxx(\"qmu.root\",${MH},1,\"x\")"
root -l -q -b "extractJCPHypothesisTestResults.C(\"0M\",true)"
