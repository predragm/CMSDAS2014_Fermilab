#! /bin/bash

### Inputs: 
# 1: directory with cards
# 2: card name (signal model1 + signal model2)
# 3: action type
if [ $# -lt 2 ]
    then
    echo "Need at least two arguments: "
    echo "    1) directory with cards"
    echo "    2) card name (signal models)"
    echo "    3) action type (not mandatory, default=1)"
    exit 0
fi

cardDir=$1
card1=$2

cd $cardDir
outDir="output_combine/"

if [ -d $outDir ]
    then
    echo "Output directory ${cardDir}/${outDir}/ already exisiting. I will not overwrite. Please remove it and try again."
    #exit 2
fi

mkdir $outDir

action=1
if [ $# -ge 3 ]
    then
    action=$3
fi

# mass of the signal hypothesis
MH=126

if [ $action -eq 1 ]
    then 
### FIXED MU: 
    text2workspace.py -m $MH $card1 -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs -o fixedMu.root
    combine -M MultiDimFit fixedMu.root --algo=grid --points 100  -m $MH -n 1D --setPhysicsModelParameters x=0 --toysFreq -t -1
elif [ $action -eq 2 ]
    then
### FLOAT MU:
    text2workspace.py -m $MH $card1 -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=muFloating -o floatMu.root
    combine -M MultiDimFit floatMu.root --algo=grid --points 100  -m $MH -n 1D --setPhysicsModelParameters x=0 --toysFreq -t -1
### to draw the output: limit->Draw("2*deltaNLL:x", "deltaNLL > 0","PL")
else
    echo "Requested to perform and unrecognized action: "${action}
    echo "action can be 1:1D scan, fixed mu  ;   2:1D scan, mu floated"
    echo "Exiting."
    exit 3
fi
