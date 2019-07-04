#!/bin/bash

export FULLJRDIR=$PWD
#export FASTJETINSTALLPATH="/Users/cfmcginn/Packages/FastJet/"
export FASTJETINSTALLPATH="/afs/cern.ch/work/c/cmcginn/public/Fastjet/"

#export ROOUNFPATH=/Users/cfmcginn/Packages/RooUnfold/RooUnfold-1.1.1
#export ROOUNFPATH=/afs/cern.ch/work/c/cmcginn/public/Packages/RooUnfold/RooUnfold_20180820_ROOT6_15_01/

#This root version is from CMSSW_9_4_8
#export ROOUNFPATH=/afs/cern.ch/work/c/cmcginn/public/Packages/RooUnfold/RooUnfold_20180820_ROOT6_10_09/

#This root version is latest/greatest i can get to build
export ROOUNFPATH=/afs/cern.ch/work/c/cmcginn/public/Packages/RooUnfold/RooUnfold_20190103_ROOT6_12_07_CMSSW1031p2/

echo "SETTING NEW ENVIRONMENT VARIABLE FULLJRDIR AS '$PWD'"

echo "All environment variables set"
