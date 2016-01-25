#!/bin/bash

dir=$1
nevt=$2
seed=$3
tarball=$4
label=$5
LHE=$6

cp $dir/generate.sh .
./generate.sh $dir $nevt $seed $tarball $LHE
cd CMSSW_7_1_20/src/
eval `scramv1 runtime -sh`
cd -
cmsRun $LHE
cd /afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_7_4_14/src/
eval `scramv1 runtime -sh`
cd -
cmsRun makeBacon.py
option=1
cd /afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_7_4_12_patch1/src/
eval `scramv1 runtime -sh`
cd -

runGen -1  MonoJ.root ${option}  
cmsMkdir               /store/cmst3/user/pharris/mc_v2/${label}
cmsRm                  /store/cmst3/user/pharris/mc_v2/${label}/${label}_${seed}.root
cmsStage Output.root   /store/cmst3/user/pharris/mc_v2/${label}/${label}_${seed}.root
 
