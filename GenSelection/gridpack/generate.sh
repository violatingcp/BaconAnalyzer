#!/bin/bash

dir=$1
nevt=$2
seed=$3
tarball=$4 
LHE=$5
cp $dir/gpacks/$tarball           .
cp $dir/$LHE               .
cp $dir/Hadronizer_TuneCUETP8M1_13TeV_generic_LHE_pythia8_cff.py .
cp $dir/Hadronizer_Tune4C_8TeV_aMCatNLO_LHE_pythia8_cff_Z.py .
#cp $dir/Hadronizer_Tune4C_8TeV_aMCatNLO_LHE_pythia8_cff_A.py .
cp $dir/makeBacon.py .

tar xvf $tarball   > /dev/null
#cd $dir
scramv1 project CMSSW CMSSW_7_1_20
cd CMSSW_7_1_20/src/
eval `scramv1 runtime -sh`
cd -

./runcmsgrid.sh $nevt $seed 1