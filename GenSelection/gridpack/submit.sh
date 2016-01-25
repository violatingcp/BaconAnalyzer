#!/bin/bash

label=$1
LHE=$2
dir=`pwd`
nevt=120000
tarball=$3
cmsMkdir /store/cmst3/user/pharris/mc/$label
for x in `seq 200 1200`; do 
    bsub -q 2nd -o out.%J runbatch.sh $dir $nevt $x $tarball $label $LHE
done
