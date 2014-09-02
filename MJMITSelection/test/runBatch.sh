#!/bin/bash

filedir=$1
gen=$2
json=$3
dimu=$4
xs=$5
name=$6
begin=$7
basedir=$8

cp $basedir/$name/files_${begin} .

cd $basedir
eval `scramv1 runtime -sh`
cd -

rm *.root
cp $basedir/rescale.C .
for x in `cat files_${begin}`; do 
    runMJMIT -1 $filedir/$x $gen $basedir/$json $dimu $xs
    mv Output.root $x
    
    root -b -q rescale.C\(\"${filedir}/$x\"\)
    mv Output.root $basedir/$name/Hist$x
    filename=$x
done

hadd Output.root *.root 
#cmsRm  /store/group/phys_jetmet/pharris/production/01/ntuples/${name}/${filename}
#cmsStage Output.root /store/group/phys_jetmet/pharris/production/01/ntuples/${name}/${filename}
mv Output.root $basedir/$name/$filename

