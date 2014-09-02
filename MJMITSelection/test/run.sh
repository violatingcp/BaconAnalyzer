#!/bin/bash

eosdir=$1
gen=$2
name=$3
json=$4
dimu=$5
xs=$6
unit=$7
mkdir $name
echo $eosdir

dir=/store/cmst3/user/pharris/production/03/ntuples/${name}
cmsMkdir ${dir}
ytot=`/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select ls -ltr $eosdir  | wc -l`
echo "tot "$ytot
#tot=`expr $y`
tot=-1
count=-1
outputname=files_0
rm ${name}/${outputname}
for x in `/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select ls  $eosdir `; do
    echo $x
    y=`ls $name/Hist$x | wc -l`
    if [ $y -eq 1 ]; then
	echo " found "$x
	continue
    fi
    if [ $count -eq $unit ]; then
	outputname=files_${tot}
        count=0
        rm ${name}/${outputname}
    fi
    tot=`expr $tot + 1`
    echo $x >> ${name}/${outputname}
    count=`expr $count + 1`
done
for x in `seq 0 $unit $tot`; do
    #echo    runMJMIT -1 root://eoscms.cern.ch/$eosdir/$x $gen $json $dimu $xs 
    #    mv Output.root $name/$x
    #    root -b -q rescale.C\(\"root://eoscms.cern.ch/$eosdir/$x\"\)
    #mv Output.root $name/Hist$x
    ./runBatch.sh root://eoscms.cern.ch/$eosdir $gen $json $dimu $xs $name $x $PWD 
    #bsub -q 1nh  -o out.%J  runBatch.sh root://eoscms.cern.ch/$eosdir $gen $json $dimu $xs $name $x $PWD
    exit
done
