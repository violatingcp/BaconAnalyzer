# BaconAnalyzer
#Analysis code for Bacon

cd $CMSSW_BASE/src

git clone https://github.com/violatingcp/BaconAnalyzer.git

scram b

runMJMIT -1 $file $isgen $json $dimu $xs

dimu=0,2 (no muons, two muons)
