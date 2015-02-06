#!/bin/bash

#runMJMIT NEvents file IsGen(0=false,1=true) JSON(dataOnly) Selection(-1=>photon,0=>MET,1=>single & dimuon) XS(sample xs)

runMJMIT 10000 root://eoscms.cern.ch//store/user/ksung/production/phys14test/Phys14-PU20bx25_GJets_HT-600toInf_Tune4C/Phys14-PU20bx25_GJets_HT-600toInf_Tune4C_bacon.root 0 ../data/json/a06.json -1 1
