#!/bin/bash

#runMJMIT NEvents file IsGen(0=false,1=true) JSON(dataOnly) Selection(-1=>photon,0=>MET,1=>single & dimuon) XS(sample xs)

runMJMIT -1 root://eoscms.cern.ch//store/user/ksung/production/phys14test/Phys14-PU20bx25_ZJetsToNuNu_HT-600toInf_Tune4C/Output_99_1_m2V.root 1 ../data/json/a06.json 0 1