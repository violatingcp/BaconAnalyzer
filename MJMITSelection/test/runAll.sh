#!/bin/bash

###Gamma+Jets
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_G_Pt-50to80/                                    2 G050080     ../data/json/a06.json -1 3322.309     1
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_G_Pt-80to120/                                   2 G080120     ../data/json/a06.json -1 558.2865     1    
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_G_Pt-120to170/                                  2 G120170     ../data/json/a06.json -1 108.0068     1
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_G_Pt-170to300/                                  2 G170300     ../data/json/a06.json -1 30.12207     1
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_G_Pt-300to470/                                  2 G300470     ../data/json/a06.json -1 2.138632     1
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_G_Pt-470to800/                                  2 G470800     ../data/json/a06.json -1 0.2119244    1
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_G_Pt-800to1400/                                 2 G8001400    ../data/json/a06.json -1 0.007077847  1
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_G_Pt-1400to1800/                                2 G14001800   ../data/json/a06.json -1 0.000045103  1
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_G_Pt-1800/                                      2 G1800       ../data/json/a06.json -1 0.000001867  1
##GJ Data
#./run.sh /store/cmst3/group/monojet/production/03/Photon_2012A-22Jan2013/             0 PhoA     ../data/json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt -1 1 1
#./run.sh /store/cmst3/group/monojet/production/03/Photon_2012B-22Jan2013/             0 PhoB     ../data/json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt -1 1 1
#./run.sh /store/cmst3/group/monojet/production/03/SinglePhoton_2012C-22Jan2013/       0 PhoC     ../data/json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt -1 1 1
#./run.sh /store/cmst3/group/monojet/production/03/SinglePhotonParked_2012D-22Jan2013/ 0 PhoD     ../data/json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt -1 1 1

###MET Selection 
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_GluGlu_HToInvisible_M-125/                     1  GGH         ../data/json/a06.json 0 19.27   1
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_VBF_HToInvisible_M-125/                        1  VBF         ../data/json/a06.json 0 1.578   1
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_WH_HToMuMu_M-125/                              1  WH          ../data/json/a06.json 1 0.7046  1
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_ZH_HToMuMu_M-125/                              1  ZH          ../data/json/a06.json 1 0.4153  1
#./run.sh  /store/cmst3/group/monojet/production/03/Summer13-ZH_ZToLL_HToBB_M-125_8T/                       1  ZLLHBB      ../data/json/a06.json 1 0.08306 20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer13-ZH_ZToNuNu_HToBB_M-125/                        1  ZNNHBB      ../data/json/a06.json 0 0.08306 20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer13-WW_TuneZ2star_8TeV/                            1  WW          ../data/json/a06.json 0 57.2    20  
#./run.sh  /store/cmst3/group/monojet/production/03/Summer13-WZ_TuneZ2star_8TeV/                            1  WZ          ../data/json/a06.json 0 22.44   20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer13-ZZ_TuneZ2star_8TeV/                            1  ZZ          ../data/json/a06.json 0 9.03    20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer13-T_s-channel_TuneZ2star_8TeV                    1  Ts          ../data/json/a06.json 0 3.89    20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer13-T_t-channel_TuneZ2star_8TeV                    1  Tt          ../data/json/a06.json 0 55.531  20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer13-T_tW-channel-DR_TuneZ2star_8TeV                1  TW          ../data/json/a06.json 0 11.1    20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer13-TTJets_SemiLeptMGDecays_8TeV                   1  TTSem   ../data/json/a06.json 0 109.281 20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer13-TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Bacon 1  TTLep   ../data/json/a06.json 0 26.1975 20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_ZJetsToNuNu_50_HT_100                          1  DYNN050100  ../data/json/a06.json 0 381.2   1
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_ZJetsToNuNu_100_HT_200                         1  DYNN100200  ../data/json/a06.json 0 160.3   1
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_JetsToNuNu_200_HT_400                          1  DYNN200400  ../data/json/a06.json 0 41.49   1
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_ZJetsToNuNu_400_HT_inf                         1  DYNN400inf  ../data/json/a06.json 0 5.274   1
#./run.sh  /store/cmst3/group/monojet/production/03/Sumer13-WJetsToLNu_TuneZ2Star_8TeV                      1  W           ../data/json/a06.json 0 30400.0 20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer13-WJetsToLNu_HT-250To300_8TeV/                   1  W250300     ../data/json/a06.json 0 48.01   20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer13-WJetsToLNu_HT-300To400_8TeV/                   1  W300400     ../data/json/a06.json 0 38.3    20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer13-WJetsToLNu_HT-400Toinf_8TeV/                   1  W400inf     ../data/json/a06.json 0 25.22   20
#exit
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_DYJetsToLL_M-50                                1  RDY          ../data/json/a06.json 0 3354    20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_DYJetsToLL_PtZ-50To70_ext                      1  RDY5070      ../data/json/a06.json 0 52.31   20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_DYJetsToLL_PtZ-70To100_ext                     1  RDY70100     ../data/json/a06.json 0 34.1    20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_DYJetsToLL_PtZ-100_ext                         1  RDY100       ../data/json/a06.json 0 93.8    20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer13-DYJetsToLL_PtZ-50To70                          1  RDY5070      ../data/json/a06.json 0 52.31   20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer13-DYJetsToLL_PtZ-70To100                         1  RDY70100     ../data/json/a06.json 0 34.1    20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer13-DYJetsToLL_PtZ-100                             1  RDY100       ../data/json/a06.json 0 93.8    20

##MET Data
./run.sh  /store/cmst3/group/monojet/production/03/MET_2012A-22Jan2013/          0 META     ../data/json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt 0 1 1
./run.sh  /store/cmst3/group/monojet/production/03/MET_2012B-22Jan2013/          0 METB     ../data/json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt 0 1 1
./run.sh  /store/cmst3/group/monojet/production/03/MET_2012C-22Jan2013/          0 METC     ../data/json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt 0 1 20
./run.sh  /store/cmst3/group/monojet/production/03/METParked_2012D-22Jan2013/    0 METD     ../data/json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt 0 1 20

##Lepton Selection
./run.sh  /store/cmst3/group/monojet/production/03/Sumer13-WJetsToLNu_TuneZ2Star_8TeV                    1 WL           ../data/json/a06.json 2 30400.0 20
./run.sh  /store/cmst3/group/monojet/production/03/Summer13-WJetsToLNu_HT-250To300_8TeV/                 1 W250300L     ../data/json/a06.json 2 48.01   20
./run.sh  /store/cmst3/group/monojet/production/03/Summer13-WJetsToLNu_HT-300To400_8TeV/                 1 W300400L     ../data/json/a06.json 2 38.3    20
./run.sh  /store/cmst3/group/monojet/production/03/Summer13-WJetsToLNu_HT-400Toinf_8TeV/                 1 W400infL     ../data/json/a06.json 2 25.22   20
./run.sh  /store/cmst3/group/monojet/production/03/Summer13-DYJetsToLL_M-50_TuneZ2Star_8TeV/             1 DYL          ../data/json/a06.json 2 3354    20
./run.sh  /store/cmst3/group/monojet/production/03/Summer13-DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV        1 DY5070L      ../data/json/a06.json 2 52.31   20
./run.sh  /store/cmst3/group/monojet/production/03/Summer13-DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV       1 DY70100L     ../data/json/a06.json 2 34.1    20
./run.sh  /store/cmst3/group/monojet/production/03/Summer13-DYJetsToLL_PtZ-100_TuneZ2star_8TeV           1 DY100L       ../data/json/a06.json 2 93.8    20
./run.sh  /store/cmst3/group/monojet/production/03/Summer13-TTJets_SemiLeptMGDecays_8TeV                 1 TTSemL ../data/json/a06.json 2 109.281 20
./run.sh  /store/cmst3/group/monojet/production/03/Summer13-TTJets_FullLeptMGDecays_8TeV-madgraph-tauola_Bacon 1 TTLepL ../data/json/a06.json 2 26.1975 20
./run.sh  /store/cmst3/group/monojet/production/03/Summer13-T_s-channel_TuneZ2star_8TeV                  1 TsL          ../data/json/a06.json 2 3.89    20
./run.sh  /store/cmst3/group/monojet/production/03/Summer13-T_t-channel_TuneZ2star_8TeV                  1 TtL          ../data/json/a06.json 2 55.531  20
./run.sh  /store/cmst3/group/monojet/production/03/Summer13-T_tW-channel-DR_TuneZ2star_8TeV              1 TWL          ../data/json/a06.json 2 11.1    20        
./run.sh  /store/cmst3/group/monojet/production/03/Summer13-WW_TuneZ2star_8TeV/                          1 WWL          ../data/json/a06.json 2 57.2    20 
./run.sh  /store/cmst3/group/monojet/production/03/Summer13-WZ_TuneZ2star_8TeV/                          1 WZL          ../data/json/a06.json 2 22.44   20   
./run.sh  /store/cmst3/group/monojet/production/03/Summer13-ZZ_TuneZ2star_8TeV/                          1 ZZL          ../data/json/a06.json 2 9.03    20  
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_DYJetsToLL_PtZ-50To70_ext                    1 DY5070L      ../data/json/a06.json 2 52.31   20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_DYJetsToLL_PtZ-70To100_ext                   1 DY70100L     ../data/json/a06.json 2 34.1    20
#./run.sh  /store/cmst3/group/monojet/production/03/Summer12_DYJetsToLL_PtZ-100_ext                       1 DY100L       ../data/json/a06.json 2 93.8    20

#Lepton Data
./run.sh  /store/cmst3/group/monojet/production/03/MET_2012A-22Jan2013/          0 METAL     ../data/json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt 2 1 1
./run.sh  /store/cmst3/group/monojet/production/03/MET_2012B-22Jan2013/          0 METBL     ../data/json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt 2 1 1
./run.sh  /store/cmst3/group/monojet/production/03/MET_2012C-22Jan2013/          0 METCL     ../data/json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt 2 1 20
./run.sh  /store/cmst3/group/monojet/production/03/METParked_2012D-22Jan2013/    0 METDL     ../data/json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt 2 1 20



