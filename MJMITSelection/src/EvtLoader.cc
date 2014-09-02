#include "TFile.h"
#include "TMatrixD.h"
#include "../include/EvtLoader.hh"
#include <iostream>

using namespace baconhep;

EvtLoader::EvtLoader(TTree *iTree,std::string iName,std::string iHLTFile,std::string iPUWeight) { 
  fEvt      = new TEventInfo();
  iTree->SetBranchAddress("Info",       &fEvt);
  fEvtBr    = iTree->GetBranch("Info");
  fTrigger  = new TTrigger(iHLTFile);
  
  fVertices = new TClonesArray("baconhep::TVertex");
  iTree->SetBranchAddress("PV",       &fVertices);
  fVertexBr   = iTree->GetBranch("PV");
  
  TFile *lFile = new TFile(iPUWeight.c_str()); 
  fPUWeightHist = (TH1F*) lFile->FindObjectAny("pileup");
  
  fSample = (char*) (iName.c_str());
}
EvtLoader::~EvtLoader() { 
  delete  fEvt;
  delete  fEvtBr;
  delete  fVertices;
  delete  fVertexBr;
}
void EvtLoader::reset() { 
  fRun       = 0;
  fEvtV      = 0; 
  fLumi      = 0; 
  fRho       = 0; 
  fITrigger  = 0; 
  fHLTMatch  = 0; 
  fNVtx      = 0; 
  fNPU       = 0;
  fNPUP      = 0;
  fNPUM      = 0;
  fPUWeight  = 0; 
  
  fRawSumEt  = 0;
  fMet       = 0; 
  fMetPhi    = 0; 
  fMetRaw    = 0; 
  fMetRawPhi = 0; 

  fMetCov00  = 0; 
  fMetCov01  = 0; 
  fMetCov10  = 0; 
  fMetCov11  = 0; 
  fMetSig    = 0;

  fTKMet     = 0; 
  fTKMetPhi  = 0; 
  fMVAMet    = 0; 
  fMVAMetPhi = 0; 
  
  fMVAMetUnity    = 0; 
  fMVAMetUnityPhi = 0; 

  fMVAMetCov00 = 0; 
  fMVAMetCov01 = 0; 
  fMVAMetCov10 = 0; 
  fMVAMetCov11 = 0; 
  fMVAMetSig   = 0;

  fMVAMetUnityCov00 = 0; 
  fMVAMetUnityCov01 = 0; 
  fMVAMetUnityCov10 = 0; 
  fMVAMetUnityCov11 = 0; 
  fMVAMetUnitySig   = 0;

  fMtTrue         = 0;
  fRawMtTrue      = 0;
  fTKMtTrue       = 0;
  fMVAMtTrue      = 0;
  fMVAMtUnityTrue = 0;
}
void EvtLoader::setupTree(TTree *iTree,float iWeight) { 
  reset();
  fTree = iTree;
  fTree->Branch("sample"        ,fSample         ,"fSample/C",1024);
  fTree->Branch("run"           ,&fRun           ,"fRun/i");
  fTree->Branch("lumi"          ,&fLumi          ,"fLumi/i");
  fTree->Branch("event"         ,&fEvtV          ,"fEvtV/i");
  fTree->Branch("trigger"       ,&fITrigger      ,"fITrigger/i");
  fTree->Branch("hltmatch"      ,&fHLTMatch      ,"fHLTMatch/i");
  fTree->Branch("puweight"      ,&fPUWeight      ,"fPUWeight/F");
  fTree->Branch("npu"           ,&fNPU           ,"fNPU/i");
  fTree->Branch("npuPlusOne"    ,&fNPUP          ,"fNPUP/i");
  fTree->Branch("npuMinusOne"   ,&fNPUM          ,"fNPUM/i");
  fTree->Branch("nvtx"          ,&fNVtx          ,"fNVtx/i");
  fTree->Branch("metFiltersWord",&fMetFilters    , "fMetFilters/I");
  fTree->Branch("scale1fb"      ,&fScale         ,"fScale/F");  fScale = iWeight;
  fTree->Branch("rho"           ,&fRho           ,"fRho/F");
  fTree->Branch("sumet"         ,&fRawSumEt      ,"fRawSumEt/F");

  fTree->Branch("metRaw"        ,&fMetRaw        ,"fMetRaw/F");
  fTree->Branch("metRawPhi"     ,&fMetRawPhi     ,"fMetRawPhi/F");
  fTree->Branch("met"           ,&fMet           ,"fMet/F");
  fTree->Branch("metphi"        ,&fMetPhi        ,"fMetPhi/F");
  fTree->Branch("tkmet"         ,&fTKMet         ,"fTKMet/F");
  fTree->Branch("tkmetphi"      ,&fTKMetPhi      ,"fTKMetPhi/F");
  fTree->Branch("mvamet"        ,&fMVAMet        ,"fMVAMet/F");
  fTree->Branch("mvametphi"     ,&fMVAMetPhi     ,"fMVAMetPhi/F");
  fTree->Branch("mvametu"       ,&fMVAMetUnity   ,"fMVAMet/F");
  fTree->Branch("mvametuphi"    ,&fMVAMetUnityPhi,"fMVAMetPhi/F");

  fTree->Branch("mttrue"         ,&fMtTrue        ,"fMtTrue/F");
  fTree->Branch("rawmttrue"      ,&fRawMtTrue     ,"fRawMtTrue/F");
  fTree->Branch("tkmttrue"       ,&fTKMtTrue      ,"fTKMtTrue/F");
  fTree->Branch("mvamttrue"      ,&fMVAMtTrue     ,"fMVAMtTrue/F");
  fTree->Branch("mvamtutrue"     ,&fMVAMtUnityTrue,"fMVAMtUnityTrue/F");

  fTree->Branch("metCov00"      ,&fMetCov00      ,"fMetCov00/F");
  fTree->Branch("metCov10"      ,&fMetCov10      ,"fMetCov10/F");
  fTree->Branch("metCov01"      ,&fMetCov01      ,"fMetCov01/F");
  fTree->Branch("metCov11"      ,&fMetCov11      ,"fMetCov11/F");
  fTree->Branch("metSig"        ,&fMetSig        ,"fMetSig/F");

  fTree->Branch("mvaCov00"      ,&fMVAMetCov00   ,"fMVAMetCov00/F");
  fTree->Branch("mvaCov10"      ,&fMVAMetCov10   ,"fMVAMetCov10/F");
  fTree->Branch("mvaCov01"      ,&fMVAMetCov01   ,"fMVAMetCov01/F");
  fTree->Branch("mvaCov11"      ,&fMVAMetCov11   ,"fMVAMetCov11/F");
  fTree->Branch("mvaMetSig"     ,&fMVAMetSig     ,"fMVAMetSig/F");

  fTree->Branch("mvaUCov00"     ,&fMVAMetUnityCov00   ,"fMVAMetCovUnity00/F");
  fTree->Branch("mvaUCov10"     ,&fMVAMetUnityCov10   ,"fMVAMetUnityCov10/F");
  fTree->Branch("mvaUCov01"     ,&fMVAMetUnityCov01   ,"fMVAMetUnityCov01/F");
  fTree->Branch("mvaUCov11"     ,&fMVAMetUnityCov11   ,"fMVAMetUnityCov11/F");
  fTree->Branch("mvaMetUSig"    ,&fMVAMetUnitySig     ,"fMVAMetUnitySig/F");
}
void EvtLoader::load(int iEvent) { 
  fVertices ->Clear();
  fEvtBr    ->GetEntry(iEvent);
  fVertexBr ->GetEntry(iEvent);
}
//MonoCentralPFJet80_PFMETnoMu
//HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v
void EvtLoader::addTrigger(std::string iName) { 
  fTrigString.push_back(iName);
}
bool EvtLoader::passFilter() { 
  return (fEvt->metFilterFailBits == 0);
}
bool EvtLoader::passTrigger() {
  bool lPass = false;
  for(unsigned int i0 = 0; i0 < fTrigString.size(); i0++) { 
    if(fTrigger->pass(fTrigString[i0],fEvt->triggerBits)) lPass = true;
  }
  return lPass;
}
bool EvtLoader::passTrigger(std::string iTrigger) { 
  return fTrigger->pass(iTrigger,fEvt->triggerBits);
}
bool EvtLoader::passSkim() { 
  //if(fMet > 200) std::cout << "==> Trig " << fITrigger << " -- " << (fITrigger % 2) << " -- " << fMet << " -- " << (fMetFilters % 2 == 0) << std::endl;
  //bool lTrig   = (fITrigger % 2 == 1);
  bool lMet    = fMet > 80;
  bool lFilter = fMetFilters % 2 == 0; 
  return (lMet && lFilter); 
}
unsigned int EvtLoader::triggerBit() {
  unsigned int lBit = 0;
  unsigned int lId = 0; 
  for(unsigned int i0 = 0; i0 < fTrigString.size(); i0++) { 
    if(fTrigger->pass(fTrigString[i0],fEvt->triggerBits))  lBit |= 1 << lId;
    lId++;
    if(i0 == 0) lId--;
  }
  return lBit;
}
void EvtLoader::fillEvent(std::vector<TLorentzVector> &iVecCorr) { 
  reset();
  //fSample     = iSample;
  fNPU        = fEvt->nPUmean;
  fNPUM       = fEvt->nPUmeanm;
  fNPUP       = fEvt->nPUmeanp;
  fNVtx       = nVtx();
  fITrigger   = triggerBit();
  fPUWeight   = puWeight(fNPU); 
  fMetFilters = metFilter(fEvt->metFilterFailBits);
  fRun        = fEvt->runNum;
  fLumi       = fEvt->lumiSec;
  fEvtV       = fEvt->evtNum;
  fRho        = fEvt->rhoJet;
  fRawSumEt   = -1;//fEvt->pfSumET;
  fMet        = fEvt->pfMETC;
  fMetPhi     = fEvt->pfMETCphi;
  fMetRaw     = fEvt->pfMET;
  fMetRawPhi  = fEvt->pfMETphi;
  fMetCov00   = fEvt->pfMETCov00;
  fMetCov01   = fEvt->pfMETCov01;
  fMetCov10   = fEvt->pfMETCov01;
  fMetCov11   = fEvt->pfMETCov11;
  fMetSig     = metSig(fMet,fMetPhi,fMetCov00,fMetCov01,fMetCov10,fMetCov11);

  fTKMet      = fEvt->trkMET;
  fTKMetPhi   = fEvt->trkMETphi;
  fMVAMet     = fEvt->mvaMET0;
  fMVAMetPhi  = fEvt->mvaMET0phi;
  fMVAMetCov00     = fEvt->mvaMET0Cov00;
  fMVAMetCov01     = fEvt->mvaMET0Cov01;
  fMVAMetCov10     = fEvt->mvaMET0Cov01;
  fMVAMetCov11     = fEvt->mvaMET0Cov11;
  fMVAMetSig       = metSig(fMVAMet,fMVAMetPhi,fMVAMetCov00,fMVAMetCov01,fMVAMetCov10,fMVAMetCov11);

  fMVAMetUnity          = fEvt->mvaMETU;
  fMVAMetUnityPhi       = fEvt->mvaMETUphi;
  fMVAMetUnityCov00     = fEvt->mvaMETUCov00;
  fMVAMetUnityCov01     = fEvt->mvaMETUCov01;
  fMVAMetUnityCov10     = fEvt->mvaMETUCov01;
  fMVAMetUnityCov11     = fEvt->mvaMETUCov11;
  fMVAMetUnitySig       = metSig(fMVAMetUnity,fMVAMetUnityPhi,fMVAMetUnityCov00,fMVAMetUnityCov01,fMVAMetUnityCov10,fMVAMetUnityCov11);

  TLorentzVector lCorr(0,0,0,0);
  for(unsigned int i0 =0; i0 < iVecCorr.size(); i0++) { 
    TLorentzVector pVec; pVec.SetPtEtaPhiM(iVecCorr[i0].Pt(),0,iVecCorr[i0].Phi(),0);
    lCorr += pVec;
  }
  fMtTrue                = mT(fMet,        fMetPhi,        lCorr);
  fRawMtTrue             = mT(fMetRaw,     fMetRawPhi,     lCorr);
  fTKMtTrue              = mT(fTKMet,      fTKMetPhi,      lCorr);
  fMVAMtTrue             = mT(fMVAMet,     fMVAMetPhi,     lCorr);
  fMVAMtUnityTrue        = mT(fMVAMetUnity,fMVAMetUnityPhi,lCorr);

  //Apply Corrections to the MET 
  if(iVecCorr.size() > 0) { 
    correctMet(fMet        ,fMetPhi,        lCorr);
    correctMet(fMetRaw     ,fMetRawPhi,     lCorr);
    correctMet(fTKMet      ,fTKMetPhi,      lCorr);
    correctMet(fMVAMet     ,fMVAMetPhi,     lCorr);
    correctMet(fMVAMetUnity,fMVAMetUnityPhi,lCorr);
    fMetSig               = metSig(fMet,fMetPhi,fMetCov00,fMetCov01,fMetCov10,fMetCov11);
    fMVAMetSig            = metSig(fMVAMet,fMVAMetPhi,fMVAMetCov00,fMVAMetCov01,fMVAMetCov10,fMVAMetCov11);
    fMVAMetUnitySig       = metSig(fMVAMetUnity,fMVAMetUnityPhi,fMVAMetUnityCov00,fMVAMetUnityCov01,fMVAMetUnityCov10,fMVAMetUnityCov11);
  }
  return;
}
void  EvtLoader::correctMet(float &iMet,float &iMetPhi,TLorentzVector &iCorr) { 
  TLorentzVector lVec;  lVec.SetPtEtaPhiM(iMet,0,iMetPhi,0);
  //Add back the correction vector
  lVec += iCorr;
  iMet    = lVec.Pt();
  iMetPhi = lVec.Phi();
}
float EvtLoader::metSig(float iMet,float iMetPhi,float iCov00,float iCov01,float iCov10,float iCov11) { 
  TMatrixD lInv(2,2);
  TLorentzVector lUVec; lUVec.SetPtEtaPhiM(iMet,0,iMetPhi,0);
  lInv(0,0) = iCov00; lInv(1,1) = iCov11; lInv(1,0) = iCov10; lInv(0,1) = iCov01;
  if(lInv.Determinant() != 0) lInv.Invert();
  double lSignificance = TMath::Sqrt(lUVec.Px()*lUVec.Px()*(lInv)(0,0) + 2.*lUVec.Px()*lUVec.Py()*(lInv)(1,0) + lUVec.Py()*lUVec.Py()*(lInv)(1,1));
  return lSignificance;
}
float EvtLoader::puWeight(double iPU) { 
  return fPUWeightHist->GetBinContent(fPUWeightHist->GetXaxis()->FindBin(iPU));
}
unsigned int EvtLoader::metFilter(unsigned int iMetFilter) { 
  unsigned int lWord = 0; 
  //Mapping of filters form Bacon to Bambu 
  lWord |=  (iMetFilter &  1)   << 0; //HBHE Noise
  lWord |=  (iMetFilter &  8)   << 1; //Ecal Dead Cell
  lWord |=  (iMetFilter & 16)   << 2; //Tracking Failure
  lWord |=  (iMetFilter & 32)   << 3; //EEBadSc
  lWord |=  (iMetFilter & 64)   << 4; //EcalLaser
  lWord |=  (iMetFilter & 128)  << 5; //tkManyStrip
  lWord |=  (iMetFilter & 256)  << 6; //tktooManyStrip
  lWord |=  (iMetFilter & 512)  << 7; //tkLogError
  lWord |=  (iMetFilter & 2)    << 8; //CSCTight
  //CSCLoose Halo Filter missing from our bitset
  lWord |=  (iMetFilter & 4)    << 10; //Hcal Laser
  return lWord;
}
unsigned int EvtLoader::nVtx() { 
  unsigned int lNVertex = 0; 
  for(int i0 = 0; i0 < fVertices->GetEntries(); i0++) { 
    TVertex *pVertex = (TVertex*) ((*fVertices)[i0]);
    if(fabs(pVertex->z) > 24) continue;
    if(pVertex->ndof    <  4) continue;
    float pX = pVertex->x;
    float pY = pVertex->y;
    float pRho  = sqrt(pX*pX+pY*pY);
    if(pRho             > 2.) continue;
    lNVertex++;
  }
  return lNVertex;
}
float  EvtLoader::mT(float &iMet,float &iMetPhi,TLorentzVector &iVec) { 
  float lDPhi = fabs(iMetPhi-iVec.Phi());
  if(fabs(lDPhi) > TMath::Pi()*2.-lDPhi) lDPhi = TMath::Pi()*2.-lDPhi;
  float lMt = sqrt(2.0*(iVec.Pt()*iMet*(1.0-cos(lDPhi))));
  return lMt;
}
