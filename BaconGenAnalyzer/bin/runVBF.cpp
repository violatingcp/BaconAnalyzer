#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TRandom2.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "../include/Tools.hh"

//Flat Ntuple
TTree *fTree;
TRandom *fRandom = new TRandom2(0xDEADBEEF);

#include "../include/Gen.hh"
#include "../include/Jets.hh"
#include "../include/Met.hh"
#include "../include/LeptonsWW.hh"

Gen  *fGen          = 0;
Jets *fJets         = 0;  
Met  *fMet          = 0;
LeptonsWW *fLeptons = 0;

void clear(int i0) {
  fGen ->clearGen (i0);
  fMet ->clearMet ();
  fJets->clearJets();
}
/*
bool passKineCuts() {
  //if(fabs(fId1)  == 17 && fabs(fEta1) > 2.1) return false;
  if(fDecay == 0 && (fPt1 < 45  || fPt2 < 45) ) return false; //CMS Default => Cerntainly can be switched
  if(fDecay == 1 && (fPt1 < 20  || fPt2 < 20) ) return false; //Currently use 20 30 for mu Tau and 24 30 for eTau ==> public values are 20 20 and 24 20
  if(fDecay == 2 && (fPt1 < 20  || fPt2 < 10) ) return false; //ee and mm are super tough, but can be made at little bit less sensitive than emu
  return true;
}
bool passMetBasedCuts() { 
  if(fDecay == 0  && fPtH   < 100    ) return false; // Higgs PT is the key cut for tau tau
  if(fDecay == 1  && fMt1   > 2000    ) return false; // Cut to remove W+Jets
  //if(fDecay == 2  && fPZeta > -2500  ) return false; // Cut to remove the ttbar
  if(isnan(fPt2) || isinf(fPt2)   ) return false;
  return true;
}
*/
int select() { 
  fLeptons->findRealLep(*fGen);
  int lId = -1;
  if(fGen->fPt1 > 0) lId = 0;
  if(fGen->fPt2 > 0) lId = 1;
  if(fGen->fPt3 > 0) lId = 2;
  if(fGen->fPt3 > 0) lId = 3;
  //ApplyMetCuts(lId);
  return lId;
}
bool fill(TFile *iFile) { 
  fGen ->fillGen   ();
  fJets->fillJets  (*fGen);
  fMet ->fillMet   (*fGen,*fJets);
  //if(!passMetBasedCuts()) return false;
  iFile->cd();
  fTree->Fill();
  return true;
}
int main(int argc, char **argv) { 
  //gROOT->ProcessLine("#include <vector>");          
  std::string iName  = argv[1];//"test.root";
  double      iXS    = atof(argv[2]);//1.0;
  double      iNFile = atof(argv[3]);

  int lNMCEvents = 0; 
  TFile *lIFile  = 0;
  TTree *lTree   = Tools::load(iName,lIFile,lNMCEvents); 
  TFile *lFile   = new TFile("Output.root","RECREATE");
  fTree     = new TTree("Events","Events");
  fGen      = new Gen(lTree,fTree);
  fLeptons  = new LeptonsWW(fTree);
  fJets     = new Jets     (fTree);
  fMet      = new Met      (fTree);
  
  double lTot = lTree->GetEntries();
  cout << "===> Events " << lNMCEvents << endl;
  fGen->fMCWeight = iXS/double(iNFile);
  lFile->cd();
  for(int i0 = 0; i0 < lTot; i0++) { 
    if(i0 % 1000 == 0) std::cout  << " ----> " << i0/lTot << " -- " << fGen->fMCWeight << std::endl;
    clear(i0);
    if(!select()) continue;
    fill(lFile);
  }
  lFile->cd();
  fTree->Write();
  lFile->Close();
}
