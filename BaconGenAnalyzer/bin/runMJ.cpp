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

#include "../include/GenJet.hh"
#include "../include/Jets.hh"
#include "../include/MetDM.hh"

GenJet  *fGen          = 0;
Jets    *fJets         = 0;  
MetDM   *fMet          = 0;

void clear(int i0) {
  fGen ->clearGen (i0);
  fMet ->clearMet ();
  fJets->clearJets();
}
int select() { 
  int lId = -1;
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
  std::string iName  = argv[1];
  double      iXS    = atof(argv[2]);//1.0;
  double      iNFile = atof(argv[3]);

  int lNMCEvents = 0; 
  TFile *lIFile  = 0;
  TTree *lTree   = Tools::load(iName,lIFile,lNMCEvents); 
  TFile *lFile   = new TFile("Output.root","RECREATE");
  fTree     = new TTree("Events","Events");
  fGen      = new GenJet(lTree,fTree);
  fJets     = new Jets        (fTree);
  fMet      = new MetDM       (fTree);
  
  double lTot = lTree->GetEntries();
  cout << "===> Events " << lNMCEvents << endl;
  fGen->fMCWeight = iXS/double(iNFile);
  lFile->cd();
  for(int i0 = 0; i0 < lTot; i0++) { 
    if(i0 % 1000 == 0) std::cout  << " ----> " << i0/lTot << " -- " << fGen->fMCWeight << std::endl;
    clear(i0);
    //if(!select()) continue;
    fill(lFile);
  }
  lFile->cd();
  fTree->Write();
  lFile->Close();
}
