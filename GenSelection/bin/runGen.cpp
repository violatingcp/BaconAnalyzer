#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TRandom2.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

//Flat Ntuple
TTree *fTree;
TRandom *fRandom = new TRandom2(0xDEADBEEF);
#include "../include/Jets.hh"
#include "../include/GenSelect.hh"
#include "../include/GenWeight.hh"
Jets       *fGenJet       = 0;
GenSelect  *fGen          = 0;
GenWeight  *fGenW         = 0;

void clear(int i0) {
  fGen   ->clearGen (i0);
  fGenJet->clearJets(i0);
  fGenW ->clearGen (i0);
}
int select() { 
  int lId = -1;
  return lId;
}
TTree* load(std::string iName) { 
  TFile *lFile = TFile::Open(iName.c_str());
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  return lTree;
}
bool fill(TFile *iFile,float iXS,float iTot) { 
  std::vector<TLorentzVector> lVec;
  fGen    ->fillGen   (double(iTot),lVec);
  fGenJet ->fillJets  (lVec);
  fGen->fMCWeight  = iXS/iTot;
  fGen->fXS        = iXS;
  fGenW->fill();
  iFile->cd();
  fTree->Fill();
  return true;
}
int main(int argc, char **argv) { 
  //gROOT->ProcessLine("#include <vector>");          
  std::string iName    = argv[2];
  double      iNEvents = atof(argv[1]);
  double      iXS      = atof(argv[3]);

  TTree *lTree   = load(iName);
  TFile *lFile   = new TFile("Output.root","RECREATE");
  fTree     = new TTree("Events","Events");
  fGen      = new GenSelect(lTree,fTree);
  fGenJet   = new Jets     (lTree,fTree);
  fGenW     = new GenWeight(lTree,fTree);
  double lTot = lTree->GetEntries();
  if(iNEvents < 0) iNEvents = lTot;
  fGen->fMCWeight = iXS/double(iNEvents);
  lFile->cd();
  for(int i0 = 0; i0 < lTot; i0++) { 
    if(i0 % 1000 == 0) std::cout  << " ----> " << i0/lTot << " -- " << fGen->fMCWeight << std::endl;
    clear(i0);
    fill(lFile,iXS,iNEvents);
  }
  lFile->cd();
  fTree->Write();
  lFile->Close();
}
