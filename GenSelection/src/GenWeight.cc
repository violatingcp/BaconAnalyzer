#include "../include/GenWeight.hh"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iostream>

GenWeight::GenWeight(TTree *iTree,TTree *iOTree) { 
  loadGen(iTree);
  setupGenTree(iOTree);
}
GenWeight::~GenWeight() { }
void GenWeight::loadGen(TTree *iTree) { 
  baconhep::TLHEWeight   ::Class() ->IgnoreTObjectStreamer();
  fGenWeight        = new TClonesArray("baconhep::TLHEWeight");
  iTree->SetBranchAddress("LHEWeight"    ,&fGenWeight);   fGenWeightBr   = iTree->GetBranch("LHEWeight"  ); 
}
void GenWeight::setupGenTree(TTree *iTree) {
  iTree->Branch("scale12",&fScale12,"fScale12/F");
  iTree->Branch("scale10",&fScale10,"fScale10/F");
  iTree->Branch("scale21",&fScale21,"fScale21/F");
  iTree->Branch("scale22",&fScale22,"fScale22/F");
  iTree->Branch("scale20",&fScale20,"fScale20/F");
  iTree->Branch("scale01",&fScale01,"fScale01/F");
  iTree->Branch("scale02",&fScale02,"fScale02/F");
  iTree->Branch("scale00",&fScale00,"fScale00/F");
  iTree->Branch("pdf"    ,&fPDF    ,"fPDF[100]/F");
  iTree->Branch("pdfUp"  ,&fPDFUp  ,"fPDFUp/F");
  iTree->Branch("pdfDown",&fPDFDown,"fPDFDown/F");
}
void GenWeight::clearGen(int i0) {
  fGenWeightBr->GetEntry(i0);
  fScale10 = 1.;
  fScale12 = 1.;
  fScale20 = 1.;
  fScale21 = 1.;
  fScale22 = 1.;
  fScale00 = 1.;
  fScale01 = 1.;
  fScale02 = 1.;
  fPDFUp   = 0.; 
  fPDFDown = 0.;
  for(int i0 = 0; i0 < 100; i0++) fPDF[i0] = 0.;
}
void GenWeight::fill() { 
  double lBaseWeight = 1.;
  int lNUp   = 0; 
  int lNDown = 0; 
  for(int i0 = 0; i0 < fGenWeight->GetEntriesFast(); i0++) { 
    baconhep::TLHEWeight *pWeight = (baconhep::TLHEWeight*) fGenWeight->At(i0);
    if(i0 == 0) lBaseWeight = pWeight->weight;
    if(i0 == 1) fScale10 = pWeight->weight/lBaseWeight;
    if(i0 == 2) fScale12 = pWeight->weight/lBaseWeight;
    if(i0 == 3) fScale20 = pWeight->weight/lBaseWeight;
    if(i0 == 4) fScale21 = pWeight->weight/lBaseWeight;
    if(i0 == 5) fScale22 = pWeight->weight/lBaseWeight;
    if(i0 == 6) fScale00 = pWeight->weight/lBaseWeight;
    if(i0 == 7) fScale01 = pWeight->weight/lBaseWeight;
    if(i0 == 8) fScale02 = pWeight->weight/lBaseWeight;
    if(i0 >  8 && i0 < 109) fPDF[i0-9] = pWeight->weight/lBaseWeight;
    if(pWeight->weight/lBaseWeight > 1.) fPDFUp   += (pWeight->weight/lBaseWeight-1.)*(pWeight->weight/lBaseWeight-1.);
    if(pWeight->weight/lBaseWeight < 1.) fPDFDown += (pWeight->weight/lBaseWeight-1.)*(pWeight->weight/lBaseWeight-1.);
    if(pWeight->weight/lBaseWeight > 1.) lNUp++;
    if(pWeight->weight/lBaseWeight < 1.) lNDown ++;
  }
  fPDFUp   = 1.+sqrt(fPDFUp)  /lNUp;
  fPDFDown = 1.-sqrt(fPDFDown)/lNDown;
}
