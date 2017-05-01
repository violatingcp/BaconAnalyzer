#ifndef GENWEIGHT_HH
#define GENWEIGHT_HH
#include <vector>
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TLHEWeight.hh"

class GenWeight { 
public:
  GenWeight(TTree *iTree,TTree *iOTree);
  ~GenWeight();
  void loadGen     (TTree *iTree);
  void setupGenTree(TTree *iTree);
  void clearGen(int i0);
  void fill();
  
  float fScale10;
  float fScale12;
  float fScale20;
  float fScale21;
  float fScale22;
  float fScale00;
  float fScale01;
  float fScale02;
  float fPDFUp  ;
  float fPDFDown;
  float fPDF[100];

private:
  TBranch      *fGenWeightBr;
  TClonesArray *fGenWeight;
  bool          fNoWeight;
};
#endif
