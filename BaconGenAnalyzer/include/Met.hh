#include "TLorentzVector.h"
#include "RecoilCorrector.hh"
#include "NTFit.hh"
#include "Gen.hh"
#include "Jets.hh"

class Met{ 
public:
  Met(TTree *iTree);
  ~Met();
  void setupMetTree(TTree *iTree);
  void clearMet();
  void subtractGenJets(TLorentzVector &iRecoil,Jets &iJets,bool iSubtract=true);  
  void subtractJets   (TLorentzVector &iRecoil,Jets &iJets,bool iSubtract=true);
  void fillMet(Gen &iGen, Jets &iJets,bool iFake=false);
  
private:
  float fMet;
  float fMetPhi;
  float fMetCov00;
  float fMetCov01;
  float fMetCov10;
  float fMetCov11;
  float fMt1;
  float fMt2;
  float fPZetaVis;
  float fPZetaMiss;
  float fPZeta;
  float fPtH;

  //System variables
  float fMNT   ;
  float fPtNT  ;
  float fPhiNT ;
  float fYNT   ;

  float fMVis  ;
  float fPhiVis;
  float fMR;
  
  RecoilCorrector *fCorrector;
  NTFit           *fNTFitWW;
  NTFit           *fNTFitWZ;
};
