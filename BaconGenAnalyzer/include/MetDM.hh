#include "TLorentzVector.h"
#include "RecoilCorrector.hh"
#include "GenJet.hh"
#include "Jets.hh"

class MetDM{ 
public:
  MetDM(TTree *iTree);
  ~MetDM();
  void setupMetTree(TTree *iTree);
  void clearMet();
  void subtractGenJets(TLorentzVector &iRecoil,Jets &iJets,bool iSubtract=true);  
  void subtractJets   (TLorentzVector &iRecoil,Jets &iJets,bool iSubtract=true);
  void fillMet(GenJet &iGen, Jets &iJets,bool iFake=false);
  
private:
  float fMet;
  float fMetPhi;
  float fMetCov00;
  float fMetCov01;
  float fMetCov10;
  float fMetCov11;
  float fMt1;
  float fMt2;
  
  RecoilCorrector *fCorrector;
};
