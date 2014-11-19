#include <utility>
#include "TF1.h"
#include "TRandom2.h"
#include "TLorentzVector.h"
#include "Gen.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

class LeptonsWW { 
public:
 typedef std::pair<TLorentzVector,int> lepton;
  LeptonsWW(TTree *iTree);
  ~LeptonsWW();
  float  smearLep(int iId, float &iPt);
  double effWeight(double iPt,double iEff=0.9,double iSigma=10.,double iMed=15.);
  void   fillWeight(int iId,double iPt);
  lepton convert(baconhep::TGenParticle *iPart);
  int    parent(baconhep::TGenParticle *iGen, TClonesArray *iGenParts);
  void   getAllLeps(std::vector<lepton> &iLeptons,TClonesArray *iArr,bool iFromBoson);
  void   reorderZ(std::vector<lepton> &lLeptons,std::vector<int> iOrder);
  void   setupSelector();
  void   fill(float &iPt,float &iPhi,float &iEta,float &iM,int &iId,lepton &iLepton);
  bool   findRealLep(Gen &iGen);

private:
  TRandom2* fRandom;
  TF1*      fFakeRate;
  float     fEffWeight;
};
