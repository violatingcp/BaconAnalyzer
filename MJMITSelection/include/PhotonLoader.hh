#include "TTree.h"
#include "TLorentzVector.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TPhoton.hh"

using namespace baconhep;

class PhotonLoader { 
public:
  PhotonLoader(TTree *iTree);
  ~PhotonLoader();
  void reset();
  void setupTree(TTree *iTree);
  void load (int iEvent);
  bool selectPhotons(std::vector<TLorentzVector> &iVetoes,float iRho);
  bool selectPhoton (std::vector<TLorentzVector> &iVetoes,float iRho);
  void fillVetoes(std::vector<TLorentzVector> &iVec);  
  //Vetos
  bool vetoPhoton(float iRho);
  //Ids
  double effArea(double eta,int iType);
  bool passMedium      (TPhoton *iPhoton,float iRho);
  bool passCiCPhotonPre(TPhoton *p,float iRho);  
  bool passCiCPFIso    (TPhoton *p,float iRho);
protected: 
  TClonesArray *fPhotons;
  TBranch      *fPhotonBr;
  TTree        *fTree;
  int             fNPhotons;
  TLorentzVector *fPtr1;
};
