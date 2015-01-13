#ifndef GEN_HH
#define GEN_HH
#include <vector>
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TGenJet.hh"

class Gen { 
public:
  Gen(TTree *iTree,TTree *iOTree);
  ~Gen();
  void loadGen     (TTree *iTree);
  void setupGenTree(TTree *iTree);
  void clearGen(int i0);
  void replace(int iId,baconhep::TGenParticle *iGen,std::vector<baconhep::TGenParticle*> &iBoson,bool iRemove);
  void findBosons(std::vector<baconhep::TGenParticle*> & iBoson);  
  void findBosonsPythia6(std::vector<baconhep::TGenParticle*> & iBoson);  
  void fillGen(bool iPythia6=false);

  TClonesArray *fGenParts;
  TClonesArray *fGenJets;

  //Decay Mode
  int   fDecay;  //0=>tt  1=>lt  2=>ll
  float fMCWeight;
  float fEffWeight;
  
  //Flat ntuple variables
  float fVPt;
  float fVEta;
  float fVPhi;
  float fVM;

  float fV1Pt;
  float fV1Eta;
  float fV1Phi;
  float fV1M;
  int   fV1Id;

  float fV2Pt;
  float fV2Eta;
  float fV2Phi;
  float fV2M;
  int   fV2Id;
  
  float fPt1;
  float fEta1;
  float fPhi1;
  float fM1;
  float fCosPhi1;
  int   fId1;
  
  float fPt2;
  float fEta2;
  float fPhi2;
  float fM2;
  float fCosPhi2;
  int   fId2;

  float fPt3;
  float fEta3;
  float fPhi3;
  float fM3;
  int   fId3;

  float fPt4;
  float fEta4;
  float fPhi4;
  float fM4;
  int   fId4;

private:
  TBranch      *fGenBr;
  TBranch      *fGenPartBr;
  TBranch      *fGenJetBr;
  baconhep::TGenEventInfo     *fGen;
};
#endif
