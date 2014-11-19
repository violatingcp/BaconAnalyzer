#ifndef GENJET_HH
#define GENJET_HH
#include <vector>
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TGenJet.hh"

class GenJet { 
public:
  GenJet(TTree *iTree,TTree *iOTree);
  ~GenJet();
  void loadGen     (TTree *iTree);
  void setupGenTree(TTree *iTree);
  void clearGen(int i0);
  void replace(int iId,baconhep::TGenParticle *iGen,std::vector<baconhep::TGenParticle*> &iBoson,bool iRemove);
  void findBosons(std::vector<baconhep::TGenParticle*> & iBoson,int iIdMin,int iIdMax);
  void fillGen();

  TClonesArray *fGenParts;
  TClonesArray *fGenJets;

  //Decay Mode
  int   fDecay; 
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

  float fJ1Pt;
  float fJ1Eta;
  float fJ1Phi;
  float fJ1M;
  int   fJ1Id;

private:
  TBranch      *fGenBr;
  TBranch      *fGenPartBr;
  TBranch      *fGenJetBr;
  baconhep::TGenEventInfo     *fGen;
};
#endif
