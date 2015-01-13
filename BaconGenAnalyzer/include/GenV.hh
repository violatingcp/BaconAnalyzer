#ifndef GENV_HH
#define GENV_HH
#include <vector>
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TGenJet.hh"

class GenV { 
public:
  GenV(TTree *iTree,TTree *iOTree);
  ~GenV();
  void loadGen     (TTree *iTree);
  void setupGenTree(TTree *iTree);
  void clearGen(int i0);
  void replace(int iId,baconhep::TGenParticle *iGen,std::vector<baconhep::TGenParticle*> &iBoson,bool iRemove);
  void findBosons(std::vector<baconhep::TGenParticle*> & iBoson,int iIdMin,int iIdMax,int iStatus);
  void findDaughters(baconhep::TGenParticle *iParticle,std::vector<baconhep::TGenParticle*> & iDaughter,int iIdMin,int iIdMax);
  void fillGen(int iPythia6=0);

  TClonesArray *fGenParts;
  TClonesArray *fGenJets;
  TClonesArray *fGenFatJets;

  //Decay Mode
  int   fDecay;  //0=>tt  1=>lt  2=>ll
  float fMCWeight;
  float fEffWeight;
  
  //Flat ntuple variables
  float fDPt;
  float fDEta;
  float fDPhi;
  float fDM;
  int   fDId;

  //Reference for the MET
  float fV1Pt;
  float fV1Eta;
  float fV1Phi;
  float fV1M;
  int   fV1Id;

  float fD1Pt;
  float fD1Eta;
  float fD1Phi;
  float fD1M;
  int   fD1Id;

  float fD2Pt;
  float fD2Eta;
  float fD2Phi;
  float fD2M;
  int   fD2Id;

  float fQ1Pt;
  float fQ1Eta;
  float fQ1Phi;
  float fQ1M;
  int   fQ1Id;

  float fQ2Pt;
  float fQ2Eta;
  float fQ2Phi;
  float fQ2M;
  int   fQ2Id;

  float fBJPt;
  float fBJEta;
  float fBJPhi;
  float fBJM;
  int   fBJId;

  float fBPt;
  float fBEta;
  float fBPhi;
  float fBM;
  int   fBId;
private:
  TBranch      *fGenBr;
  TBranch      *fGenPartBr;
  TBranch      *fGenJetBr;
  TBranch      *fGenFatJetBr;
  baconhep::TGenEventInfo     *fGen;
};
#endif
