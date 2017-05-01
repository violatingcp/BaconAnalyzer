#ifndef SMEAR_HH
#define SMEAR_HH
#include "TClonesArray.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom2.h"
#include "../include/JetRes.hh"
#include "BaconAna/DataFormats/interface/TGenJet.hh"

class Smear{ 
public:
  Smear(TTree *iTree,TTree *iOTree);
  ~Smear();
  
  void setupTree(TTree *iTree);
  void clear(int i0);
  void smear(double iMet,double iMetPhi,double iHT,std::vector<baconhep::TGenJet*> &iJets);
  JetRes *fJetRes;
  float fGJetPt1;
  float fGJetEta1;
  float fGJetPhi1;
  float fGJetM1;
  float fGJetId1;
  float fGJetPt2;
  float fGJetEta2;
  float fGJetPhi2;
  float fGJetM2;
  float fGJetId2;
  float fGJCPt;
  float fGJCEta;
  float fGJCPhi;
  float fGJCM;
  float fGJCId;
  float fHT;
  float fMet;
  float fMetPhi;
};
#endif
