#ifndef JETS_HH
#define JETS_HH
#include "TClonesArray.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom2.h"
#include "BaconAna/DataFormats/interface/TGenJet.hh"

class Jets{ 
public:
  Jets(TTree *iTree,TTree *iOTree);
  ~Jets();
  
  void setupJetTree(TTree *iTree);
  void clearJets(int i0);
  void selectJets(std::vector<TLorentzVector> &iVec);
  void computeFinal();
  void fillJets(std::vector<TLorentzVector> &iVec);
  void insert(baconhep::TGenJet *iJet,std::vector<baconhep::TGenJet*> &lJets);
  
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

  float fGFJPt;
  float fGFJEta;
  float fGFJPhi;
  float fGFJM;
  float fGFJMTrim;
  float fFJT2T1;
  int   fGFJId;
  float fNJets;
  float fNJets50;

private: 
  TBranch      *fGenJetBr;
  TClonesArray *fGenJets;
  TBranch      *fGenFatJetBr;
  TClonesArray *fGenFatJets;
  float fMJJ;
  float fPtJJ;
  float fEtaJJ;
  float fPhiJJ;
  float fJDEta;
  float fJDPhi;  
};
#endif
