#ifndef TOPJETS_HH
#define TOPJETS_HH
#include "TClonesArray.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom2.h"
#include "BaconAna/DataFormats/interface/TGenJet.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

class TopJets{ 
public:
  TopJets(TTree *iTree,TTree *iOTree);
  ~TopJets();
  
  void setupJetTree(TTree *iTree);
  void clearJets(int i0);
  void selectJets (std::vector<TLorentzVector> &iVec,std::vector<baconhep::TGenParticle*>  &iParts);  
  void selectBJets(std::vector<TLorentzVector> &iVec,std::vector<baconhep::TGenParticle*>  &iParts);  
  void computeFinal();
  void fillJets(std::vector<TLorentzVector> &iVec,std::vector<baconhep::TGenParticle*>  &iParts);
  void insert(baconhep::TGenJet *iJet,std::vector<baconhep::TGenJet*> &lJets);
  int  matchTop(double iEta,double iPhi,std::vector<baconhep::TGenParticle*>  &iParts);

  float fGJetPt1;
  float fGJetEta1;
  float fGJetPhi1;
  float fGJetM1;
  float fGJetId1;
  float fTGJetId1;
  float fGJetPt2;
  float fGJetEta2;
  float fGJetPhi2;
  float fGJetM2;
  float fGJetId2;
  float fTGJetId2;
  float fGJCPt;
  float fGJCEta;
  float fGJCPhi;
  float fGJCM;
  float fGJCId;
  float fTGJCId;

  float fGFJPt;
  float fGFJEta;
  float fGFJPhi;
  float fGFJM;
  float fGFJMSD;
  float fGFJT2T1;
  float fGFJT3T2;
  int   fGFJId;
  float fNJets;
  float fNJets50;

  float fBGJetPt1;
  float fBGJetEta1;
  float fBGJetPhi1;
  float fBGJetM1;
  float fBGJetId1;
  float fBTGJetId1;
  float fBGJetPt2;
  float fBGJetEta2;
  float fBGJetPhi2;
  float fBGJetM2;
  float fBGJetId2;
  float fBTGJetId2;
  float fBGJCPt;
  float fBGJCEta;
  float fBGJCPhi;
  float fBGJCM;
  float fBGJCId;
  float fBTGJCId;
  float fNBJets;
  float fNBJets50;
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
