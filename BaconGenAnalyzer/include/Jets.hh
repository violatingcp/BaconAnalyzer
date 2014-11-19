#ifndef JETS_HH
#define JETS_HH
#include "TClonesArray.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom2.h"
#include "Gen.hh"
#include "GenJet.hh"

class Jets{ 
public:
  Jets(TTree *iTree);
  ~Jets();
  
  void setupJetTree(TTree *iTree);
  void clearJets();
  void buildLeptonVec(Gen &iGen);
  void selectJets(Gen    &iGen);
  void selectJets(GenJet &iGen);
  void computeFinal();
  void setSmear(int iBin);
  void smear(float &iJPt,float &iJEta,float &iJPhi,float &iJM);
  void smear();
  void fillJets(Gen    &iGen);
  void fillJets(GenJet &iGen);
  void insert(baconhep::TGenJet *iJet,std::vector<baconhep::TGenJet*> &lJets);
  std::vector<TLorentzVector> fLeptons;
  
  float fJetPt1;
  float fJetEta1;
  float fJetPhi1;
  float fJetM1;
  int   fJetId1;
  
  float fJetPt2;
  float fJetEta2;
  float fJetPhi2;
  float fJetM2;
  int   fJetId2;

  float fJCPt;
  float fJCEta;
  float fJCPhi;
  float fJCM;
  int   fJCId;
  
  float fGJetPt1;
  float fGJetEta1;
  float fGJetPhi1;
  float fGJetM1;
  float fGJetPt2;
  float fGJetEta2;
  float fGJetPhi2;
  float fGJetM2;
  float fGJCPt;
  float fGJCEta;
  float fGJCPhi;
  float fGJCM;
  float fNJets;
private: 
  
  float fMJJ;
  float fJDEta;
  float fJDPhi;
  
  TF1 * fPtSmear;
  TF1 * fPhiSmear;
  TF1 * fEtaSmear;
  double *fPtPar0 ,*fPtPar1 ,*fPtPar2 ,*fPtPar3;
  double *fEtaPar0,*fEtaPar1,*fEtaPar2,*fEtaPar3,*fEtaPar4;
  double *fPhiPar0,*fPhiPar1,*fPhiPar2,*fPhiPar3,*fPhiPar4;
  
  TRandom2 *fRandom;
  
};
#endif
