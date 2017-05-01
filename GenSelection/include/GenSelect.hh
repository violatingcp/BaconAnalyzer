#ifndef GENSELECT_HH
#define GENSELECT_HH
#include <vector>
#include "TLorentzVector.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"

class GenSelect { 
public:
  GenSelect(TTree *iTree,TTree *iOTree);
  ~GenSelect();
  void loadGen     (TTree *iTree);
  void setupGenTree(TTree *iTree);
  void clearGen(int i0);
  void replace(int iId,baconhep::TGenParticle *iGen,std::vector<baconhep::TGenParticle*> &iBoson,bool iRemove);
  void add(baconhep::TGenParticle *iGen,std::vector<baconhep::TGenParticle*> &iBoson);
  void findBosons(std::vector<baconhep::TGenParticle*> & iBoson,int iIdMin,int iIdMax,int iStatusMin,int iStatusMax,bool iLast=true);
  void npartons();
  float iso(double iEta,double iPhi,int iIdMin,int iIdMax,double iStatusMin,double iStatusMax,double iDRMax);
  void fillGen(double iNum,std::vector<TLorentzVector> &iVec);
  void getTop(std::vector<baconhep::TGenParticle*> &iParts);
  void met();
  TClonesArray *fGenParts;

  //Decay Mode
  int   fNPartons;
  int   fDecay;  //0=>tt  1=>lt  2=>ll
  float fQ2;
  float fEvtWeight;
  float fMCWeight;
  float fMCWeight2;
  float fXS;
  float fXS2;
  
  //Flat ntuple variables
  float fMPt;
  float fMEta;
  float fMPhi;
  float fMM;
  int   fMId;
  int   fMStatus;

  float fM2Pt;
  float fM2Eta;
  float fM2Phi;
  float fM2M;
  int   fM2Id;
  int   fM2Status;

  float fVPt;
  float fVEta;
  float fVPhi;
  float fVM;
  int   fVId;
  int   fVStatus;
  float fVIso;
  float fDVIso;

  float fVMEPt;
  float fVMEEta;
  float fVMEPhi;
  float fVMEM;
  int   fVMEId;
  int   fVMEStatus;
  float fVMEIso;
  float fDVMEIso;
  
  float fPt1;
  float fEta1;
  float fPhi1;
  float fM1;
  float fCosPhi1;
  int   fId1;
  int   fStatusId1;
  
  float fPt2;
  float fEta2;
  float fPhi2;
  float fM2;
  float fCosPhi2;
  int   fId2;
  int   fStatusId2;

  float fTPt;
  float fTEta;
  float fTPhi;
  float fTM;

  float fT1Pt;
  float fT1Eta;
  float fT1Phi;
  float fT1M;

  float fMet;
  float fMetPhi;

private:
  TBranch      *fGenBr;
  TBranch      *fGenPartBr;
  baconhep::TGenEventInfo     *fGen;
};
#endif
