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
  void findBosons(std::vector<baconhep::TGenParticle*> & iBoson,int iIdMin,int iIdMax,int iStatusMin,int iStatusMax);
  void npartons();
  void fillGen(double iNum,std::vector<TLorentzVector> &iVec);
  TClonesArray *fGenParts;

  //Decay Mode
  int   fNPartons;
  int   fDecay;  //0=>tt  1=>lt  2=>ll
  float fEvtWeight;
  float fMCWeight;
  float fMCWeight2;
  float fXS;
  float fXS2;
  
  //Flat ntuple variables
  float fVPt;
  float fVEta;
  float fVPhi;
  float fVM;
  int   fVId;
  
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

private:
  TBranch      *fGenBr;
  TBranch      *fGenPartBr;
  baconhep::TGenEventInfo     *fGen;
};
#endif
