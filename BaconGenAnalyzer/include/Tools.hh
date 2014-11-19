#ifndef TUPLETOOLS_HH
#define TUPLETOOLS_HH

#include <string>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"

class Tools{
public:
  Tools();
  ~Tools();
  static  double deltaPhi(double iPhi1,double iPhi2);
  static  double deltaR  (double iEta1,double iPhi1,double iEta2,double iPhi2) ;
  static  TTree* load    (std::string iName,TFile *iFile,int &iNMCEvents);
  static  double mT      (double iE1,double iPhi1,double iE2,double iPhi2);
  static  double pZeta   (int iId,double iPt1,double iEta1,double iPhi1,double iPt2,double iEta2,double iPhi2,double iMet,double iMetPhi,double iM2=0);
  static  double pVis    (int iId,double iPt1,double iEta1,double iPhi1,double iPt2,double iEta2,double iPhi2,double iM2=0);
  static  double pMet    (int iId,double iPt1,double iEta1,double iPhi1,double iPt2,double iEta2,double iPhi2,double iMet,double iMetPhi,double iM2=0); 
  static  TLorentzVector system(int iId,double iPt1,double iEta1,double iPhi1,double iPt2,double iEta2,double iPhi2,double iMet,double iMetPhi,double iM2 = 0);
  static  double collinear(int iId,double iPt1,double iEta1,double iPhi1,double iPt2,double iEta2,double iPhi2,double iMet,double iMetPhi,double iM2 = 0);
  static  TLorentzVector getVisMass(int iId,double iPt1,double iEta1,double iPhi1,double iPt2,double iEta2,double iPhi2,double iM2 = 0);
  static  double u1(double iPt,double iPhi,double iMet,double iMPhi,double iTruePhi);
  static  double u2(double iPt,double iPhi,double iMet,double iMPhi,double iTruePhi);
  static  double uphi(double iPt,double iPhi,double iMet,double iMPhi);
};

#endif
