#ifndef JetRes_hh
#define JetRes_hh
#include <vector>
#include <string>
#include "RecoilCorrector.hh"
#include "TF1.h"
#include "TGraph.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "BaconAna/DataFormats/interface/TGenJet.hh"

class JetRes { 
public:
  struct TRes { 
    TGraph *func;
    double etaMin;
    double etaMax;
    TRes(TGraph *iF,double iEtaMin,double iEtaMax) { 
      func   = iF;
      etaMin = iEtaMin;
      etaMax = iEtaMax;
    }
  };
  JetRes(bool iData=true);
  ~JetRes() {}
  std::vector<TLorentzVector> smearEvent(double iGMet,double iGMetPhi,double iHT,std::vector<baconhep::TGenJet*> &iJets);
  double* getId(std::string iName);
  double res(double iCSV,double iPt,double iEta,double iSig=0.);
  double smear(double iCSV,double iPt,double iEta,double iSig=0.);
  double findSig(std::vector<TRes*> &iBase,std::vector<TRes*> &iSig1,std::vector<TRes*> &iSig2,double iPt,double iEta);
  TRes*  find(std::vector<TRes*> &iBase,double iEta);
  TLorentzVector smearedJet(baconhep::TGenJet *iJet,double iCSV,double iSig);
protected:
  std::vector<TRes*> fBase;
  std::vector<TRes*> fHBase;
  std::vector<TRes*> fRec;
  std::vector<TRes*> fS1Base;
  std::vector<TRes*> fS1HBase;
  std::vector<TRes*> fS1Rec;
  std::vector<TRes*> fS2Base;
  std::vector<TRes*> fS2HBase;
  std::vector<TRes*> fS2Rec;
  TRandom * fRandom;
  bool      fData;
  RecoilCorrector *fCorrector;
};
#endif
