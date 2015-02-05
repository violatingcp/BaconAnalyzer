#include "TTree.h"
#include "TLorentzVector.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "BaconAna/DataFormats/interface/TJet.hh"
//#include "BaconAna/DataFormats/interface/TAddJet.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"

using namespace baconhep;

class JetLoader { 
public:
  JetLoader(TTree *iTree,bool iData=false,std::string iHLTFile="/afs/cern.ch/user/p/pharris/pharris/public/bacon/CMSSW_5_3_13/src/BaconAna/DataFormats/data/HLTFile_v0");
  ~JetLoader();
  void reset();
  void reset(TJet &iJet);

  void setupTree(TTree *iTree);
  void load (int iEvent);
  
  bool selectJets(std::vector<TLorentzVector> &iVetoes,double iRho);
  void fillVars(TJet *iJet,double iRho,TLorentzVector *iPtr,TJet &iSaveJet);
  //Selectors
  bool vetoJet();
  bool passLoose      (TJet *iJet);
  bool passTight      (TJet *iJet);
  bool passVeto       (TJet *iJet);
  bool passPUId       (TJet *iJet);
  TJet    *fatJet     (TJet *iJet);
  float    mva        (TJet *iJet);
  //Trigger Stuff
  void addTrigger (std::string iName);
  bool passTrigObj(TJet *iJet,int iId);
protected: 
  TClonesArray *fJets;
  TBranch      *fJetBr;

  TClonesArray *fFatJets;
  TBranch      *fFatJetBr;
  
  TTree        *fTree;
  std::vector<std::string> fTrigString;
  TTrigger     *fTrigger;

  unsigned int fNJets;
  unsigned int fNBTags;
  unsigned int fNBTags10;
  unsigned int fNQTags;
  float        fJDPhi;
  float        fJDEta;
  float        fJDPullY;
  float        fJDPullPhi;

  TLorentzVector *fPtr1;
  TLorentzVector *fPtr2;
  TLorentzVector *fPtr3;
  TLorentzVector *fPtr4;
  TLorentzVector *fDiJet;

  TLorentzVector *fFPtr1;
  TLorentzVector *fFPtr2;

  TJet fJet1;
  TJet fJet2;
  TJet fJet3;
  TJet fJet4;
  TJet fFJet1;
  TJet fFJet2;

  unsigned int fHLTMatch;
  unsigned int fNoiseClean;
};
