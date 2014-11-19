#include "TClonesArray.h"
#include "TBranch.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector2.h"
#include "../include/Tools.hh"
#include "../include/MetDM.hh"

MetDM::MetDM(TTree *iTree) { 
  setupMetTree(iTree);
}
MetDM::~MetDM() { 
}
void MetDM::setupMetTree(TTree *iTree) { 
  iTree->Branch("met"      ,&fMet     ,"fMet/F");
  iTree->Branch("metphi"   ,&fMetPhi  ,"fMetPhil/F");
  iTree->Branch("metcov00" ,&fMetCov00,"fMetCov00/F");
  iTree->Branch("metcov01" ,&fMetCov01,"fMetCov01/F");
  iTree->Branch("metcov10" ,&fMetCov10,"fMetCov10/F");
  iTree->Branch("metcov11" ,&fMetCov11,"fMetCov1/F");
  iTree->Branch("mt1"      ,&fMt1     ,"fMt1/F");
  //iTree->Branch("m_r"      ,&fMR      ,"fMR/F");    // Mass using a kinematic liklhood fit

  fCorrector = new RecoilCorrector("files/recoilfit_zmm53XRR_2012_njet.root");
  fCorrector->fFactor = 1.3;
}
void MetDM::clearMet(){
  fMet      = 0;
  fMetPhi   = 0;
  fMetCov00 = 0;
  fMetCov01 = 0;
  fMetCov10 = 0;
  fMetCov11 = 0;
  fMt1      = 0;
  fMt2      = 0;
  //fPtH      = 0;
}
void MetDM::subtractGenJets(TLorentzVector &iRecoil,Jets &iJets,bool iSubtract) { 
  TLorentzVector lJet1; if(iJets.fGJetPt1 > 0) lJet1.SetPtEtaPhiM(iJets.fGJetPt1,iJets.fGJetEta1,iJets.fGJetPhi1,iJets.fGJetM1);
  TLorentzVector lJet2; if(iJets.fGJetPt2 > 0) lJet2.SetPtEtaPhiM(iJets.fGJetPt2,iJets.fGJetEta2,iJets.fGJetPhi2,iJets.fGJetM2);
  TLorentzVector lJet3; if(iJets.fGJCPt   > 0) lJet3.SetPtEtaPhiM(iJets.fGJCPt  ,iJets.fGJCEta  ,iJets.fGJCPhi  ,iJets.fGJCM);
  if(iSubtract) { 
    if(iJets.fGJetPt1 > 0) iRecoil -= lJet1;
    //if(iJets.fGJetPt2 > 0) iRecoil -= lJet2;
    //if(iJets.fGJCPt   > 0) iRecoil -= lJet3;
  } else { 
    if(iJets.fGJetPt1 > 0) iRecoil -= lJet1;
    //if(iJets.fGJetPt2 > 0) iRecoil -= lJet2;
    //if(iJets.fGJCPt   > 0) iRecoil -= lJet3;
  }
}
void MetDM::subtractJets(TLorentzVector &iRecoil,Jets &iJets,bool iSubtract) { 
  TLorentzVector lJet1; if(iJets.fJetPt1 > 0) lJet1.SetPtEtaPhiM(iJets.fJetPt1,iJets.fJetEta1,iJets.fJetPhi1,iJets.fJetM1);
  TLorentzVector lJet2; if(iJets.fJetPt2 > 0) lJet2.SetPtEtaPhiM(iJets.fJetPt2,iJets.fJetEta2,iJets.fJetPhi2,iJets.fJetM2);
  TLorentzVector lJet3; if(iJets.fJCPt   > 0) lJet3.SetPtEtaPhiM(iJets.fJCPt  ,iJets.fJCEta  ,iJets.fJCPhi  ,iJets.fJCM);
  if(iSubtract) { 
    if(iJets.fJetPt1 > 0) iRecoil -= lJet1;
    //if(iJets.fJetPt2 > 0) iRecoil -= lJet2;
    //if(iJets.fJCPt   > 0) iRecoil -= lJet3;
  } else { 
    if(iJets.fJetPt1 > 0) iRecoil -= lJet1;
    //if(iJets.fJetPt2 > 0) iRecoil -= lJet2;
    //if(iJets.fJCPt   > 0) iRecoil -= lJet3;
  }
}
void MetDM::fillMet(GenJet &iGen,Jets &iJets,bool iFake) { //Assume Fake is part of the recoil
  TLorentzVector lRecoil; lRecoil.SetPtEtaPhiM(iGen.fV1Pt,0,iGen.fV1Phi,0); lRecoil.RotateZ(TMath::Pi());
  subtractGenJets(lRecoil,iJets);
  TLorentzVector lMet;   lMet  .SetPtEtaPhiM(0   ,0,    0,0);
  lMet -= lRecoil;
  
  double pU1,pU2  = 0;
  double lVMet    = lMet.Pt(); 
  double lVMetPhi = lMet.Phi();
  int    lNJets   = 0;//=> Just assume we correct unclustered energy
  if(lRecoil.Pt() > 50) lNJets = 1.;
  lRecoil.RotateZ(TMath::Pi()); //=> Convert it to modified boson pT 
  
  //Smear the MET from code in my thesis using official CMS Resolutions
  //fCorrector->CorrectAll(lVMet, lVMetPhi, fVPt, fVPhi, (lLep1+lLep2).Pt(), (lLep1+lLep2).Phi(), pU1, pU2, 0, 0, lNJets);
  fCorrector->CorrectAll(lVMet, lVMetPhi, lRecoil.Pt(), lRecoil.Phi(), 0., 0., pU1, pU2, 0, 0, lNJets);
  lMet.SetPtEtaPhiM(lVMet,0,lVMetPhi,0);
  subtractJets(lMet,iJets);
  
  TLorentzVector lJet1; if(iJets.fJetPt1 > 0) lJet1.SetPtEtaPhiM(iJets.fJetPt1,iJets.fJetEta1,iJets.fJetPhi1,iJets.fJetM1);
  TLorentzVector lJet2; if(iJets.fJetPt2 > 0) lJet2.SetPtEtaPhiM(iJets.fJetPt2,iJets.fJetEta2,iJets.fJetPhi2,iJets.fJetM2);
  fMet    = lMet.Pt();
  fMetPhi = lMet.Phi();
  fMt1    = sqrt(2.0*(lJet1.Pt()*lMet.Pt()*(1.0-cos(Tools::deltaPhi(lJet1.Phi(),lMet.Phi())))));
  fMt2    = sqrt(2.0*(lJet2.Pt()*lMet.Pt()*(1.0-cos(Tools::deltaPhi(lJet2.Phi(),lMet.Phi())))));
  //fPtH    = (lMet + lLep1 + lLep2).Pt();
  //Get MET Covariance matrix in a lazy manner => b/c I have it saved here
  lNJets = 1;
  fCorrector->CorrectAll(lVMet, lVMetPhi, iGen.fVPt, iGen.fVPhi,0.,0., pU1, pU2, 0, 0, lNJets);
  fMetCov00 = (*(fCorrector->fCov))(0,0);
  fMetCov10 = (*(fCorrector->fCov))(1,0);
  fMetCov01 = (*(fCorrector->fCov))(0,1);
  fMetCov11 = (*(fCorrector->fCov))(1,1);
  //fMetCov00 = 400.; fMetCov11 = 400.; fMetCov10 = 20.; fMetCov01 = 20.;
}
