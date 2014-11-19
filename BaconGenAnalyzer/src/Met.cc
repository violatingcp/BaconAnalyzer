#include "TClonesArray.h"
#include "TBranch.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector2.h"
#include "../include/Tools.hh"
#include "../include/Met.hh"

Met::Met(TTree *iTree) { 
  setupMetTree(iTree);
}
Met::~Met() { 
}
void Met::setupMetTree(TTree *iTree) { 
  iTree->Branch("met"      ,&fMet     ,"fMet/F");
  iTree->Branch("metphi"   ,&fMetPhi  ,"fMetPhil/F");
  iTree->Branch("metcov00" ,&fMetCov00,"fMetCov00/F");
  iTree->Branch("metcov01" ,&fMetCov01,"fMetCov01/F");
  iTree->Branch("metcov10" ,&fMetCov10,"fMetCov10/F");
  iTree->Branch("metcov11" ,&fMetCov11,"fMetCov1/F");
  iTree->Branch("mt1"      ,&fMt1     ,"fMt1/F");
  iTree->Branch("mt2"      ,&fMt2     ,"fMt2/F");
  iTree->Branch("pzeta"    ,&fPZeta    ,"fPZeta/F");
  iTree->Branch("pzetavis" ,&fPZetaVis ,"fPZetaVis/F");
  iTree->Branch("pzetamiss",&fPZetaMiss,"fPZetaMiss/F");
  iTree->Branch("pth"      ,&fPtH      ,"fPtH/F");
  
  //System variables
  iTree->Branch("m_nt"     ,&fMNT     ,"fMNT/F");    // Mass using a kinematic liklhood fit
  iTree->Branch("pt_nt"    ,&fPtNT    ,"fPtNT/F");   // Pt using a kinematic likelihood fit
  iTree->Branch("phi_nt"   ,&fPhiNT   ,"fPhiNT/F");  // Phi using a kinematic likelihood fit
  iTree->Branch("y_nt"     ,&fYNT     ,"fYNT/F");    // Rapdity
  iTree->Branch("m_vis"    ,&fMVis    ,"fMVis/F");   // visible mass
  iTree->Branch("phi_vis"  ,&fPhiVis  ,"fPhiVis/F"); // phlYi vis 
  
  iTree->Branch("m_r"      ,&fMR      ,"fMR/F");    // Mass using a kinematic liklhood fit


  fCorrector = new RecoilCorrector("files/recoilfit_zmm53XRR_2012_njet.root");
  //fNTFitWZ   = new NTFit(-1,"files/params.root",NTFit::kME,NTFit::kMean,1000);
  //fNTFitWW   = new NTFit( 0,"files/params.root",NTFit::kME,NTFit::kMean,1000);

  fNTFitWZ   = new NTFit(-1,"files/params.root",NTFit::kME,NTFit::kMean,10000,true);
  fNTFitWW   = new NTFit( 0,"files/params.root",NTFit::kME,NTFit::kMean,10000,true);

  fCorrector->fFactor = 1.3;
}
void Met::clearMet(){
  fMet      = 0;
  fMetPhi   = 0;
  fMetCov00 = 0;
  fMetCov01 = 0;
  fMetCov10 = 0;
  fMetCov11 = 0;
  fMt1      = 0;
  fMt2      = 0;
  fPtH      = 0;
  
  //System variables
  fMNT    = 0;
  fPtNT   = 0;
  fPhiNT  = 0;
  fYNT    = 0;
  
  fMVis   = 0;
  fPhiVis = 0;
  
  fPZetaMiss = 0;
  fPZetaVis  = 0;
  fPZeta     = 0;
}
void Met::subtractGenJets(TLorentzVector &iRecoil,Jets &iJets,bool iSubtract) { 
  TLorentzVector lJet1; if(iJets.fGJetPt1 > 0) lJet1.SetPtEtaPhiM(iJets.fGJetPt1,iJets.fGJetEta1,iJets.fGJetPhi1,iJets.fGJetM1);
  TLorentzVector lJet2; if(iJets.fGJetPt2 > 0) lJet2.SetPtEtaPhiM(iJets.fGJetPt2,iJets.fGJetEta2,iJets.fGJetPhi2,iJets.fGJetM2);
  TLorentzVector lJet3; if(iJets.fGJCPt   > 0) lJet3.SetPtEtaPhiM(iJets.fGJCPt  ,iJets.fGJCEta  ,iJets.fGJCPhi  ,iJets.fGJCM);
  if(iSubtract) { 
    if(iJets.fGJetPt1 > 0) iRecoil -= lJet1;
    if(iJets.fGJetPt2 > 0) iRecoil -= lJet2;
    if(iJets.fGJCPt   > 0) iRecoil -= lJet3;
  } else { 
    if(iJets.fGJetPt1 > 0) iRecoil -= lJet1;
    if(iJets.fGJetPt2 > 0) iRecoil -= lJet2;
    if(iJets.fGJCPt   > 0) iRecoil -= lJet3;
  }
}
void Met::subtractJets(TLorentzVector &iRecoil,Jets &iJets,bool iSubtract) { 
  TLorentzVector lJet1; if(iJets.fJetPt1 > 0) lJet1.SetPtEtaPhiM(iJets.fJetPt1,iJets.fJetEta1,iJets.fJetPhi1,iJets.fJetM1);
  TLorentzVector lJet2; if(iJets.fJetPt2 > 0) lJet2.SetPtEtaPhiM(iJets.fJetPt2,iJets.fJetEta2,iJets.fJetPhi2,iJets.fJetM2);
  TLorentzVector lJet3; if(iJets.fJCPt   > 0) lJet3.SetPtEtaPhiM(iJets.fJCPt  ,iJets.fJCEta  ,iJets.fJCPhi  ,iJets.fJCM);
  if(iSubtract) { 
    if(iJets.fJetPt1 > 0) iRecoil -= lJet1;
    if(iJets.fJetPt2 > 0) iRecoil -= lJet2;
    if(iJets.fJCPt   > 0) iRecoil -= lJet3;
  } else { 
    if(iJets.fJetPt1 > 0) iRecoil -= lJet1;
    if(iJets.fJetPt2 > 0) iRecoil -= lJet2;
    if(iJets.fJCPt   > 0) iRecoil -= lJet3;
  }
}
void Met::fillMet(Gen &iGen,Jets &iJets,bool iFake) { //Assume Fake is part of the recoil
  TLorentzVector lRecoil; lRecoil.SetPtEtaPhiM(iGen.fVPt,0,iGen.fVPhi,0); lRecoil.RotateZ(TMath::Pi());
  subtractGenJets(lRecoil,iJets);
  TLorentzVector lLep1;  lLep1 .SetPtEtaPhiM(iGen.fPt1,iGen.fEta1,iGen.fPhi1,0.1);//iGen.fM1);
  TLorentzVector lLep2;  lLep2 .SetPtEtaPhiM(iGen.fPt2,iGen.fEta2,iGen.fPhi2,0.1);//iGen.fM2);
  for(int i0 = 2; i0 < int(iJets.fLeptons.size()); i0++) lLep2 += iJets.fLeptons[i0]; 
  TLorentzVector lMet;   lMet  .SetPtEtaPhiM(0   ,0,    0,0);
  lMet -= lRecoil;
  lMet -= lLep1;
  if(!iFake) lMet -= lLep2;
  
  double pU1,pU2  = 0;
  double lVMet    = lMet.Pt(); 
  double lVMetPhi = lMet.Phi();
  int    lNJets   = 0;//=> Just assume we correct unclustered energy
  lRecoil.RotateZ(TMath::Pi()); //=> Convert it to modified boson pT 
  //Smear the MET from code in my thesis using official CMS Resolutions
  //fCorrector->CorrectAll(lVMet, lVMetPhi, fVPt, fVPhi, (lLep1+lLep2).Pt(), (lLep1+lLep2).Phi(), pU1, pU2, 0, 0, lNJets);
  fCorrector->CorrectAll(lVMet, lVMetPhi, lRecoil.Pt(), lRecoil.Phi(), (lLep1+lLep2).Pt(), (lLep1+lLep2).Phi(), pU1, pU2, 0, 0, lNJets);
  //cout << "==> " << lMet.Pt() << " -- " << lVMet << " -- " << lRecoil.Pt() << " -- " << fVPt  << endl;
  lMet.SetPtEtaPhiM(lVMet,0,lVMetPhi,0);
  subtractJets(lMet,iJets);
  
  fMet    = lMet.Pt();
  fMetPhi = lMet.Phi();
  fMt1    = sqrt(2.0*(lLep1.Pt()*lMet.Pt()*(1.0-cos(Tools::deltaPhi(lLep1.Phi(),lMet.Phi())))));
  fMt2    = sqrt(2.0*(lLep2.Pt()*lMet.Pt()*(1.0-cos(Tools::deltaPhi(lLep2.Phi(),lMet.Phi())))));
  fPtH    = (lMet + lLep1 + lLep2).Pt();
  //Get MET Covariance matrix in a lazy manner => b/c I have it saved here
  lNJets = 2;
  fCorrector->CorrectAll(lVMet, lVMetPhi, iGen.fVPt, iGen.fVPhi, (lLep1+lLep2).Pt(), (lLep1+lLep2).Phi(), pU1, pU2, 0, 0, lNJets);
  fMetCov00 = (*(fCorrector->fCov))(0,0);
  fMetCov10 = (*(fCorrector->fCov))(1,0);
  fMetCov01 = (*(fCorrector->fCov))(0,1);
  fMetCov11 = (*(fCorrector->fCov))(1,1);
  //fMetCov00 = 400.; fMetCov11 = 400.; fMetCov10 = 20.; fMetCov01 = 20.;
  TLorentzVector lBoson = lLep1 + lLep2;

  NTFit *lNTFit = 0;
  if(iGen.fPt4   == 0) lNTFit = fNTFitWZ;
  if(iGen.fPt3   == 0) lNTFit = fNTFitWW;
  lNTFit = fNTFitWW;
  //Compute the Mass from the likelihood 
  if(std::isnan(lLep2.Pt())) return;
  if(iJets.fLeptons.size() < 4) {  
    lBoson = lNTFit->NeutrinoToys(lNJets,lNJets,
				  -1,lLep1.Pt(),lLep1.Phi(),lLep1.Eta(),lLep1.M(),
				  -1,lLep2.Pt(),lLep2.Phi(),lLep2.Eta(),lLep2.M(),
				  fMet,fMetPhi,fMetCov00,fMetCov10,fMetCov10,fMetCov11);
  } 
  TLorentzVector lVec; lVec.SetPtEtaPhiM(fMet,0,fMetPhi,0);
  //std::cout << "Mass ===> " << lBoson.M() << " -- " << (lLep1+lLep2+lVec).M() << " -- " << iGen.fVM << " -- " << iJets.fLeptons.size() << std::endl;
  //lBoson   += iJets.fLeptons[2];
  fMNT      = lBoson.M();
  fPtNT     = lBoson.Pt();
  fPhiNT    = lBoson.Phi();
  fYNT      = lBoson.Rapidity();

  fMVis     = (lLep1 + lLep2).M();
  fPhiVis   = (lLep1 + lLep2).Phi();
  double lPtVis = (lLep1 + lLep2).Pt();
  
  TVector3 lBisector(lLep1.Vect().Unit() + lLep2.Vect().Unit());
  fPZetaVis      = (lLep1+lLep2).Vect().Dot(lBisector);
  fPZetaMiss     = lMet.Vect().Dot(lBisector);
  fPZeta         = 0.85*fPZetaVis-fPZetaMiss;
  
  TVector2 lRVMet; lRVMet.SetMagPhi(fMet,fMetPhi);
  TVector2 lRVis;  lRVis .SetMagPhi(lPtVis,fPhiVis);
  fMR       = sqrt(0.5*(fMVis*fMVis-(lRVMet*lRVis)+
			sqrt((fMVis*fMVis+lPtVis*lPtVis)*(fMVis*fMVis+fMet*fMet)) ));
}
